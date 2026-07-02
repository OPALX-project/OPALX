//
// Tests for the analytic SBEND/RBEND fringe field (Enge profile + pole-face edge
// focusing) evaluated through SBendRep/RBendRep.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "BeamlineCore/RBendRep.h"
#include "BeamlineCore/SBendRep.h"
#include "BeamlineGeometry/Geometry.h"
#include "Physics/Physics.h"
#include "gtest/gtest.h"

#include <cmath>

namespace {
    using Vector3 = Vector_t<double, 3>;

    // F(0) at the nominal pole face for the OPAL default Enge profile.
    const double edgeScale = 1.0 / (1.0 + std::exp(0.478959));

    class BendRepTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }
        static void TearDownTestSuite() { ippl::finalize(); }
    };
}  // namespace

// The SBEND field support extends one Enge fringe half width (5 gaps) past each
// pole face, and the on-axis dipole follows the OPAL Enge profile as a function of
// the arc length s along the curved reference orbit: F(0)=0.3825 at the entrance
// face, F=1 in the body, F=0.3825 at the exit face, 0 outside. Because the body is
// curved, the sample points are taken on the design arc (framePosition), not at a
// rigid Cartesian z (which for a 45 deg bend is not the same as s).
TEST_F(BendRepTest, SBendFringeSupportProducesExpectedOpalEngeProfile) {
    constexpr double bodyLength = 1.0;
    constexpr double angle      = Physics::pi / 4.0;
    constexpr double halfGap    = 0.02;
    const double curvature      = angle / bodyLength;

    SBendRep bend("SBEND");
    bend.getGeometry() = Geometry::makeSBend(bodyLength, curvature);
    bend.getGeometry().setElementLength(bodyLength);
    bend.getGeometry().setBendAngle(angle);
    bend.setFringeHalfGap(halfGap);
    bend.setFringeIntegral(0.5);
    bend.setB(-1.0);

    // Local Cartesian point on the reference orbit at arc length s: the circular arc
    // inside the body, the straight entrance/exit tangents outside it. This is the
    // forward map inverted by SBend::bendCoords.
    auto pointAtArc = [&](double s) -> Vector3 {
        if (s <= 0.0) {
            return Vector3(0.0, 0.0, s);  // straight entrance tangent (+z)
        }
        const double sBody = std::min(s, bodyLength);
        const double phi   = curvature * sBody;
        Vector3 p((std::cos(phi) - 1.0) / curvature, 0.0, std::sin(phi) / curvature);
        if (s > bodyLength) {  // straight exit tangent past the exit face
            const Vector3 tangent(-std::sin(angle), 0.0, std::cos(angle));
            p += (s - bodyLength) * tangent;
        }
        return p;
    };

    double fieldBegin = 0.0;
    double fieldEnd   = 0.0;
    bend.getFieldExtent(fieldBegin, fieldEnd);
    EXPECT_NEAR(fieldBegin, -10.0 * halfGap, 1.0e-12);              // -5 * (2*halfGap)
    EXPECT_NEAR(fieldEnd, bodyLength + 10.0 * halfGap, 1.0e-12);

    Vector3 E(0.0);
    Vector3 B(0.0);

    // Just outside the entrance support: not selected, no field.
    EXPECT_FALSE(bend.applyToReferenceParticle(pointAtArc(-0.21), Vector3(0.0), 0.0, E, B));
    EXPECT_NEAR(B(1), 0.0, 1.0e-12);

    B = Vector3(0.0);
    bend.applyToReferenceParticle(pointAtArc(0.0), Vector3(0.0), 0.0, E, B);
    EXPECT_NEAR(B(1), -edgeScale, 1.0e-12);  // entrance face

    B = Vector3(0.0);
    bend.applyToReferenceParticle(pointAtArc(0.5), Vector3(0.0), 0.0, E, B);
    EXPECT_NEAR(B(1), -1.0, 1.0e-12);  // body interior

    B = Vector3(0.0);
    bend.applyToReferenceParticle(pointAtArc(bodyLength), Vector3(0.0), 0.0, E, B);
    EXPECT_NEAR(B(1), -edgeScale, 1.0e-12);  // exit face

    B = Vector3(0.0);
    EXPECT_FALSE(bend.applyToReferenceParticle(pointAtArc(1.21), Vector3(0.0), 0.0, E, B));
    EXPECT_NEAR(B(1), 0.0, 1.0e-12);

    // The longitudinal field Bz = B0 F'(s) y at the entrance face (s = 0, where the
    // arc-tangent basis coincides with the entrance frame).
    const double expAtEdge  = std::exp(0.478959);
    const double profileGap = 2.0 * halfGap;
    const double dFdz       = 1.911289 * expAtEdge / (profileGap * std::pow(1.0 + expAtEdge, 2.0));
    const double yOffset    = 1.0e-3;
    B                       = Vector3(0.0);
    Vector3 entryWithY      = pointAtArc(0.0);
    entryWithY(1)           = yOffset;
    bend.applyToReferenceParticle(entryWithY, Vector3(0.0), 0.0, E, B);
    EXPECT_NEAR(B(2), -dFdz * yOffset, 1.0e-12);
}

// With no gap the fringe collapses to a hard edge: extent is the plain body, the
// dipole is full inside and absent outside. The exit edge is sampled at its true arc
// position on the curved reference orbit (arc length s = bodyLength), not at z = L.
TEST_F(BendRepTest, SBendWithoutGapIsHardEdge) {
    constexpr double bodyLength = 1.0;
    constexpr double curvature  = 0.2;

    SBendRep bend("SBEND");
    bend.getGeometry() = Geometry::makeSBend(bodyLength, curvature);
    bend.getGeometry().setElementLength(bodyLength);
    bend.getGeometry().setBendAngle(curvature * bodyLength);
    bend.setB(-1.0);

    double fieldBegin = 0.0;
    double fieldEnd   = 0.0;
    bend.getFieldExtent(fieldBegin, fieldEnd);
    EXPECT_NEAR(fieldBegin, 0.0, 1.0e-15);
    EXPECT_NEAR(fieldEnd, bodyLength, 1.0e-15);

    // Body interior at arc length 0.5.
    const double phiMid = curvature * 0.5;
    Vector3 E(0.0);
    Vector3 B(0.0);
    bend.applyToReferenceParticle(
            Vector3((std::cos(phiMid) - 1.0) / curvature, 0.0, std::sin(phiMid) / curvature),
            Vector3(0.0), 0.0, E, B);
    EXPECT_NEAR(B(1), -1.0, 1.0e-15);  // full dipole, no fringe scaling

    // Exit edge (arc length = bodyLength) is exclusive: not selected, no field.
    const double phiExit = curvature * bodyLength;
    const Vector3 exitFace((std::cos(phiExit) - 1.0) / curvature, 0.0, std::sin(phiExit) / curvature);
    B = Vector3(0.0);
    EXPECT_FALSE(bend.applyToReferenceParticle(exitFace, Vector3(0.0), 0.0, E, B));
    EXPECT_NEAR(B(1), 0.0, 1.0e-15);
}

// The RBEND is evaluated in its straight box frame (+z along the box axis), so the
// field support is the plain box length plus one Enge fringe half width past each face.
// The faces are perpendicular to the box axis (tilted only by an explicit E1/E2), so with
// E1=E2=0 there is no projection: the extent is [-halfWidth, L + halfWidth]. The intrinsic
// angle/2 is an orbit-to-face angle used for edge focusing, not a face tilt in this frame.
TEST_F(BendRepTest, RBendFringeSupportProjectsExitByFaceAngle) {
    constexpr double bodyLength = 1.0;
    constexpr double angle      = 0.2;
    constexpr double halfGap    = 0.02;

    RBendRep bend("RBEND");
    bend.getGeometry() = Geometry::makeRBend(bodyLength, angle);
    bend.getGeometry().setElementLength(bodyLength);
    bend.getGeometry().setBendAngle(angle);
    bend.setFringeHalfGap(halfGap);
    bend.setFringeIntegral(0.5);
    bend.setB(-1.0);

    const double profileGap = 2.0 * halfGap;
    const double halfWidth   = 5.0 * profileGap;

    double fieldBegin = 0.0;
    double fieldEnd   = 0.0;
    bend.getFieldExtent(fieldBegin, fieldEnd);
    EXPECT_NEAR(fieldBegin, -halfWidth, 1.0e-12);              // entrance face is perpendicular
    EXPECT_NEAR(fieldEnd, bodyLength + halfWidth, 1.0e-12);    // box length, exit perpendicular

    Vector3 E(0.0);
    Vector3 B(0.0);
    bend.applyToReferenceParticle(Vector3(0.0, 0.0, 0.0), Vector3(0.0), 0.0, E, B);
    EXPECT_NEAR(B(1), -edgeScale, 1.0e-12);  // entrance face (box z = 0)

    B = Vector3(0.0);
    bend.applyToReferenceParticle(Vector3(0.0, 0.0, 0.5), Vector3(0.0), 0.0, E, B);
    EXPECT_NEAR(B(1), -1.0, 1.0e-12);  // body interior
}

// Inside the entrance fringe the vertical edge focusing appears as a horizontal
// field Bx proportional to y, and vanishes when FINT = 0 with a perpendicular face.
TEST_F(BendRepTest, SBendEntryFringeAddsHorizontalEdgeField) {
    constexpr double bodyLength = 1.0;
    constexpr double angle      = 0.2;
    constexpr double halfGap    = 0.02;
    constexpr double entryAngle = 0.15;
    constexpr double yOffset    = 1.0e-3;

    SBendRep bend("SBEND");
    bend.getGeometry() = Geometry::makeSBend(bodyLength, angle / bodyLength);
    bend.getGeometry().setElementLength(bodyLength);
    bend.getGeometry().setBendAngle(angle);
    bend.getGeometry().setEntranceAngle(entryAngle);
    bend.getGeometry().setExitAngle(0.0);
    bend.setFringeHalfGap(halfGap);
    bend.setFringeIntegral(0.5);
    bend.setB(-1.0);

    // A point inside the entrance fringe (upstream of the body) with a vertical offset.
    Vector3 E(0.0);
    Vector3 B(0.0);
    bend.applyToReferenceParticle(Vector3(0.0, yOffset, -0.05), Vector3(0.0), 0.0, E, B);
    EXPECT_GT(std::abs(B(0)), 0.0);  // edge focusing present

    // With FINT = 0 and a perpendicular entrance face the edge angle is zero, so no Bx.
    SBendRep straight("SBEND");
    straight.getGeometry() = Geometry::makeSBend(bodyLength, angle / bodyLength);
    straight.getGeometry().setElementLength(bodyLength);
    straight.getGeometry().setBendAngle(angle);
    straight.getGeometry().setEntranceAngle(0.0);
    straight.getGeometry().setExitAngle(0.0);
    straight.setFringeHalfGap(halfGap);
    straight.setFringeIntegral(0.0);
    straight.setB(-1.0);

    B = Vector3(0.0);
    straight.applyToReferenceParticle(Vector3(0.0, yOffset, -0.05), Vector3(0.0), 0.0, E, B);
    EXPECT_NEAR(B(0), 0.0, 1.0e-12);
}
