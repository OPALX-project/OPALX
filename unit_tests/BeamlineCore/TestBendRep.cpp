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
#include "BeamlineCore/MultipoleRep.h"
#include "BeamlineGeometry/Geometry.h"
#include "PartBunch/BunchStateHandler.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"
#include "gtest/gtest.h"

#include <cmath>

extern Inform* gmsg;

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
            gmsg = new Inform(nullptr, -1);
        }
        static void TearDownTestSuite() {
            delete gmsg;
            gmsg = nullptr;
            ippl::finalize();
        }
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
    bend.setFullGap(2.0 * halfGap);
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
    EXPECT_NEAR(fieldBegin, -10.0 * halfGap, 1.0e-12);  // -5 * (2*halfGap)
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
    const Vector3 exitFace(
            (std::cos(phiExit) - 1.0) / curvature, 0.0, std::sin(phiExit) / curvature);
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
    bend.setFullGap(2.0 * halfGap);
    bend.setFringeIntegral(0.5);
    bend.setB(-1.0);

    const double profileGap = 2.0 * halfGap;
    const double halfWidth  = 5.0 * profileGap;

    double fieldBegin = 0.0;
    double fieldEnd   = 0.0;
    bend.getFieldExtent(fieldBegin, fieldEnd);
    EXPECT_NEAR(fieldBegin, -halfWidth, 1.0e-12);            // entrance face is perpendicular
    EXPECT_NEAR(fieldEnd, bodyLength + halfWidth, 1.0e-12);  // box length, exit perpendicular

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
    bend.setFullGap(2.0 * halfGap);
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
    straight.setFullGap(2.0 * halfGap);
    straight.setFringeIntegral(0.0);
    straight.setB(-1.0);

    B = Vector3(0.0);
    straight.applyToReferenceParticle(Vector3(0.0, yOffset, -0.05), Vector3(0.0), 0.0, E, B);
    EXPECT_NEAR(B(0), 0.0, 1.0e-12);
}

// A particle exactly on a shared face must receive the downstream field only.
// Exercise the production device apply(), both host field APIs and spatial
// membership at identical points. Far transverse particles inside the same
// longitudinal slab must not receive a remote ring magnet's field. Existing
// gathered fields must remain additive; this field query must not mark losses.
TEST_F(BendRepTest, DeviceAndHostHardEdgeFieldsUseIdenticalSpatialSupport) {
    ippl::NDIndex<3> domain;
    for (unsigned d = 0; d < 3; ++d) domain[d] = ippl::Index(8);
    ippl::UniformCartesian<double, 3> mesh(domain, Vector3(0.25), Vector3(-1));
    std::array<bool, 3> decomp{true, true, true};
    ippl::FieldLayout<3> layout(MPI_COMM_WORLD, domain, decomp, false);
    auto pc = std::make_shared<ParticleContainer_t>(mesh, layout);
    pc->setBunchStateHandler(std::make_shared<BunchStateHandler>());
    constexpr unsigned count = 9;
    pc->allocateParticles(count);
    pc->createParticles(count);

    const double sampleS[count] = {-1e-9, 0., 1e-9, 0.5, 1. - 1e-9, 1., 1. + 1e-9, 0.5, 0.5};
    const bool inside[count] = {false, true, true, true, true, false, false, false, false};
    const Vector3 initialE(0.1, 0.2, 0.3), initialB(0.4, 0.5, 0.6);
    for (unsigned kind = 0; kind < 3; ++kind) {
        SBendRep sector("sector");
        RBendRep rectangular("rectangular");
        MultipoleRep multipole("multipole");
        constexpr double curvature = 0.2;
        sector.getGeometry() = Geometry::makeSBend(1., curvature);
        sector.getGeometry().setElementLength(1.);
        sector.getGeometry().setBendAngle(curvature);
        sector.setB(-1.);
        rectangular.getGeometry() = Geometry::makeRBend(1., 0.2);
        rectangular.getGeometry().setElementLength(1.);
        rectangular.setB(-1.);
        multipole.getGeometry().setElementLength(1.);
        // The normal dipole setter stores half its argument.
        multipole.setNormalComponent(0, -2.);
        ElementBase& element = kind == 0 ? static_cast<ElementBase&>(sector)
                               : kind == 1 ? static_cast<ElementBase&>(rectangular)
                                           : static_cast<ElementBase&>(multipole);
        SCOPED_TRACE(element.getName());
        element.setAperture(ApertureType::RECTANGULAR, {1., 1.});
        auto r = Kokkos::create_mirror_view(pc->R.getView());
        for (unsigned i = 0; i < count; ++i) {
            const double s = sampleS[i];
            r(i) = Vector3(0, 0, s);
            if (kind == 0 && s > 0) {
                const double phi = curvature * std::min(s, 1.);
                r(i) = Vector3((std::cos(phi) - 1.) / curvature, 0., std::sin(phi) / curvature);
                if (s > 1.) r(i) += (s - 1.) * Vector3(-std::sin(curvature), 0, std::cos(curvature));
            }
            if (i == 7) {
                // Shift in the bend's radial direction while preserving arc s.
                const double phi = kind == 0 ? curvature * s : 0.;
                r(i) += 2. * Vector3(std::cos(phi), 0, std::sin(phi));
            }
            if (i == 8) r(i)(1) = 2.;
        }
        Kokkos::deep_copy(pc->R.getView(), r);
        Kokkos::deep_copy(pc->E.getView(), initialE);
        Kokkos::deep_copy(pc->B.getView(), initialB);
        element.apply(pc);
        const auto e = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->E.getView());
        const auto b = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->B.getView());
        const auto invalid = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->InvalidMask.getView());
        for (unsigned i = 0; i < count; ++i) {
            SCOPED_TRACE(i);
            EXPECT_EQ(element.isInside(r(i)), inside[i]);
            Vector3 hostE(initialE), hostB(initialB), refE(initialE), refB(initialB);
            element.apply(r(i), Vector3(0, 0, 1), 0., hostE, hostB);
            EXPECT_EQ(element.applyToReferenceParticle(r(i), Vector3(0, 0, 1), 0., refE, refB), i >= 7);
            EXPECT_FALSE(invalid(i));
            const Vector3 expectedB = initialB + Vector3(0, inside[i] ? -1. : 0., 0);
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_DOUBLE_EQ(e(i)(d), initialE(d));
                EXPECT_NEAR(b(i)(d), expectedB(d), 1e-14);
                EXPECT_NEAR(hostB(d), expectedB(d), 1e-14);
                EXPECT_NEAR(refB(d), expectedB(d), 1e-14);
                EXPECT_DOUBLE_EQ(hostE(d), initialE(d));
                EXPECT_DOUBLE_EQ(refE(d), initialE(d));
            }
        }
    }
}
