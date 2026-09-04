//
// Tests for the combined-function multipole with a tanh fringe, straight and
// curved, evaluated through MultipoleTRep.
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

#include "BeamlineCore/MultipoleTRep.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "AbsBeamline/MultipoleTFieldModel.h"
#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/SplineTimeDependence.h"
#include "BeamlineGeometry/Geometry.h"
#include "Beamlines/Beamline.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

#include <cmath>
#include <memory>
#include <vector>

namespace {
    using Vector3 = Vector_t<double, 3>;

    class MultipoleTRepTest : public ::testing::Test {
    protected:
        // The element allocates a Kokkos view for its tanh table in the constructor.
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }
        static void TearDownTestSuite() { ippl::finalize(); }
    };

    /// A straight magnet of the given profile and fringe.
    MultipoleTRep makeStraight(
            const std::vector<double>& profile, const double length, const double fringe) {
        MultipoleTRep magnet("MULTIPOLET");
        magnet.getGeometry() = Geometry::makeStraight(length);
        magnet.setMaxFOrder(5);
        magnet.setTransverseProfile(profile);
        magnet.setFringeLengths(fringe, fringe);
        return magnet;
    }

    /// A curved magnet of the given profile, arc length, bend angle and fringe.
    MultipoleTRep makeCurved(
            const std::vector<double>& profile, const double length, const double angle,
            const double fringe) {
        MultipoleTRep magnet("MULTIPOLET");
        magnet.getGeometry() = Geometry::makeSBend(length, angle / length);
        magnet.setMaxFOrder(5);
        magnet.setTransverseProfile(profile);
        magnet.setFringeLengths(fringe, fringe);
        return magnet;
    }

    /// Local Cartesian point on the design orbit at arc length s: the arc inside the
    /// body, the straight entrance and exit tangents outside it. This is the forward
    /// map that MultipoleT::bendCoords inverts.
    Vector3 pointAtArc(const double s, const double length, const double angle) {
        if (angle == 0.0) {
            return Vector3(0.0, 0.0, s);
        }
        const double curvature = angle / length;
        if (s <= 0.0) {
            return Vector3(0.0, 0.0, s);
        }
        const double sBody = std::min(s, length);
        const double phi   = curvature * sBody;
        Vector3 p((std::cos(phi) - 1.0) / curvature, 0.0, std::sin(phi) / curvature);
        if (s > length) {
            const Vector3 tangent(-std::sin(angle), 0.0, std::cos(angle));
            p += (s - length) * tangent;
        }
        return p;
    }

    /// The field the element reports at a local point.
    Vector3 fieldAt(MultipoleTRep& magnet, const Vector3& R, const double t = 0.0) {
        Vector3 E(0.0);
        Vector3 B(0.0);
        magnet.apply(R, Vector3(0.0), t, E, B);
        return B;
    }
}  // namespace

// ---------------------------------------------------------------------------
// Field extent
// ---------------------------------------------------------------------------

// The support is the body plus the fringe reach on each side.
TEST_F(MultipoleTRepTest, FieldExtentCoversTheFringe) {
    const double length = 1.0;
    const double fringe = 0.02;
    MultipoleTRep magnet = makeStraight({1.0}, length, fringe);

    double zBegin = 0.0;
    double zEnd   = 0.0;
    magnet.getFieldExtent(zBegin, zEnd);
    EXPECT_NEAR(zBegin, -MultipoleTFieldModel::FringeReach * fringe, 1.0e-12);
    EXPECT_NEAR(zEnd, length + MultipoleTFieldModel::FringeReach * fringe, 1.0e-12);
}

// A zero fringe length is a hard edge, so the support is the plain body.
TEST_F(MultipoleTRepTest, FieldExtentOfAHardEdgeIsTheBody) {
    const double length = 1.0;
    MultipoleTRep magnet = makeStraight({1.0}, length, 0.0);

    double zBegin = 0.0;
    double zEnd   = 0.0;
    magnet.getFieldExtent(zBegin, zEnd);
    EXPECT_DOUBLE_EQ(zBegin, 0.0);
    EXPECT_DOUBLE_EQ(zEnd, length);
}

// ---------------------------------------------------------------------------
// The straight magnet
// ---------------------------------------------------------------------------

// The dipole is the flat-top field in the body, half of it at each face, and gone
// beyond the fringe reach.
TEST_F(MultipoleTRepTest, StraightDipoleAlongTheAxis) {
    const double length = 1.0;
    const double fringe = 0.02;
    MultipoleTRep magnet = makeStraight({-0.5}, length, fringe);

    EXPECT_NEAR(fieldAt(magnet, Vector3(0.0, 0.0, 0.5 * length))(1), -0.5, 1.0e-9);
    EXPECT_NEAR(fieldAt(magnet, Vector3(0.0, 0.0, 0.0))(1), -0.25, 1.0e-9);
    EXPECT_NEAR(fieldAt(magnet, Vector3(0.0, 0.0, length))(1), -0.25, 1.0e-9);
    // Outside the support the element is not selected at all.
    EXPECT_NEAR(fieldAt(magnet, Vector3(0.0, 0.0, -0.2))(1), 0.0, 1.0e-12);
    EXPECT_NEAR(fieldAt(magnet, Vector3(0.0, 0.0, length + 0.2))(1), 0.0, 1.0e-12);
}

// A hard-edge magnet has the flat top right up to the faces; the exit face is
// exclusive, as everywhere else in the tracker.
TEST_F(MultipoleTRepTest, StraightHardEdgeDipole) {
    const double length = 1.0;
    MultipoleTRep magnet = makeStraight({1.0}, length, 0.0);

    EXPECT_NEAR(fieldAt(magnet, Vector3(0.0, 0.0, 0.0))(1), 1.0, 1.0e-12);
    EXPECT_NEAR(fieldAt(magnet, Vector3(0.0, 0.0, 0.5 * length))(1), 1.0, 1.0e-12);
    EXPECT_NEAR(fieldAt(magnet, Vector3(0.0, 0.0, length))(1), 0.0, 1.0e-12);
}

// The gradient is linear in x with the usual sign: a positive TP[1] raises By at
// positive x.
TEST_F(MultipoleTRepTest, StraightGradientIsLinearInX) {
    const double length = 1.0;
    MultipoleTRep magnet = makeStraight({0.0, 2.0}, length, 0.02);

    for (const double x : {-0.05, -0.01, 0.0, 0.01, 0.05}) {
        EXPECT_NEAR(fieldAt(magnet, Vector3(x, 0.0, 0.5 * length))(1), 2.0 * x, 1.0e-9);
    }
}

// ---------------------------------------------------------------------------
// The curved magnet
// ---------------------------------------------------------------------------

// On the design orbit of a curved dipole the field stays the flat-top value and
// purely vertical all the way round, and the element stays selected: this is the
// rotation from the local Frenet basis into the entrance frame.
TEST_F(MultipoleTRepTest, CurvedDipoleFollowsTheDesignOrbit) {
    const double length = 1.0;
    const double angle  = Physics::pi / 6.0;
    MultipoleTRep magnet = makeCurved({-1.0}, length, angle, 0.02);

    for (const double fraction : {0.25, 0.5, 0.75}) {
        const double s      = fraction * length;
        const Vector3 point = pointAtArc(s, length, angle);
        EXPECT_TRUE(magnet.isInside(point)) << "at s = " << s;
        const Vector3 B = fieldAt(magnet, point);
        EXPECT_NEAR(B(1), -1.0, 1.0e-9) << "at s = " << s;
        EXPECT_NEAR(B(0), 0.0, 1.0e-12) << "at s = " << s;
        EXPECT_NEAR(B(2), 0.0, 1.0e-12) << "at s = " << s;
    }
}

// The centre of curvature is on the -x side, like SBEND: the design orbit at arc
// length s sits at x = (cos(h s) - 1) / h, i.e. at negative x.
TEST_F(MultipoleTRepTest, CurvedBodyBendsTowardsNegativeX) {
    const double length = 1.0;
    const double angle  = Physics::pi / 6.0;
    MultipoleTRep magnet = makeCurved({-1.0}, length, angle, 0.02);

    const Vector3 exit = pointAtArc(length, length, angle);
    EXPECT_LT(exit(0), 0.0);
    EXPECT_NEAR(exit(0), (std::cos(angle) - 1.0) * length / angle, 1.0e-12);
    // The geometry agrees with the field frame: same curvature, same bend angle.
    EXPECT_NEAR(magnet.getGeometry().getCurvature(), angle / length, 1.0e-12);
    EXPECT_NEAR(magnet.getGeometry().getBendAngle(), angle, 1.0e-12);
    EXPECT_TRUE(magnet.getGeometry().isBend());
}

// Off the mid-plane in the exit fringe the local field has a longitudinal part,
// and what the element reports is that local field rotated into the entrance
// frame by the arc angle. Getting this wrong is what gave a curved MULTIPOLET a
// spurious vertical kick.
TEST_F(MultipoleTRepTest, CurvedFieldIsRotatedIntoTheEntranceFrame) {
    const double length  = 1.0;
    const double angle   = Physics::pi / 3.0;
    const double fringe  = 0.05;
    const double yOffset = 1.0e-3;
    MultipoleTRep magnet = makeCurved({-1.0}, length, angle, fringe);

    const double s      = 0.98 * length;  // inside the exit fringe
    Vector3 point       = pointAtArc(s, length, angle);
    point(1)            = yOffset;
    const Vector3 B     = fieldAt(magnet, point);

    // Rebuild the expectation from the model plus the rotation.
    MultipoleTFieldModel::FieldInputs in{};
    in.profile[0]  = -1.0;
    in.maxFOrder   = magnet.getMaxFOrder();
    in.bodyLength  = length;
    in.curvature   = angle / length;
    in.fringeLeft  = fringe;
    in.fringeRight = fringe;

    Kokkos::View<double**, Kokkos::HostSpace> table(
            "table", 2 * in.maxFOrder + 2, 2 * in.maxFOrder + 3);
    MultipoleTFieldModel::tanhDerivativeTable(2 * in.maxFOrder + 1, table);

    const Vector3 arc   = GeometryHelper::toBendArcCoords(point, in.curvature, in.bodyLength);
    const Vector3 local = MultipoleTFieldModel::field(arc, in, table);
    const Vector3 expected =
            GeometryHelper::rotateArcFieldToEntry(local, arc(2), in.curvature, in.bodyLength);

    // The local field really does have a longitudinal part here, so the rotation matters.
    EXPECT_GT(std::abs(local(2)), 1.0e-9);
    for (unsigned int d = 0; d < 3; ++d) {
        EXPECT_NEAR(B(d), expected(d), 1.0e-12) << "component " << d;
    }
    // A pure rotation preserves the magnitude.
    EXPECT_NEAR(
            std::hypot(B(0), B(1), B(2)), std::hypot(local(0), local(1), local(2)), 1.0e-12);
}

// ---------------------------------------------------------------------------
// Containment and the reference particle
// ---------------------------------------------------------------------------

// Containment follows the arc, so a curved body stays selected as the orbit turns
// through it even though its straight-frame x runs far from the axis.
TEST_F(MultipoleTRepTest, ContainmentFollowsTheArcNotTheStraightFrame) {
    const double length = 1.0;
    const double angle  = Physics::pi / 2.0;
    MultipoleTRep magnet = makeCurved({-1.0}, length, angle, 0.02);
    magnet.setAperture(ApertureType::ELLIPTICAL, {0.05, 0.05});

    // On the design orbit near the exit: far from the entrance axis in x, but on the
    // orbit in arc coordinates.
    const Vector3 onOrbit = pointAtArc(0.95 * length, length, angle);
    EXPECT_GT(std::abs(onOrbit(0)), 0.05);
    EXPECT_TRUE(magnet.isInside(onOrbit));

    // Past the exit fringe it is no longer selected.
    EXPECT_FALSE(magnet.isInside(pointAtArc(length + 0.5, length, angle)));
}

// The reference particle sees the same field; it reports a loss only when it is
// inside the longitudinal window but outside the aperture.
TEST_F(MultipoleTRepTest, ReferenceParticleSeesTheFieldAndReportsLosses) {
    const double length = 1.0;
    MultipoleTRep magnet = makeStraight({-0.5}, length, 0.02);
    magnet.setAperture(ApertureType::RECTANGULAR, {0.02, 0.01});

    Vector3 E(0.0);
    Vector3 B(0.0);
    EXPECT_FALSE(magnet.applyToReferenceParticle(
            Vector3(0.0, 0.0, 0.5 * length), Vector3(0.0), 0.0, E, B));
    EXPECT_NEAR(B(1), -0.5, 1.0e-9);

    // Outside the aperture, inside the window: lost, and no field added.
    B = Vector3(0.0);
    EXPECT_TRUE(magnet.applyToReferenceParticle(
            Vector3(0.05, 0.0, 0.5 * length), Vector3(0.0), 0.0, E, B));
    EXPECT_NEAR(B(1), 0.0, 1.0e-12);

    // Outside the window: not selected, not a loss.
    B = Vector3(0.0);
    EXPECT_FALSE(
            magnet.applyToReferenceParticle(Vector3(0.05, 0.0, -1.0), Vector3(0.0), 0.0, E, B));
    EXPECT_NEAR(B(1), 0.0, 1.0e-12);
}

// ---------------------------------------------------------------------------
// The scaling model
// ---------------------------------------------------------------------------

namespace {
    /// An empty line, only so a DefaultVisitor can be built to accept() the magnet.
    class DummyBeamline final : public Beamline {
    public:
        DummyBeamline() : Beamline("dummy") {}

        ElementType getType() const override { return ElementType::BEAMLINE; }
        Geometry& getGeometry() override { return geometry_m; }
        const Geometry& getGeometry() const override { return geometry_m; }
        void accept(BeamlineVisitor& visitor) const override { visitor.visitBeamline(*this); }
        ElementBase* clone() const override { return new DummyBeamline(*this); }
        void iterate(BeamlineVisitor&, bool) const override {}

    private:
        Geometry geometry_m{Geometry::makeNull()};
    };
}  // namespace

// A named time dependence scales the whole field; the name is resolved in accept().
TEST_F(MultipoleTRepTest, ScalingModelScalesTheField) {
    // A step from 1 to 2 between t = 1 and t = 1.01.
    AbstractTimeDependence::setTimeDependence(
            "STEPFUNCTION", std::make_shared<SplineTimeDependence>(
                                    1, std::vector<double>{0.0, 1.0, 1.01, 2.0},
                                    std::vector<double>{1.0, 1.0, 2.0, 2.0}));

    const double length = 1.0;
    MultipoleTRep magnet = makeStraight({1.0}, length, 0.02);
    const Vector3 centre(0.0, 0.0, 0.5 * length);

    DummyBeamline line;
    DefaultVisitor visitor(line, false, false);
    magnet.setScalingName("stepfunction");
    EXPECT_EQ(magnet.getScalingName(), "STEPFUNCTION");
    magnet.accept(visitor);
    EXPECT_NEAR(fieldAt(magnet, centre, 0.0)(1), 1.0, 1.0e-6);
    EXPECT_NEAR(fieldAt(magnet, centre, 1.2)(1), 2.0, 1.0e-6);

    // An empty name is no scaling at all.
    magnet.setScalingName("");
    magnet.accept(visitor);
    EXPECT_NEAR(fieldAt(magnet, centre, 1.2)(1), 1.0, 1.0e-6);
}

// ---------------------------------------------------------------------------
// Housekeeping
// ---------------------------------------------------------------------------

// A clone carries the field state, the geometry and the resolved scaling model.
TEST_F(MultipoleTRepTest, CloneCarriesTheState) {
    const double length = 1.0;
    const double angle  = Physics::pi / 6.0;
    MultipoleTRep magnet = makeCurved({-1.0, 0.5}, length, angle, 0.03);
    magnet.setAperture(ApertureType::RECTANGULAR, {0.02, 0.01});
    magnet.setScalingName("SCALING");

    std::unique_ptr<ElementBase> copy(magnet.clone());
    auto* clone = dynamic_cast<MultipoleTRep*>(copy.get());
    ASSERT_NE(clone, nullptr);

    EXPECT_EQ(clone->getMaxFOrder(), magnet.getMaxFOrder());
    EXPECT_DOUBLE_EQ(clone->getFringeLeft(), 0.03);
    EXPECT_DOUBLE_EQ(clone->getFringeRight(), 0.03);
    EXPECT_DOUBLE_EQ(clone->getTransverseProfile()[0], -1.0);
    EXPECT_DOUBLE_EQ(clone->getTransverseProfile()[1], 0.5);
    EXPECT_EQ(clone->getScalingName(), "SCALING");
    EXPECT_EQ(clone->getType(), ElementType::MULTIPOLET);
    EXPECT_NEAR(clone->getGeometry().getBendAngle(), angle, 1.0e-12);
    EXPECT_EQ(clone->getAperture().first, ApertureType::RECTANGULAR);

    // The clone reports the same field.
    const Vector3 point = pointAtArc(0.5 * length, length, angle);
    EXPECT_NEAR(fieldAt(*clone, point)(1), fieldAt(magnet, point)(1), 1.0e-12);
}

// The dipole coefficient is the "BY" channel attribute.
TEST_F(MultipoleTRepTest, DipoleChannel) {
    MultipoleTRep magnet = makeStraight({-0.75}, 1.0, 0.02);
    EXPECT_DOUBLE_EQ(magnet.getB(), -0.75);
    EXPECT_DOUBLE_EQ(magnet.getAttribute("BY"), -0.75);

    magnet.setAttribute("BY", 0.25);
    EXPECT_DOUBLE_EQ(magnet.getB(), 0.25);
    EXPECT_NEAR(fieldAt(magnet, Vector3(0.0, 0.0, 0.5))(1), 0.25, 1.0e-9);
}

// The setters reject what the field model cannot represent.
TEST_F(MultipoleTRepTest, SettersRejectOutOfRangeValues) {
    MultipoleTRep magnet("MULTIPOLET");
    EXPECT_THROW(magnet.setMaxFOrder(0), OpalException);
    EXPECT_THROW(
            magnet.setMaxFOrder(MultipoleTFieldModel::MaxFOrder + 1), OpalException);
    EXPECT_NO_THROW(magnet.setMaxFOrder(MultipoleTFieldModel::MaxFOrder));
    EXPECT_THROW(magnet.setFringeLengths(-0.1, 0.1), OpalException);
    EXPECT_THROW(
            magnet.setTransverseProfile(
                    std::vector<double>(MultipoleTFieldModel::NumPoles + 1, 1.0)),
            OpalException);
    EXPECT_NO_THROW(
            magnet.setTransverseProfile(
                    std::vector<double>(MultipoleTFieldModel::NumPoles, 1.0)));
}
