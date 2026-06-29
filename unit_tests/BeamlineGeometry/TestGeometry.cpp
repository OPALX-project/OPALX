/**
 * \file TestGeometry.cpp
 * \brief Pinning tests for the unified Geometry class.
 *
 * These tests assert that the new single Geometry class reproduces, bit-for-bit
 * (within tolerance), the behaviour of the legacy geometry subclasses
 * (StraightGeometry, PlanarArcGeometry, RBendGeometry, NullGeometry) and the
 * former ElementBase / BendBase edge-transform conventions. They are the
 * regression oracle for the geometry-cleanup refactor: the legacy classes still
 * exist at this stage, so each new value is compared against the value the old
 * code produced.
 */

#include "BeamlineGeometry/Geometry.h"
#include "BeamlineGeometry/NullGeometry.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "BeamlineGeometry/RBendGeometry.h"
#include "BeamlineGeometry/StraightGeometry.h"

#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/Quaternion.hpp"

#include "gtest/gtest.h"

#include <cmath>

namespace {
    // Same conversion the legacy BendBase used to turn a Euclid3D body frame into
    // the modern CoordinateSystemTrafo. Used here to build the reference value.
    CoordinateSystemTrafo toCST(const Euclid3D& frame) {
        matrix3x3_t rotation;
        const Rotation3D& euclidRotation = frame.getRotation();
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < 3; ++col) {
                rotation(row, col) = euclidRotation(row, col);
            }
        }
        const Vector3D& d = frame.getVector();
        const Vector_t<double, 3> origin(d.getX(), d.getY(), d.getZ());
        return CoordinateSystemTrafo(origin, Quaternion(rotation).conjugate());
    }

    void expectTrafoEq(
            const CoordinateSystemTrafo& got, const CoordinateSystemTrafo& ref, double tol = 1e-12) {
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(got.getOrigin()(i), ref.getOrigin()(i), tol) << "origin component " << i;
        }
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                EXPECT_NEAR(got.getRotationMatrix()(i, j), ref.getRotationMatrix()(i, j), tol)
                        << "rotation (" << i << "," << j << ")";
            }
        }
    }

    CoordinateSystemTrafo identityTrafo() {
        return CoordinateSystemTrafo(Vector_t<double, 3>({0.0, 0.0, 0.0}), Quaternion(1, 0, 0, 0));
    }
}  // namespace

TEST(GeometryTest, StraightMatchesLegacyAndDefaults) {
    const double L = 0.75;
    StraightGeometry oldG(L);
    Geometry g = Geometry::makeStraight(L);

    EXPECT_DOUBLE_EQ(g.getElementLength(), oldG.getElementLength());
    EXPECT_DOUBLE_EQ(g.getArcLength(), oldG.getArcLength());
    EXPECT_DOUBLE_EQ(g.getOrigin(), oldG.getOrigin());
    EXPECT_DOUBLE_EQ(g.getEntrance(), oldG.getEntrance());
    EXPECT_DOUBLE_EQ(g.getExit(), oldG.getExit());
    EXPECT_DOUBLE_EQ(g.getBendAngle(), 0.0);

    // Straight edge transforms match the former ElementBase defaults.
    expectTrafoEq(g.getEdgeToBegin(), identityTrafo());
    expectTrafoEq(
            g.getEdgeToEnd(),
            CoordinateSystemTrafo(Vector_t<double, 3>({0.0, 0.0, L}), Quaternion(1, 0, 0, 0)));
}

TEST(GeometryTest, NullIsZeroExtentIdentity) {
    NullGeometry oldG;
    Geometry g = Geometry::makeNull();

    EXPECT_DOUBLE_EQ(g.getElementLength(), oldG.getElementLength());
    EXPECT_DOUBLE_EQ(g.getArcLength(), oldG.getArcLength());
    EXPECT_DOUBLE_EQ(g.getElementLength(), 0.0);

    expectTrafoEq(g.getEdgeToBegin(), identityTrafo());
    expectTrafoEq(g.getEdgeToEnd(), identityTrafo());
}

TEST(GeometryTest, ArcMatchesPlanarArcGeometry) {
    const double L     = 1.2;
    const double angle = 0.35;
    const double h     = angle / L;  // curvature

    PlanarArcGeometry oldG(L, h);
    Geometry g = Geometry::makeArc(L, h);

    EXPECT_DOUBLE_EQ(g.getElementLength(), oldG.getElementLength());
    EXPECT_DOUBLE_EQ(g.getArcLength(), oldG.getArcLength());
    EXPECT_DOUBLE_EQ(g.getBendAngle(), oldG.getBendAngle());
    EXPECT_DOUBLE_EQ(g.getCurvature(), oldG.getCurvature());

    // Edge transforms reproduce the legacy bend convention (centred arc frames).
    expectTrafoEq(g.getEdgeToBegin(), toCST(oldG.getEntranceFrame()));
    expectTrafoEq(g.getEdgeToEnd(), toCST(oldG.getExitFrame()));

    // Chord matches the legacy frame-distance computation.
    const Vector3D delta = oldG.getExitFrame().getVector() - oldG.getEntranceFrame().getVector();
    const double chordRef =
            std::sqrt(delta.getX() * delta.getX() + delta.getZ() * delta.getZ());
    EXPECT_NEAR(g.getChordLength(), chordRef, 1e-12);
}

TEST(GeometryTest, RBendMatchesRBendGeometry) {
    const double L     = 0.9;   // straight body length
    const double angle = 0.40;  // full bend angle

    RBendGeometry oldG(L, angle);
    Geometry g = Geometry::makeRBend(L, angle);

    EXPECT_DOUBLE_EQ(g.getElementLength(), oldG.getElementLength());
    EXPECT_DOUBLE_EQ(g.getArcLength(), oldG.getArcLength());  // len * (a/2)/sin(a/2)
    EXPECT_DOUBLE_EQ(g.getBendAngle(), oldG.getBendAngle());

    expectTrafoEq(g.getEdgeToBegin(), toCST(oldG.getEntranceFrame()));
    expectTrafoEq(g.getEdgeToEnd(), toCST(oldG.getExitFrame()));
}

TEST(GeometryTest, SetElementLengthRecomputesArcAngle) {
    Geometry g = Geometry::makeArc(1.0, 0.2);  // angle = 0.2
    EXPECT_DOUBLE_EQ(g.getBendAngle(), 0.2);
    g.setElementLength(2.0);  // angle = h * len = 0.2 * 2.0
    EXPECT_DOUBLE_EQ(g.getBendAngle(), 0.4);
    EXPECT_DOUBLE_EQ(g.getCurvature(), 0.2);
}

TEST(GeometryTest, DesignPathEndpointsMatchArcFrames) {
    const double L = 1.0;
    const double h = 0.3;
    Geometry g     = Geometry::makeArc(L, h);

    auto path = g.getDesignPath(32);
    ASSERT_GE(path.size(), 32u);

    PlanarArcGeometry oldG(L, h);
    const Euclid3D entrance = oldG.getTransform(oldG.getEntrance());
    const Euclid3D exit     = oldG.getTransform(oldG.getExit());

    EXPECT_NEAR(path.front()(0), entrance.getX(), 1e-12);
    EXPECT_NEAR(path.front()(2), entrance.getZ(), 1e-12);
    EXPECT_NEAR(path.back()(0), exit.getX(), 1e-12);
    EXPECT_NEAR(path.back()(2), exit.getZ(), 1e-12);
}
