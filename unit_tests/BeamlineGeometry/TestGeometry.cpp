/**
 * \file TestGeometry.cpp
 * \brief Pinning tests for the unified Geometry class.
 *
 * The expected numbers were captured from the original geometry implementation
 * (StraightGeometry / PlanarArcGeometry / RBendGeometry / NullGeometry and the
 * BendBase edge-transform conventions) and are hardcoded here so the test is the
 * regression oracle for the geometry-cleanup refactor: it pins the lengths, bend
 * parameters, chord, design path and CoordinateSystemTrafo edge transforms,
 * independent of the now-deleted legacy classes.
 */

#include "BeamlineGeometry/Geometry.h"

#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/Quaternion.hpp"

#include "gtest/gtest.h"

#include <array>

namespace {
    // Compare a CoordinateSystemTrafo against an expected origin (x,y,z) and an
    // expected 3x3 rotation matrix in row-major order.
    void expectTrafo(
            const CoordinateSystemTrafo& got, const std::array<double, 3>& origin,
            const std::array<double, 9>& rot, double tol = 1e-12) {
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(got.getOrigin()(i), origin[i], tol) << "origin component " << i;
        }
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                EXPECT_NEAR(got.getRotationMatrix()(i, j), rot[3 * i + j], tol)
                        << "rotation (" << i << "," << j << ")";
            }
        }
    }

    const std::array<double, 9> kIdentity = {1, 0, 0, 0, 1, 0, 0, 0, 1};
}  // namespace

TEST(GeometryTest, Straight) {
    const double L = 0.75;
    Geometry g     = Geometry::makeStraight(L);

    EXPECT_DOUBLE_EQ(g.getElementLength(), L);
    EXPECT_DOUBLE_EQ(g.getArcLength(), L);
    EXPECT_DOUBLE_EQ(g.getOrigin(), L / 2.0);
    EXPECT_DOUBLE_EQ(g.getEntrance(), -L / 2.0);
    EXPECT_DOUBLE_EQ(g.getExit(), L / 2.0);
    EXPECT_DOUBLE_EQ(g.getBendAngle(), 0.0);
    EXPECT_DOUBLE_EQ(g.getChordLength(), L);

    expectTrafo(g.getEdgeToBegin(), {0.0, 0.0, 0.0}, kIdentity);
    expectTrafo(g.getEdgeToEnd(), {0.0, 0.0, L}, kIdentity);
}

TEST(GeometryTest, NullIsZeroExtentIdentity) {
    Geometry g = Geometry::makeNull();

    EXPECT_DOUBLE_EQ(g.getElementLength(), 0.0);
    EXPECT_DOUBLE_EQ(g.getArcLength(), 0.0);

    expectTrafo(g.getEdgeToBegin(), {0.0, 0.0, 0.0}, kIdentity);
    expectTrafo(g.getEdgeToEnd(), {0.0, 0.0, 0.0}, kIdentity);
}

TEST(GeometryTest, SectorArc) {
    const double L     = 1.2;
    const double angle = 0.35;
    Geometry g         = Geometry::makeArc(L, angle / L);

    EXPECT_DOUBLE_EQ(g.getElementLength(), L);
    EXPECT_DOUBLE_EQ(g.getArcLength(), L);
    EXPECT_DOUBLE_EQ(g.getBendAngle(), angle);
    EXPECT_DOUBLE_EQ(g.getCurvature(), angle / L);
    EXPECT_NEAR(g.getChordLength(), 1.1938843720703722, 1e-12);

    expectTrafo(
            g.getEdgeToBegin(), {-0.052366152325942467, 0.0, -0.59694218603518612},
            {0.98472653890493334, 0, -0.17410813759359647, 0, 1, 0, 0.17410813759359647, 0,
             0.98472653890493334});
    expectTrafo(
            g.getEdgeToEnd(), {-0.052366152325942467, 0.0, 0.59694218603518612},
            {0.98472653890493334, 0, 0.17410813759359647, 0, 1, 0, -0.17410813759359647, 0,
             0.98472653890493334});
}

TEST(GeometryTest, RectangularBend) {
    const double L     = 0.9;   // straight body length
    const double angle = 0.40;  // full bend angle
    Geometry g         = Geometry::makeRBend(L, angle);

    EXPECT_DOUBLE_EQ(g.getElementLength(), L);
    EXPECT_NEAR(g.getArcLength(), 0.90602811858102206, 1e-12);  // L*(a/2)/sin(a/2)
    EXPECT_DOUBLE_EQ(g.getBendAngle(), angle);
    EXPECT_DOUBLE_EQ(g.getChordLength(), L);

    expectTrafo(
            g.getEdgeToBegin(), {0.0, 0.0, -0.45},
            {0.98006657784124163, 0, -0.19866933079506133, 0, 1, 0, 0.19866933079506133, 0,
             0.98006657784124163});
    expectTrafo(
            g.getEdgeToEnd(), {0.0, 0.0, 0.45},
            {0.98006657784124163, 0, 0.19866933079506133, 0, 1, 0, -0.19866933079506133, 0,
             0.98006657784124163});
}

TEST(GeometryTest, SetElementLengthRecomputesArcAngle) {
    Geometry g = Geometry::makeArc(1.0, 0.2);  // angle = 0.2
    EXPECT_DOUBLE_EQ(g.getBendAngle(), 0.2);
    g.setElementLength(2.0);  // angle = h * len = 0.2 * 2.0
    EXPECT_DOUBLE_EQ(g.getBendAngle(), 0.4);
    EXPECT_DOUBLE_EQ(g.getCurvature(), 0.2);
}

TEST(GeometryTest, DesignPathIsCentredAndOnTheArc) {
    Geometry g = Geometry::makeArc(1.2, 0.35 / 1.2);

    auto path = g.getDesignPath(32);
    ASSERT_GE(path.size(), 32u);

    // Endpoints sit on the centred arc (z runs from -chord/2 to +chord/2).
    EXPECT_NEAR(path.front()(0), -0.052366152325942467, 1e-12);
    EXPECT_NEAR(path.front()(2), -0.59694218603518612, 1e-12);
    EXPECT_NEAR(path.back()(0), -0.052366152325942467, 1e-12);
    EXPECT_NEAR(path.back()(2), 0.59694218603518612, 1e-12);
}
