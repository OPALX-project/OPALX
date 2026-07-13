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
#include <cmath>

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
    Geometry g         = Geometry::makeSBend(L, angle / L);

    EXPECT_DOUBLE_EQ(g.getElementLength(), L);
    EXPECT_DOUBLE_EQ(g.getArcLength(), L);
    EXPECT_DOUBLE_EQ(g.getBendAngle(), angle);
    EXPECT_DOUBLE_EQ(g.getCurvature(), angle / L);
    EXPECT_NEAR(g.getChordLength(), 1.1938843720703722, 1e-12);

    // Entrance-anchored: the entrance edge is the stored frame (identity); the exit edge sits
    // at framePosition(len) with the tangent turned by the full bend angle.
    expectTrafo(g.getEdgeToBegin(), {0.0, 0.0, 0.0}, kIdentity);
    expectTrafo(
            g.getEdgeToEnd(), {-0.20786498452327237, 0.0, 1.1756496255615476},
            {0.9393727128473789, 0, 0.3428978074554514, 0, 1, 0, -0.3428978074554514, 0,
             0.9393727128473789});
}

TEST(GeometryTest, RectangularBend) {
    const double L     = 0.9;   // straight body length
    const double angle = 0.40;  // full bend angle
    Geometry g         = Geometry::makeRBend(L, angle);

    EXPECT_DOUBLE_EQ(g.getElementLength(), L);
    EXPECT_NEAR(g.getArcLength(), 0.90602811858102206, 1e-12);  // L*(a/2)/sin(a/2)
    EXPECT_DOUBLE_EQ(g.getBendAngle(), angle);
    EXPECT_DOUBLE_EQ(g.getChordLength(), L);

    // Entrance-anchored: a rectangular bend is a straight box with parallel faces, so both edge
    // transforms are face-parallel — begin is the identity and end is a pure +z shift by the box
    // length. (The reference orbit meets the faces at half the bend angle, but that orbit tangent
    // is not part of the body edge geometry.)
    expectTrafo(g.getEdgeToBegin(), {0.0, 0.0, 0.0}, kIdentity);
    expectTrafo(g.getEdgeToEnd(), {0.0, 0.0, L}, kIdentity);
}

TEST(GeometryTest, SetElementLengthRecomputesArcAngle) {
    Geometry g = Geometry::makeSBend(1.0, 0.2);  // angle = 0.2
    EXPECT_DOUBLE_EQ(g.getBendAngle(), 0.2);
    g.setElementLength(2.0);  // angle = h * len = 0.2 * 2.0
    EXPECT_DOUBLE_EQ(g.getBendAngle(), 0.4);
    EXPECT_DOUBLE_EQ(g.getCurvature(), 0.2);
}

TEST(GeometryTest, DesignPathIsEntranceAnchoredAndOnTheArc) {
    Geometry g = Geometry::makeSBend(1.2, 0.35 / 1.2);

    auto path = g.getDesignPath(32);
    ASSERT_GE(path.size(), 32u);

    // Entrance-anchored: the path starts at the entrance origin and ends at the exit.
    EXPECT_NEAR(path.front()(0), 0.0, 1e-12);
    EXPECT_NEAR(path.front()(2), 0.0, 1e-12);
    EXPECT_NEAR(path.back()(0), -0.20786498452327237, 1e-12);
    EXPECT_NEAR(path.back()(2), 1.1756496255615476, 1e-12);
}

/* ===================== GeometryHelper (device-callable functions) =============== */

TEST(GeometryHelperTest, SafeAbsCosFloorsNearPerpendicular) {
    EXPECT_DOUBLE_EQ(GeometryHelper::safeAbsCos(0.0), 1.0);
    EXPECT_NEAR(GeometryHelper::safeAbsCos(M_PI / 2.0), 1.0e-6, 1.0e-9);  // floored
    EXPECT_NEAR(GeometryHelper::safeAbsCos(M_PI), 1.0, 1.0e-12);          // |cos|
}

TEST(GeometryHelperTest, RotateAboutYRoundTripAndKnownValue) {
    const Vector_t<double, 3> v(0.3, -0.2, 0.9);
    const Vector_t<double, 3> back =
            GeometryHelper::rotateAboutY(GeometryHelper::rotateAboutY(v, 0.4), -0.4);
    for (unsigned d = 0; d < 3; ++d) {
        EXPECT_NEAR(back(d), v(d), 1e-14);
    }
    // +z rotates towards +x for a positive angle; y is untouched.
    const Vector_t<double, 3> r = GeometryHelper::rotateAboutY(Vector_t<double, 3>(0, 1, 1), 0.5);
    EXPECT_NEAR(r(0), std::sin(0.5), 1e-14);
    EXPECT_DOUBLE_EQ(r(1), 1.0);
    EXPECT_NEAR(r(2), std::cos(0.5), 1e-14);
}

// toBendArcCoords inverts the entrance-anchored design arc: points on the arc map
// to (0, y, s); the straight tangents continue s linearly outside the body.
TEST(GeometryHelperTest, ToBendArcCoordsInvertsTheDesignArc) {
    const double h = 0.7;  // curvature
    const double L = 1.1;  // body (arc) length

    // On-arc point at arc length s.
    auto onArc = [&](double s) {
        return Vector_t<double, 3>((std::cos(h * s) - 1.0) / h, 0.25, std::sin(h * s) / h);
    };
    for (const double s : {0.15, 0.5 * L, 0.97 * L}) {
        const Vector_t<double, 3> arc = GeometryHelper::toBendArcCoords(onArc(s), h, L);
        EXPECT_NEAR(arc(0), 0.0, 1e-12);
        EXPECT_DOUBLE_EQ(arc(1), 0.25);
        EXPECT_NEAR(arc(2), s, 1e-12);
    }

    // A radially displaced point keeps its offset and arc length.
    const double s0 = 0.6, dx = 0.02;
    const Vector_t<double, 3> displaced(
            (1.0 / h + dx) * std::cos(h * s0) - 1.0 / h, 0.0, (1.0 / h + dx) * std::sin(h * s0));
    const Vector_t<double, 3> arcDisplaced = GeometryHelper::toBendArcCoords(displaced, h, L);
    EXPECT_NEAR(arcDisplaced(0), dx, 1e-12);
    EXPECT_NEAR(arcDisplaced(2), s0, 1e-12);

    // Upstream of the entrance: straight tangent, s = z.
    const Vector_t<double, 3> upstream =
            GeometryHelper::toBendArcCoords(Vector_t<double, 3>(0.0, 0.0, -0.3), h, L);
    EXPECT_NEAR(upstream(2), -0.3, 1e-15);

    // Past the exit: s continues along the straight exit tangent.
    const double t = 0.2;
    const Vector_t<double, 3> pastExit =
            onArc(L) + t * Vector_t<double, 3>(-std::sin(h * L), 0.0, std::cos(h * L));
    const Vector_t<double, 3> arcPast = GeometryHelper::toBendArcCoords(pastExit, h, L);
    EXPECT_NEAR(arcPast(0), 0.0, 1e-12);
    EXPECT_NEAR(arcPast(2), L + t, 1e-12);

    // Zero curvature: identity (straight body).
    const Vector_t<double, 3> p(0.1, 0.2, 0.3);
    const Vector_t<double, 3> straight = GeometryHelper::toBendArcCoords(p, 0.0, L);
    for (unsigned d = 0; d < 3; ++d) {
        EXPECT_DOUBLE_EQ(straight(d), p(d));
    }
}

TEST(GeometryHelperTest, RotateArcFieldToEntryRotatesInPlaneOnly) {
    const double h = 0.7, L = 1.1, s = 0.6;

    // Vertical component is invariant.
    const Vector_t<double, 3> vertical =
            GeometryHelper::rotateArcFieldToEntry(Vector_t<double, 3>(0, 1, 0), s, h, L);
    EXPECT_DOUBLE_EQ(vertical(0), 0.0);
    EXPECT_DOUBLE_EQ(vertical(1), 1.0);
    EXPECT_DOUBLE_EQ(vertical(2), 0.0);

    // A radial unit vector rotates by the tangent angle h*s.
    const Vector_t<double, 3> radial =
            GeometryHelper::rotateArcFieldToEntry(Vector_t<double, 3>(1, 0, 0), s, h, L);
    EXPECT_NEAR(radial(0), std::cos(h * s), 1e-14);
    EXPECT_NEAR(radial(2), std::sin(h * s), 1e-14);

    // Upstream: no rotation. Past the exit: frozen at the exit angle.
    const Vector_t<double, 3> upstream =
            GeometryHelper::rotateArcFieldToEntry(Vector_t<double, 3>(1, 0, 0), -0.1, h, L);
    EXPECT_DOUBLE_EQ(upstream(0), 1.0);
    const Vector_t<double, 3> past =
            GeometryHelper::rotateArcFieldToEntry(Vector_t<double, 3>(1, 0, 0), L + 5.0, h, L);
    EXPECT_NEAR(past(0), std::cos(h * L), 1e-14);
    EXPECT_NEAR(past(2), std::sin(h * L), 1e-14);
}
