#include <cmath>
#include <map>
#include "AbsBeamline/ScalingFFAMagnet.h"
#include "BeamlineCore/SBendRep.h"
#include "Structure/MeshGenerator.h"
#include "gtest/gtest.h"

class MeshGeneratorTest : public ::testing::Test {
protected:
    static void checkFFA(double radius, double spiral) {
        ScalingFFAMagnet magnet("FFA");
        magnet.setR0(radius);
        magnet.setRMin(3.7);
        magnet.setRMax(4.6);
        magnet.setVerticalExtent(.25);
        magnet.setTanDelta(spiral);
        magnet.setPhiStart(.15);
        magnet.setDipoleConstant(1);
        magnet.getGeometry() = Geometry::makeSBend(.8, 1 / radius);
        const auto before    = magnet.getConfig();
        const auto mesh      = MeshGenerator::getScalingFFA(magnet);
        ASSERT_FALSE(mesh.vertices_m.empty());
        ASSERT_FALSE(mesh.triangles_m.empty());
        bool entrance = false, exit = false;
        for (const auto& vertex : mesh.vertices_m) {
            for (int d = 0; d < 3; ++d)
                ASSERT_TRUE(std::isfinite(vertex(d)));
            Vector_t<double, 5> chart;
            before.getCylindricalCoordinates(vertex, chart);
            EXPECT_GE(chart(0), 3.7 - 1e-12);
            EXPECT_LE(chart(0), 4.6 + 1e-12);
            EXPECT_NEAR(std::abs(vertex(1)), .25, 1e-14);
            EXPECT_LE(std::abs(chart(4)), .8 / (2 * std::abs(radius)) + 1e-12);
            entrance |= std::abs(chart(4) + .8 / (2 * std::abs(radius))) < 1e-12;
            exit |= std::abs(chart(4) - .8 / (2 * std::abs(radius))) < 1e-12;
        }
        EXPECT_TRUE(entrance);
        EXPECT_TRUE(exit);
        std::map<std::pair<unsigned, unsigned>, std::pair<int, int>> edges;
        double volume = 0;
        for (const auto& triangle : mesh.triangles_m) {
            for (int d = 0; d < 3; ++d)
                ASSERT_LT(triangle(d), mesh.vertices_m.size());
            const auto a = mesh.vertices_m[triangle(0)], b = mesh.vertices_m[triangle(1)],
                       c                 = mesh.vertices_m[triangle(2)];
            const Vector_t<double, 3> ab = b - a, ac = c - a;
            const Vector_t<double, 3> normal = cross(ab, ac), bc = cross(b, c);
            EXPECT_GT(dot(normal, normal), 0);
            volume += dot(a, bc) / 6;
            for (int d = 0; d < 3; ++d) {
                const unsigned from = triangle(d), to = triangle((d + 1) % 3);
                auto& edge = edges[{std::min(from, to), std::max(from, to)}];
                ++edge.first;
                edge.second += from < to ? 1 : -1;
            }
        }
        // Every edge has two oppositely directed faces: closed, consistently
        // oriented surface. Positive volume establishes outward winding.
        for (const auto& edge : edges) {
            EXPECT_EQ(edge.second.first, 2);
            EXPECT_EQ(edge.second.second, 0);
        }
        EXPECT_GT(volume, 0);
        if (spiral == 0) {
            const double exact = .5 * (4.6 * 4.6 - 3.7 * 3.7) * (.8 / std::abs(radius)) * .5;
            EXPECT_NEAR(volume, exact, .005 * exact);
        }
        EXPECT_DOUBLE_EQ(magnet.getConfig().r0_m, before.r0_m);
        EXPECT_DOUBLE_EQ(magnet.getConfig().tanDelta_m, before.tanDelta_m);
        EXPECT_DOUBLE_EQ(magnet.getGeometry().getArcLength(), .8);

        magnet.setDipoleConstant(0);
        const auto arc = MeshGenerator::getScalingFFA(magnet);
        EXPECT_TRUE(arc.triangles_m.empty());
        ASSERT_GT(arc.vertices_m.size(), 2u);
        for (const auto& vertex : arc.vertices_m) {
            EXPECT_NEAR(std::hypot(vertex(0) + radius, vertex(2)), std::abs(radius), 1e-12);
            EXPECT_DOUBLE_EQ(vertex(1), 0);
        }
    }

    static void checkFFAPlacement() {
        ScalingFFAMagnet magnet("FFA");
        magnet.setR0(4);
        magnet.setRMin(3.9);
        magnet.setRMax(4.2);
        magnet.setVerticalExtent(.1);
        magnet.setDipoleConstant(1);
        magnet.getGeometry() = Geometry::makeSBend(.4, .25);
        const auto local     = MeshGenerator::getScalingFFA(magnet);
        const Vector_t<double, 3> origin(1, 2, 3);
        magnet.setCSTrafoGlobal2Local(CoordinateSystemTrafo(origin, Quaternion()));
        MeshGenerator generator;
        generator.add(magnet);
        ASSERT_EQ(generator.elements_m.size(), 1u);
        const auto& global = generator.elements_m.front();
        ASSERT_EQ(global.vertices_m.size(), local.vertices_m.size());
        EXPECT_EQ(global.type_m, MeshGenerator::SCALINGFFA);
        for (size_t i = 0; i < local.vertices_m.size(); ++i)
            for (int d = 0; d < 3; ++d)
                EXPECT_NEAR(global.vertices_m[i](d), local.vertices_m[i](d) + origin(d), 1e-12);
        magnet.setDipoleConstant(0);
        generator.add(magnet);
        EXPECT_EQ(generator.elements_m.back().type_m, MeshGenerator::OTHER);
        EXPECT_TRUE(generator.elements_m.back().triangles_m.empty());
    }

    static void checkProbePlane() {
        // Deliberately translated, oblique global endpoints: no placement,
        // mm conversion, or lab Y/Z exchange may be applied a second time.
        const Vector_t<double, 3> a(-7, 0, 2), b(-5, 0, 3);
        MeshGenerator generator;
        generator.addProbe(a, b);
        ASSERT_EQ(generator.elements_m.size(), 1u);
        const auto& plane = generator.elements_m.front();
        EXPECT_EQ(plane.type_m, MeshGenerator::PROBE);
        ASSERT_EQ(plane.vertices_m.size(), 4u);
        ASSERT_EQ(plane.triangles_m.size(), 2u);
        for (unsigned i = 0; i < 4; ++i) {
            const auto expected = (i == 0 || i == 3) ? a : b;
            EXPECT_DOUBLE_EQ(plane.vertices_m[i](0), expected(0));
            EXPECT_DOUBLE_EQ(plane.vertices_m[i](2), expected(2));
            EXPECT_DOUBLE_EQ(plane.vertices_m[i](1), i < 2 ? -.5 : .5);
        }
        for (const auto& tri : plane.triangles_m) {
            for (int d = 0; d < 3; ++d)
                ASSERT_LT(tri(d), 4u);
            const Vector_t<double, 3> u      = plane.vertices_m[tri(1)] - plane.vertices_m[tri(0)],
                                      v      = plane.vertices_m[tri(2)] - plane.vertices_m[tri(0)];
            const Vector_t<double, 3> normal = cross(u, v);
            EXPECT_GT(dot(normal, normal), 0);
        }
    }

    static void checkBend(double horizontal, double vertical, bool finite, bool shared) {
        SBendRep bend("BEND");
        bend.getGeometry() = Geometry::makeSBend(1.0, 0.5);
        bend.setAperture(
                ApertureType::RECTANGULAR,
                finite ? std::vector<double>{horizontal, vertical} : std::vector<double>{1e6, 1e6});
        const auto aperture = bend.getAperture();
        MeshGenerator generator;
        if (shared) generator.setDriftReference(horizontal, vertical);
        generator.add(bend);
        ASSERT_EQ(generator.elements_m.size(), 1u);
        const auto& mesh = generator.elements_m.front();
        ASSERT_FALSE(mesh.vertices_m.empty());
        EXPECT_FALSE(mesh.triangles_m.empty());
        EXPECT_EQ(mesh.type_m, MeshGenerator::DIPOLE);
        EXPECT_DOUBLE_EQ(mesh.vertices_m.front()(0), horizontal);
        EXPECT_DOUBLE_EQ(mesh.vertices_m.front()(1), vertical);
        for (const auto& vertex : mesh.vertices_m)
            for (int d = 0; d < 3; ++d)
                EXPECT_TRUE(std::isfinite(vertex(d)));
        EXPECT_EQ(bend.getAperture(), aperture);
        EXPECT_DOUBLE_EQ(bend.getGeometry().getCurvature(), 0.5);
    }
};

TEST_F(MeshGeneratorTest, UnboundedApertureGetsDisplayDefault) {
    checkBend(0.05, 0.05, false, false);
}
TEST_F(MeshGeneratorTest, UnboundedApertureUsesSharedSupport) { checkBend(0.2, 0.1, false, true); }
TEST_F(MeshGeneratorTest, FiniteApertureKeepsItsDimensions) { checkBend(0.03, 0.02, true, false); }

TEST_F(MeshGeneratorTest, ScalingFFASectorsHaveClosedOutwardSurfacesForBothBendDirections) {
    for (double radius : {4., -4.})
        for (double spiral : {0., .8, -.8})
            checkFFA(radius, spiral);
}

TEST_F(MeshGeneratorTest, ScalingFFAPlacementTransformsTheMeshOnce) { checkFFAPlacement(); }

TEST_F(MeshGeneratorTest, ProbePlaneUsesGlobalEndpointsAndDisplayOnlyHeight) { checkProbePlane(); }
