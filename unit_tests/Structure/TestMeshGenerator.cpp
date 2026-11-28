#include <climits>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include "AbsBeamline/ScalingFFAMagnet.h"
#include "BeamlineCore/SBendRep.h"
#include "BeamlineCore/TravelingWaveRep.h"
#include "Structure/MeshGenerator.h"
#include "gtest/gtest.h"

class MeshGeneratorTest : public ::testing::Test {
protected:
    static void checkTravelingWaveCloneMesh() {
        TravelingWaveRep original("CLICG");
        original.getGeometry().setElementLength(.233284666);
        original.setAperture(ApertureType::ELLIPTICAL, {.05, .05});
        original.setNumCells(26);
        const std::unique_ptr<ElementBase> cloned(original.clone());
        MeshGenerator generator;
        generator.add(original);
        generator.add(*cloned);
        ASSERT_EQ(generator.elements_m.size(), 2u);
        const auto& source = generator.elements_m[0];
        const auto& copy   = generator.elements_m[1];
        ASSERT_FALSE(source.vertices_m.empty());
        ASSERT_FALSE(source.triangles_m.empty());
        ASSERT_EQ(copy.vertices_m.size(), source.vertices_m.size());
        ASSERT_EQ(copy.triangles_m.size(), source.triangles_m.size());
        EXPECT_EQ(copy.type_m, MeshGenerator::TRAVELINGWAVE);
        for (size_t i = 0; i < source.vertices_m.size(); ++i)
            for (int d = 0; d < 3; ++d) {
                EXPECT_TRUE(std::isfinite(copy.vertices_m[i](d)));
                EXPECT_DOUBLE_EQ(copy.vertices_m[i](d), source.vertices_m[i](d));
            }
        for (size_t i = 0; i < source.triangles_m.size(); ++i)
            for (int d = 0; d < 3; ++d)
                EXPECT_EQ(copy.triangles_m[i](d), source.triangles_m[i](d));
        EXPECT_EQ(cloned->getAperture(), original.getAperture());
    }

    static void checkTravelingWaveCells(int cells, double minor, double major) {
        const double length = (cells + 2) * .012;
        const auto mesh     = MeshGenerator::getTravelingWave(length, minor, major, cells);
        ASSERT_FALSE(mesh.vertices_m.empty());
        ASSERT_FALSE(mesh.triangles_m.empty());
        // One maximum-radius plane per accelerating cell; the smaller end
        // flanges do not count. Inspect geometry instead of a vertex-count formula.
        std::map<double, double> radiusByZ;
        for (const auto& v : mesh.vertices_m) {
            for (int d = 0; d < 3; ++d)
                EXPECT_TRUE(std::isfinite(v(d)));
            EXPECT_LE(std::abs(v(0)), major);
            EXPECT_LE(std::abs(v(1)), minor);
            EXPECT_GE(v(2), 0.);
            EXPECT_LE(v(2), length);
            const double radius = std::hypot(v(0) / major, v(1) / minor);
            radiusByZ[v(2)]     = std::max(radiusByZ[v(2)], radius);
        }
        EXPECT_DOUBLE_EQ(radiusByZ.begin()->first, 0.);
        EXPECT_DOUBLE_EQ(radiusByZ.rbegin()->first, length);
        double maxRadius = 0;
        for (const auto& [z, radius] : radiusByZ)
            maxRadius = std::max(maxRadius, radius);
        int bulges = 0;
        for (auto current = std::next(radiusByZ.begin()); std::next(current) != radiusByZ.end();
             ++current) {
            if (current->second > .9 * maxRadius && current->second > std::prev(current)->second
                && current->second > std::next(current)->second) {
                ++bulges;
            }
        }
        EXPECT_EQ(bulges, cells);

        // Every edge has two opposite face directions: no cracks or overlapping
        // per-cell tube walls. Positive signed volume establishes outward winding.
        std::map<std::pair<unsigned, unsigned>, std::pair<int, int>> edges;
        double volume = 0;
        for (const auto& t : mesh.triangles_m) {
            for (int d = 0; d < 3; ++d)
                ASSERT_LT(t(d), mesh.vertices_m.size());
            const auto a = mesh.vertices_m[t(0)], b = mesh.vertices_m[t(1)],
                       c                 = mesh.vertices_m[t(2)];
            const Vector_t<double, 3> ab = b - a, ac = c - a;
            const Vector_t<double, 3> normal = cross(ab, ac), bc = cross(b, c);
            EXPECT_GT(dot(normal, normal), 0.);
            volume += dot(a, bc) / 6.;
            for (int d = 0; d < 3; ++d) {
                const unsigned from = t(d), to = t((d + 1) % 3);
                auto& edge = edges[{std::min(from, to), std::max(from, to)}];
                ++edge.first;
                edge.second += from < to ? 1 : -1;
            }
        }
        for (const auto& [key, incidence] : edges) {
            EXPECT_EQ(incidence.first, 2);
            EXPECT_EQ(incidence.second, 0);
        }
        EXPECT_GT(volume, 0.);

        // The add() path must use NUMCELLS for both supported aperture shapes.
        for (auto shape : {ApertureType::ELLIPTICAL, ApertureType::RECTANGULAR}) {
            TravelingWaveRep cavity("TW");
            cavity.getGeometry().setElementLength(length);
            cavity.setAperture(shape, {minor, major});
            cavity.setNumCells(cells);
            MeshGenerator generator;
            generator.add(cavity);
            ASSERT_EQ(generator.elements_m.size(), 1u);
            EXPECT_EQ(generator.elements_m[0].vertices_m.size(), mesh.vertices_m.size());
            EXPECT_EQ(generator.elements_m[0].type_m, MeshGenerator::TRAVELINGWAVE);
            EXPECT_EQ(cavity.getNumCells(), cells);
            EXPECT_EQ(cavity.getAperture().first, shape);
            EXPECT_EQ(cavity.getAperture().second, (std::vector<double>{minor, major}));
            EXPECT_DOUBLE_EQ(cavity.getGeometry().getElementLength(), length);
        }
    }

    static void checkInvalidTravelingWaveDisplay() {
        for (int cells : {0, -1, INT_MAX}) {
            EXPECT_TRUE(MeshGenerator::getTravelingWave(1., .05, .05, cells).vertices_m.empty());
        }
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (double bad : {0., -1., nan, std::numeric_limits<double>::infinity()}) {
            EXPECT_TRUE(MeshGenerator::getTravelingWave(bad, .05, .05, 26).vertices_m.empty());
            EXPECT_TRUE(MeshGenerator::getTravelingWave(1., bad, .05, 26).vertices_m.empty());
            EXPECT_TRUE(MeshGenerator::getTravelingWave(1., .05, bad, 26).vertices_m.empty());
        }
    }

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

TEST_F(MeshGeneratorTest, TravelingWaveCloneKeepsItsRenderableBody) {
    checkTravelingWaveCloneMesh();
}

TEST_F(MeshGeneratorTest, TravelingWaveHasOneRoundedBulgePerConfiguredCell) {
    for (int cells : {1, 2, 9, 26}) {
        SCOPED_TRACE(cells);
        checkTravelingWaveCells(cells, .02, .03);
        checkTravelingWaveCells(cells, .003, .004);
    }
}

TEST_F(MeshGeneratorTest, InvalidTravelingWaveDisplayDimensionsOrCountsProduceNoMesh) {
    checkInvalidTravelingWaveDisplay();
}
