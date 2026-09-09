#include "Structure/MeshGenerator.h"
#include "BeamlineCore/SBendRep.h"
#include "gtest/gtest.h"
#include <cmath>

class MeshGeneratorTest : public ::testing::Test {
protected:
    static void checkBend(double horizontal, double vertical, bool finite, bool shared) {
        SBendRep bend("BEND");
        bend.getGeometry() = Geometry::makeSBend(1.0, 0.5);
        bend.setAperture(ApertureType::RECTANGULAR,
                         finite ? std::vector<double>{horizontal, vertical}
                                : std::vector<double>{1e6, 1e6});
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
            for (int d = 0; d < 3; ++d) EXPECT_TRUE(std::isfinite(vertex(d)));
        EXPECT_EQ(bend.getAperture(), aperture);
        EXPECT_DOUBLE_EQ(bend.getGeometry().getCurvature(), 0.5);
    }
};

TEST_F(MeshGeneratorTest, UnboundedApertureGetsDisplayDefault) {
    checkBend(0.05, 0.05, false, false);
}
TEST_F(MeshGeneratorTest, UnboundedApertureUsesSharedSupport) {
    checkBend(0.2, 0.1, false, true);
}
TEST_F(MeshGeneratorTest, FiniteApertureKeepsItsDimensions) {
    checkBend(0.03, 0.02, true, false);
}
