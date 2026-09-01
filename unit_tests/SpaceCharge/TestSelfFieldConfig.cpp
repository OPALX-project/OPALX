#include <gtest/gtest.h>

#include "SpaceCharge/ParticleSetView.h"
#include "SpaceCharge/SelfFieldConfig.h"
#include "Utilities/OpalException.h"

#include <array>

namespace opalx::spacecharge {
    namespace {

        Pic3DConfigValues p3mValues(BoundaryConditionKind boundary) {
            Pic3DConfigValues values;
            values.backend            = PoissonBackendKind::P3M;
            values.p3mCutoff          = 0.025;
            values.boundaryConditions = {boundary, boundary, boundary};
            return values;
        }

        TEST(SelfFieldConfigTest, P3MSelectsOverlapLayoutAndBoundaryMode) {
            SelfFieldConfig open(Pic3DConfig(p3mValues(BoundaryConditionKind::Open)));
            const ParticleLayoutConfig openLayout = open.particleLayoutConfig();
            EXPECT_EQ(openLayout.kind, ParticleLayoutKind::SpatialOverlap);
            EXPECT_DOUBLE_EQ(openLayout.overlapCutoff, 0.025);
            EXPECT_FALSE(openLayout.periodic);

            SelfFieldConfig periodic(Pic3DConfig(p3mValues(BoundaryConditionKind::Periodic)));
            EXPECT_TRUE(periodic.particleLayoutConfig().periodic);
        }

        TEST(SelfFieldConfigTest, P3MRejectsInvalidCombinations) {
            auto values      = p3mValues(BoundaryConditionKind::Open);
            values.p3mCutoff = 0.0;
            EXPECT_THROW({ Pic3DConfig config(values); }, OpalException);

            values                       = p3mValues(BoundaryConditionKind::Open);
            values.boundaryConditions[1] = BoundaryConditionKind::Periodic;
            EXPECT_THROW({ Pic3DConfig config(values); }, OpalException);

            values = p3mValues(BoundaryConditionKind::Open);
            values.binning.emplace(BinningConfigValues{});
            EXPECT_THROW({ Pic3DConfig config(values); }, OpalException);

            values = p3mValues(BoundaryConditionKind::Open);
            values.correction =
                    CorrectionConfig({.kind = CorrectionKind::ImageCharge, .planeZ = 0.0});
            EXPECT_THROW({ Pic3DConfig config(values); }, OpalException);
        }

        TEST(SelfFieldConfigTest, Pic2d5UsesIndependentAlgorithmAndSpatialLayout) {
            Pic2d5ConfigValues values;
            values.meshSize          = {16, 18, 20};
            values.pipeSizeX         = 0.1;
            values.pipeSizeY         = 0.2;
            values.beamRadius        = 0.01;
            values.referencePathFile = "design-path.dat";
            SelfFieldConfig config{Pic2d5Config(values)};

            EXPECT_EQ(config.algorithmKind(), SelfFieldAlgorithmKind::Pic2d5);
            EXPECT_EQ(config.particleLayoutConfig().kind, ParticleLayoutKind::Spatial);
            EXPECT_FALSE(config.particleLayoutConfig().periodic);
        }

        TEST(ParticleSelectionTest, PoliciesSelectPrimaryOrEveryTrackingActiveContainer) {
            ParticleContainerAttributes attributes;
            std::array containers{
                    ParticleContainerView("primary", attributes, true, true),
                    ParticleContainerView("secondary", attributes, false, true),
                    ParticleContainerView("inactive", attributes, false, false)};
            ParticleSetView particles(containers, 0, 0);

            particles.applySelection(ParticleSelectionPolicy::AllTrackingActive);
            EXPECT_TRUE(containers[0].selectedForSolve());
            EXPECT_TRUE(containers[1].selectedForSolve());
            EXPECT_FALSE(containers[2].selectedForSolve());

            particles.applySelection(ParticleSelectionPolicy::PrimaryOnly);
            EXPECT_TRUE(containers[0].selectedForSolve());
            EXPECT_FALSE(containers[1].selectedForSolve());
            EXPECT_FALSE(containers[2].selectedForSolve());
        }

    }  // namespace
}  // namespace opalx::spacecharge
