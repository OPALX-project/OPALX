#include <gtest/gtest.h>

#include "SpaceCharge/Pic3D/CorrectionPlan.h"
#include "SpaceCharge/SelfFieldDiagnostics.h"

namespace opalx::spacecharge {
    namespace {

        RequestedPhysics request(CorrectionKind correction, double planeZ = 0.5) {
            RequestedPhysics result;
            result.correction = {correction, planeZ};
            return result;
        }

        CorrectionConfig correctionConfig(CorrectionKind kind, std::size_t maximumSteps = 0) {
            return CorrectionConfig({.kind = kind, .planeZ = 0.5, .maximumSteps = maximumSteps});
        }

        TEST(CorrectionPlanTest, EmitsEveryDirectPassSequence) {
            CorrectionPlan none(correctionConfig(CorrectionKind::None), false);
            const PreparedCorrection direct = none.prepare(request(CorrectionKind::None), 0);
            ASSERT_EQ(direct.passCount, 1u);
            EXPECT_EQ(direct.passes[0], SolvePass::DirectPrimary);

            CorrectionPlan image(correctionConfig(CorrectionKind::ImageCharge), false);
            const PreparedCorrection directImage =
                    image.prepare(request(CorrectionKind::ImageCharge), 0);
            ASSERT_EQ(directImage.passCount, 1u);
            EXPECT_EQ(directImage.passes[0], SolvePass::DirectPrimaryWithImage);

            CorrectionPlan shifted(correctionConfig(CorrectionKind::ShiftedGreen), false);
            const PreparedCorrection directShifted =
                    shifted.prepare(request(CorrectionKind::ShiftedGreen), 0);
            ASSERT_EQ(directShifted.passCount, 1u);
            EXPECT_EQ(directShifted.passes[0], SolvePass::DirectPrimary);
            EXPECT_TRUE(directShifted.shiftedIgnored);
        }

        TEST(CorrectionPlanTest, EmitsEveryBinnedPassSequence) {
            CorrectionPlan none(correctionConfig(CorrectionKind::None), true);
            const PreparedCorrection primary = none.prepare(request(CorrectionKind::None), 0);
            ASSERT_EQ(primary.passCount, 1u);
            EXPECT_EQ(primary.passes[0], SolvePass::BinnedPrimary);

            CorrectionPlan image(correctionConfig(CorrectionKind::ImageCharge), true);
            const PreparedCorrection withImage =
                    image.prepare(request(CorrectionKind::ImageCharge), 0);
            ASSERT_EQ(withImage.passCount, 2u);
            EXPECT_EQ(withImage.passes[0], SolvePass::BinnedPrimary);
            EXPECT_EQ(withImage.passes[1], SolvePass::BinnedImage);

            CorrectionPlan shifted(correctionConfig(CorrectionKind::ShiftedGreen), true);
            const PreparedCorrection withShifted =
                    shifted.prepare(request(CorrectionKind::ShiftedGreen), 0);
            ASSERT_EQ(withShifted.passCount, 2u);
            EXPECT_EQ(withShifted.passes[0], SolvePass::BinnedPrimary);
            EXPECT_EQ(withShifted.passes[1], SolvePass::BinnedShiftedGreen);
        }

        TEST(CorrectionPlanTest, ReportsCorrectionExpirationWithoutValidatingAgain) {
            CorrectionPlan plan(correctionConfig(CorrectionKind::ImageCharge, 4), true);
            const PreparedCorrection active = plan.prepare(request(CorrectionKind::ImageCharge), 3);
            EXPECT_FALSE(active.correctionExpired);

            const PreparedCorrection expired = plan.prepare(request(CorrectionKind::None), 4);
            EXPECT_TRUE(expired.correctionExpired);
            EXPECT_EQ(expired.maximumSteps, 4u);
            EXPECT_EQ(expired.configuredCorrection.kind, CorrectionKind::ImageCharge);
            ASSERT_EQ(expired.passCount, 1u);
            EXPECT_EQ(expired.passes[0], SolvePass::BinnedPrimary);
        }

        TEST(SolvePassPolicyTest, DerivesAllBehaviorFromTheExhaustiveTag) {
            using DepositKind = PicScatterGather<double, 3>::DepositKind;

            const auto direct = solvePassPolicy<double, 3>(SolvePass::DirectPrimary, 0.5);
            EXPECT_FALSE(direct.binned);
            EXPECT_EQ(direct.depositKind, DepositKind::Primary);
            EXPECT_EQ(direct.backendRule, BackendSolveRule::Standard);
            EXPECT_FALSE(direct.suppressFieldDump);
            EXPECT_FALSE(direct.dumpDirichletPlaneAfter);

            const auto directImage =
                    solvePassPolicy<double, 3>(SolvePass::DirectPrimaryWithImage, 0.5);
            EXPECT_FALSE(directImage.binned);
            EXPECT_EQ(directImage.depositKind, DepositKind::PrimaryAndImage);
            EXPECT_TRUE(directImage.imagePolicy.enabled);
            EXPECT_DOUBLE_EQ(directImage.imagePolicy.planeZ, 0.5);
            EXPECT_TRUE(directImage.dumpDirichletPlaneAfter);

            const auto primary = solvePassPolicy<double, 3>(SolvePass::BinnedPrimary, 0.5);
            EXPECT_TRUE(primary.binned);
            EXPECT_EQ(primary.depositKind, DepositKind::Primary);
            EXPECT_TRUE(primary.suppressFieldDump);
            EXPECT_DOUBLE_EQ(primary.magneticSign, 1.0);
            EXPECT_EQ(primary.sourceRule, FieldSourceRule::Direct);

            const auto image = solvePassPolicy<double, 3>(SolvePass::BinnedImage, 0.5);
            EXPECT_TRUE(image.binned);
            EXPECT_EQ(image.depositKind, DepositKind::Image);
            EXPECT_TRUE(image.imagePolicy.enabled);
            EXPECT_TRUE(image.suppressFieldDump);
            EXPECT_DOUBLE_EQ(image.magneticSign, -1.0);
            EXPECT_EQ(image.sourceRule, FieldSourceRule::Direct);
            EXPECT_TRUE(image.dumpDirichletPlaneAfter);

            const auto shifted = solvePassPolicy<double, 3>(SolvePass::BinnedShiftedGreen, 0.5);
            EXPECT_TRUE(shifted.binned);
            EXPECT_EQ(shifted.depositKind, DepositKind::Primary);
            EXPECT_EQ(shifted.backendRule, BackendSolveRule::ShiftedGreen);
            EXPECT_FALSE(shifted.suppressFieldDump);
            EXPECT_DOUBLE_EQ(shifted.magneticSign, -1.0);
            EXPECT_EQ(shifted.sourceRule, FieldSourceRule::ShiftedGreenImageZ);
        }

        TEST(SelfFieldDiagnosticsTest, CountsCompletedWorkAndAppliesCadence) {
            SelfFieldDiagnostics diagnostics({.binTableFrequency = 4, .planeDumpFrequency = 3});
            EXPECT_EQ(diagnostics.backendSolveCount(), 0u);
            EXPECT_EQ(diagnostics.redistributionCount(), 0u);
            EXPECT_FALSE(diagnostics.hasCurrentBinCount());

            diagnostics.completeBackendSolve();
            diagnostics.completeBackendSolve();
            diagnostics.completeRedistribution();
            diagnostics.recordBinCount(7);
            EXPECT_EQ(diagnostics.backendSolveCount(), 2u);
            EXPECT_EQ(diagnostics.redistributionCount(), 1u);
            EXPECT_TRUE(diagnostics.hasCurrentBinCount());
            EXPECT_EQ(diagnostics.currentBinCount(), 7u);
            EXPECT_TRUE(diagnostics.shouldPrintBinTable(8));
            EXPECT_FALSE(diagnostics.shouldPrintBinTable(9));
            EXPECT_TRUE(diagnostics.shouldDumpPlane(6));
            EXPECT_FALSE(diagnostics.shouldDumpPlane(7));
            EXPECT_FALSE(diagnostics.shouldDumpPlane(-1));
        }

    }  // namespace
}  // namespace opalx::spacecharge
