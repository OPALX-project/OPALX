#include <gtest/gtest.h>

#include "SpaceCharge/CartesianPIC/CorrectionPassSchedule.h"
#include "SpaceCharge/SpaceChargeDiagnostics.h"

namespace opalx::spacecharge {
    namespace {

        SpaceChargeRequest request(SpaceChargeCorrectionType correction, double planeZ = 0.5) {
            SpaceChargeRequest result;
            result.correction = {correction, planeZ};
            return result;
        }

        CorrectionConfig correctionConfig(
                SpaceChargeCorrectionType kind, std::size_t maximumSteps = 0) {
            return CorrectionConfig({.kind = kind, .planeZ = 0.5, .maximumSteps = maximumSteps});
        }

        TEST(CorrectionPassScheduleTest, EmitsEveryDirectPassSequence) {
            CorrectionPassSchedule none(correctionConfig(SpaceChargeCorrectionType::None), false);
            const CorrectionPassSequence direct =
                    none.passesForStep(request(SpaceChargeCorrectionType::None), 0);
            ASSERT_EQ(direct.passCount, 1u);
            EXPECT_EQ(direct.passes[0], CartesianPICPass::DirectPrimary);

            CorrectionPassSchedule image(
                    correctionConfig(SpaceChargeCorrectionType::ImageCharge), false);
            const CorrectionPassSequence directImage =
                    image.passesForStep(request(SpaceChargeCorrectionType::ImageCharge), 0);
            ASSERT_EQ(directImage.passCount, 1u);
            EXPECT_EQ(directImage.passes[0], CartesianPICPass::DirectPrimaryWithImage);

            CorrectionPassSchedule shifted(
                    correctionConfig(SpaceChargeCorrectionType::ShiftedGreen), false);
            const CorrectionPassSequence directShifted =
                    shifted.passesForStep(request(SpaceChargeCorrectionType::ShiftedGreen), 0);
            ASSERT_EQ(directShifted.passCount, 1u);
            EXPECT_EQ(directShifted.passes[0], CartesianPICPass::DirectPrimary);
            EXPECT_TRUE(directShifted.shiftedIgnored);
        }

        TEST(CorrectionPassScheduleTest, EmitsEveryBinnedPassSequence) {
            CorrectionPassSchedule none(correctionConfig(SpaceChargeCorrectionType::None), true);
            const CorrectionPassSequence primary =
                    none.passesForStep(request(SpaceChargeCorrectionType::None), 0);
            ASSERT_EQ(primary.passCount, 1u);
            EXPECT_EQ(primary.passes[0], CartesianPICPass::BinnedPrimary);

            CorrectionPassSchedule image(
                    correctionConfig(SpaceChargeCorrectionType::ImageCharge), true);
            const CorrectionPassSequence withImage =
                    image.passesForStep(request(SpaceChargeCorrectionType::ImageCharge), 0);
            ASSERT_EQ(withImage.passCount, 2u);
            EXPECT_EQ(withImage.passes[0], CartesianPICPass::BinnedPrimary);
            EXPECT_EQ(withImage.passes[1], CartesianPICPass::BinnedImage);

            CorrectionPassSchedule shifted(
                    correctionConfig(SpaceChargeCorrectionType::ShiftedGreen), true);
            const CorrectionPassSequence withShifted =
                    shifted.passesForStep(request(SpaceChargeCorrectionType::ShiftedGreen), 0);
            ASSERT_EQ(withShifted.passCount, 2u);
            EXPECT_EQ(withShifted.passes[0], CartesianPICPass::BinnedPrimary);
            EXPECT_EQ(withShifted.passes[1], CartesianPICPass::BinnedShiftedGreen);
        }

        TEST(CorrectionPassScheduleTest, ReportsCorrectionExpirationWithoutValidatingAgain) {
            CorrectionPassSchedule plan(
                    correctionConfig(SpaceChargeCorrectionType::ImageCharge, 4), true);
            const CorrectionPassSequence active =
                    plan.passesForStep(request(SpaceChargeCorrectionType::ImageCharge), 3);
            EXPECT_FALSE(active.correctionExpired);

            const CorrectionPassSequence expired =
                    plan.passesForStep(request(SpaceChargeCorrectionType::None), 4);
            EXPECT_TRUE(expired.correctionExpired);
            EXPECT_EQ(expired.maximumSteps, 4u);
            EXPECT_EQ(expired.configuredCorrection.kind, SpaceChargeCorrectionType::ImageCharge);
            ASSERT_EQ(expired.passCount, 1u);
            EXPECT_EQ(expired.passes[0], CartesianPICPass::BinnedPrimary);
        }

        TEST(CartesianPICPassPropertiesTest, DerivesAllBehaviorFromTheExhaustiveTag) {
            using DepositKind = ParticleMeshFieldTransfer<double, 3>::DepositKind;

            const auto direct =
                    cartesianPICPassProperties<double, 3>(CartesianPICPass::DirectPrimary, 0.5);
            EXPECT_FALSE(direct.binned);
            EXPECT_EQ(direct.depositKind, DepositKind::Primary);
            EXPECT_EQ(direct.backendRule, BackendSolveRule::Standard);
            EXPECT_FALSE(direct.suppressFieldDump);
            EXPECT_FALSE(direct.dumpDirichletPlaneAfter);

            const auto directImage = cartesianPICPassProperties<double, 3>(
                    CartesianPICPass::DirectPrimaryWithImage, 0.5);
            EXPECT_FALSE(directImage.binned);
            EXPECT_EQ(directImage.depositKind, DepositKind::PrimaryAndImage);
            EXPECT_TRUE(directImage.imagePolicy.enabled);
            EXPECT_DOUBLE_EQ(directImage.imagePolicy.planeZ, 0.5);
            EXPECT_TRUE(directImage.dumpDirichletPlaneAfter);

            const auto primary =
                    cartesianPICPassProperties<double, 3>(CartesianPICPass::BinnedPrimary, 0.5);
            EXPECT_TRUE(primary.binned);
            EXPECT_EQ(primary.depositKind, DepositKind::Primary);
            EXPECT_TRUE(primary.suppressFieldDump);
            EXPECT_DOUBLE_EQ(primary.magneticSign, 1.0);
            EXPECT_EQ(primary.sourceRule, FieldSourceRule::Direct);

            const auto image =
                    cartesianPICPassProperties<double, 3>(CartesianPICPass::BinnedImage, 0.5);
            EXPECT_TRUE(image.binned);
            EXPECT_EQ(image.depositKind, DepositKind::Image);
            EXPECT_TRUE(image.imagePolicy.enabled);
            EXPECT_TRUE(image.suppressFieldDump);
            EXPECT_DOUBLE_EQ(image.magneticSign, -1.0);
            EXPECT_EQ(image.sourceRule, FieldSourceRule::Direct);
            EXPECT_TRUE(image.dumpDirichletPlaneAfter);

            const auto shifted = cartesianPICPassProperties<double, 3>(
                    CartesianPICPass::BinnedShiftedGreen, 0.5);
            EXPECT_TRUE(shifted.binned);
            EXPECT_EQ(shifted.depositKind, DepositKind::Primary);
            EXPECT_EQ(shifted.backendRule, BackendSolveRule::ShiftedGreen);
            EXPECT_FALSE(shifted.suppressFieldDump);
            EXPECT_DOUBLE_EQ(shifted.magneticSign, -1.0);
            EXPECT_EQ(shifted.sourceRule, FieldSourceRule::ShiftedGreenImageZ);
        }

        TEST(SpaceChargeDiagnosticsTest, CountsCompletedWorkAndAppliesCadence) {
            SpaceChargeDiagnostics diagnostics({.binTableFrequency = 4, .planeDumpFrequency = 3});
            EXPECT_EQ(diagnostics.backendSolveCount(), 0u);
            EXPECT_EQ(diagnostics.redistributionCount(), 0u);
            EXPECT_FALSE(diagnostics.hasCurrentBinCount());

            diagnostics.recordBackendSolve();
            diagnostics.recordBackendSolve();
            diagnostics.recordRedistribution();
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
