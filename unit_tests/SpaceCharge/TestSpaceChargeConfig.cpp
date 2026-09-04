#include <gtest/gtest.h>

#include "Attributes/Attributes.h"
#include "SpaceCharge/ParticleFieldSet.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeConfigBuilder.h"
#include "SpaceCharge/SpaceChargeRequestSchedule.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <array>

namespace opalx::spacecharge {
    namespace {

        class TestableFieldSolverCmd : public FieldSolverCmd {
        public:
            void setType(const std::string& type) {
                Attributes::setPredefinedString(itsAttr[FIELDSOLVER::TYPE], type);
            }
        };

        CartesianPICConfig::Parameters p3mValues(FieldBoundaryCondition boundary) {
            CartesianPICConfig::Parameters values;
            values.backend            = PoissonSolverType::P3M;
            values.p3mCutoff          = 0.025;
            values.boundaryConditions = {boundary, boundary, boundary};
            return values;
        }

        CartesianDomainConfig3D storageFor(const CartesianPICConfig::Parameters& values) {
            CartesianDomainConfig3D storage;
            storage.meshSize                   = values.meshSize;
            storage.decomposition              = values.parallelDimensions;
            storage.boundingBoxIncreasePercent = values.boundingBoxIncreasePercent;
            storage.periodicParticleBoundary   = std::all_of(
                    values.boundaryConditions.begin(), values.boundaryConditions.end(),
                    [](FieldBoundaryCondition kind) {
                        return kind == FieldBoundaryCondition::Periodic;
                    });
            if (values.backend == PoissonSolverType::P3M) {
                storage.layoutType    = ParticleLayoutType::SpatialOverlap;
                storage.overlapCutoff = values.p3mCutoff;
            }
            return storage;
        }

        TEST(SpaceChargeConfigTest, P3MSelectsOverlapLayoutAndBoundaryMode) {
            const auto openValues = p3mValues(FieldBoundaryCondition::Open);
            SpaceChargeConfig open(CartesianPICConfig(openValues), storageFor(openValues));
            const CartesianDomainConfig3D& openStorage = open.cartesianDomainConfig();
            EXPECT_EQ(openStorage.layoutType, ParticleLayoutType::SpatialOverlap);
            EXPECT_DOUBLE_EQ(openStorage.overlapCutoff, 0.025);
            EXPECT_FALSE(openStorage.periodicParticleBoundary);

            const auto periodicValues = p3mValues(FieldBoundaryCondition::Periodic);
            SpaceChargeConfig periodic(
                    CartesianPICConfig(periodicValues), storageFor(periodicValues));
            EXPECT_TRUE(periodic.cartesianDomainConfig().periodicParticleBoundary);
        }

        TEST(SpaceChargeConfigTest, P3MRejectsInvalidCombinations) {
            auto values      = p3mValues(FieldBoundaryCondition::Open);
            values.p3mCutoff = 0.0;
            EXPECT_THROW({ CartesianPICConfig config(values); }, OpalException);

            values                       = p3mValues(FieldBoundaryCondition::Open);
            values.boundaryConditions[1] = FieldBoundaryCondition::Periodic;
            EXPECT_THROW({ CartesianPICConfig config(values); }, OpalException);

            values = p3mValues(FieldBoundaryCondition::Open);
            values.binning.emplace(BinningConfig::Parameters{});
            EXPECT_THROW({ CartesianPICConfig config(values); }, OpalException);

            values            = p3mValues(FieldBoundaryCondition::Open);
            values.correction = CorrectionConfig(
                    {.kind = SpaceChargeCorrectionType::ImageCharge, .planeZ = 0.0});
            EXPECT_THROW({ CartesianPICConfig config(values); }, OpalException);
        }

        TEST(SpaceChargeConfigTest, FFT2D5UsesIndependentAlgorithmAndSpatialLayout) {
            FFT2D5Config::Parameters values;
            values.meshSize          = {16, 18, 20};
            values.pipeSizeX         = 0.1;
            values.pipeSizeY         = 0.2;
            values.beamRadius        = 0.01;
            values.referencePathFile = "design-path.dat";
            CartesianDomainConfig3D storage;
            storage.meshSize = values.meshSize;
            SpaceChargeConfig config{FFT2D5Config(values), storage};

            EXPECT_EQ(config.algorithmType(), SpaceChargeAlgorithmType::FFT2D5);
            EXPECT_EQ(config.cartesianDomainConfig().layoutType, ParticleLayoutType::Spatial);
            EXPECT_FALSE(config.cartesianDomainConfig().periodicParticleBoundary);
        }

        TEST(SpaceChargeConfigTest, RejectsInconsistentCartesianDomainConfig) {
            CartesianPICConfig::Parameters values;
            CartesianDomainConfig3D storage = storageFor(values);
            storage.meshSize[2] += 1;
            EXPECT_THROW(SpaceChargeConfig(CartesianPICConfig(values), storage), OpalException);
        }

        TEST(SpaceChargeConfigTest, ShiftedGreenRequiresOpenBackend) {
            CartesianPICConfig::Parameters values;
            values.backend    = PoissonSolverType::PeriodicFFT;
            values.correction = CorrectionConfig(
                    {.kind = SpaceChargeCorrectionType::ShiftedGreen, .planeZ = 0.0});
            EXPECT_THROW(static_cast<void>(CartesianPICConfig{values}), OpalException);
        }

        TEST(SpaceChargeRequestScheduleTest, ResolvesCorrectionScheduleAndBinning) {
            CartesianPICConfig::Parameters values;
            values.binning.emplace(BinningConfig::Parameters{});
            values.correction = CorrectionConfig(
                    {.kind               = SpaceChargeCorrectionType::ImageCharge,
                     .planeZ             = 0.25,
                     .planeDumpFrequency = 3,
                     .maximumSteps       = 2});
            const CartesianDomainConfig3D storage = storageFor(values);
            SpaceChargeConfig config(CartesianPICConfig(std::move(values)), storage);
            const SpaceChargeRequestSchedule schedule(config);

            const SpaceChargeRequest first = schedule.requestForStep(0);
            EXPECT_TRUE(first.useBinning);
            EXPECT_TRUE(first.writePotential);
            EXPECT_EQ(first.correction.kind, SpaceChargeCorrectionType::ImageCharge);
            EXPECT_DOUBLE_EQ(first.correction.planeZ, 0.25);

            const SpaceChargeRequest last = schedule.requestForStep(1);
            EXPECT_EQ(last.correction.kind, SpaceChargeCorrectionType::ImageCharge);

            const SpaceChargeRequest expired = schedule.requestForStep(2);
            EXPECT_TRUE(expired.useBinning);
            EXPECT_FALSE(expired.writePotential);
            EXPECT_EQ(expired.correction.kind, SpaceChargeCorrectionType::None);
            EXPECT_EQ(schedule.configuredCorrection().kind, SpaceChargeCorrectionType::ImageCharge);
            EXPECT_DOUBLE_EQ(schedule.configuredCorrection().planeZ, 0.25);
        }

        TEST(SpaceChargeConfigBuilderTest, RejectsRecognizedCGBeforeRuntimeConstruction) {
            TestableFieldSolverCmd command;
            command.setType("CG");
            EXPECT_THROW(
                    static_cast<void>(SpaceChargeConfigBuilder::build(command, {})), OpalException);
        }

    }  // namespace
}  // namespace opalx::spacecharge
