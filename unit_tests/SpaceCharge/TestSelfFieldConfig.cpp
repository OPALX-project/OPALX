#include <gtest/gtest.h>

#include "Attributes/Attributes.h"
#include "SpaceCharge/ParticleSetView.h"
#include "SpaceCharge/SelfFieldConfig.h"
#include "SpaceCharge/SelfFieldConfigBuilder.h"
#include "SpaceCharge/SelfFieldRequestPolicy.h"
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

        Pic3DConfigValues p3mValues(BoundaryConditionKind boundary) {
            Pic3DConfigValues values;
            values.backend            = PoissonBackendKind::P3M;
            values.p3mCutoff          = 0.025;
            values.boundaryConditions = {boundary, boundary, boundary};
            return values;
        }

        ParticleStorageConfig3d storageFor(const Pic3DConfigValues& values) {
            ParticleStorageConfig3d storage;
            storage.meshSize                   = values.meshSize;
            storage.decomposition              = values.parallelDimensions;
            storage.boundingBoxIncreasePercent = values.boundingBoxIncreasePercent;
            storage.periodicParticleBoundary   = std::all_of(
                    values.boundaryConditions.begin(), values.boundaryConditions.end(),
                    [](BoundaryConditionKind kind) {
                        return kind == BoundaryConditionKind::Periodic;
                    });
            if (values.backend == PoissonBackendKind::P3M) {
                storage.layoutKind    = ParticleLayoutKind::SpatialOverlap;
                storage.overlapCutoff = values.p3mCutoff;
            }
            return storage;
        }

        TEST(SelfFieldConfigTest, P3MSelectsOverlapLayoutAndBoundaryMode) {
            const auto openValues = p3mValues(BoundaryConditionKind::Open);
            SelfFieldConfig open(Pic3DConfig(openValues), storageFor(openValues));
            const ParticleStorageConfig3d& openStorage = open.particleStorageConfig();
            EXPECT_EQ(openStorage.layoutKind, ParticleLayoutKind::SpatialOverlap);
            EXPECT_DOUBLE_EQ(openStorage.overlapCutoff, 0.025);
            EXPECT_FALSE(openStorage.periodicParticleBoundary);

            const auto periodicValues = p3mValues(BoundaryConditionKind::Periodic);
            SelfFieldConfig periodic(Pic3DConfig(periodicValues), storageFor(periodicValues));
            EXPECT_TRUE(periodic.particleStorageConfig().periodicParticleBoundary);
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
            ParticleStorageConfig3d storage;
            storage.meshSize = values.meshSize;
            SelfFieldConfig config{Pic2d5Config(values), storage};

            EXPECT_EQ(config.algorithmKind(), SelfFieldAlgorithmKind::Pic2d5);
            EXPECT_EQ(config.particleStorageConfig().layoutKind, ParticleLayoutKind::Spatial);
            EXPECT_FALSE(config.particleStorageConfig().periodicParticleBoundary);
        }

        TEST(SelfFieldConfigTest, RejectsInconsistentParticleStorageSnapshot) {
            Pic3DConfigValues values;
            ParticleStorageConfig3d storage = storageFor(values);
            storage.meshSize[2] += 1;
            EXPECT_THROW(SelfFieldConfig(Pic3DConfig(values), storage), OpalException);
        }

        TEST(SelfFieldConfigTest, ShiftedGreenRequiresOpenBackend) {
            Pic3DConfigValues values;
            values.backend = PoissonBackendKind::FftPeriodic;
            values.correction =
                    CorrectionConfig({.kind = CorrectionKind::ShiftedGreen, .planeZ = 0.0});
            EXPECT_THROW(static_cast<void>(Pic3DConfig{values}), OpalException);
        }

        TEST(SelfFieldRequestPolicyTest, ResolvesCorrectionScheduleAndBinning) {
            Pic3DConfigValues values;
            values.binning.emplace(BinningConfigValues{});
            values.correction = CorrectionConfig(
                    {.kind               = CorrectionKind::ImageCharge,
                     .planeZ             = 0.25,
                     .planeDumpFrequency = 3,
                     .maximumSteps       = 2});
            const ParticleStorageConfig3d storage = storageFor(values);
            SelfFieldConfig config(Pic3DConfig(std::move(values)), storage);
            const SelfFieldRequestPolicy policy(config);

            const RequestedPhysics first = policy.forStep(0);
            EXPECT_TRUE(first.useBinning);
            EXPECT_TRUE(first.writePotential);
            EXPECT_EQ(first.correction.kind, CorrectionKind::ImageCharge);
            EXPECT_DOUBLE_EQ(first.correction.planeZ, 0.25);

            const RequestedPhysics last = policy.forStep(1);
            EXPECT_EQ(last.correction.kind, CorrectionKind::ImageCharge);

            const RequestedPhysics expired = policy.forStep(2);
            EXPECT_TRUE(expired.useBinning);
            EXPECT_FALSE(expired.writePotential);
            EXPECT_EQ(expired.correction.kind, CorrectionKind::None);
            EXPECT_EQ(policy.configuredCorrection().kind, CorrectionKind::ImageCharge);
            EXPECT_DOUBLE_EQ(policy.configuredCorrection().planeZ, 0.25);
        }

        TEST(SelfFieldConfigBuilderTest, RejectsRecognizedCGBeforeRuntimeConstruction) {
            TestableFieldSolverCmd command;
            command.setType("CG");
            EXPECT_THROW(
                    static_cast<void>(SelfFieldConfigBuilder::build(command, {})), OpalException);
        }

    }  // namespace
}  // namespace opalx::spacecharge
