#include <gtest/gtest.h>

#include "Attributes/Attributes.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeConfigBuilder.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/OpalException.h"

namespace opalx::spacecharge {
    namespace {

        class TestableFieldSolverCmd : public FieldSolverCmd {
        public:
            void setType(const std::string& type) {
                Attributes::setPredefinedString(itsAttr[FIELDSOLVER::TYPE], type);
            }
        };

        CartesianPICConfig p3mConfig(FieldBoundaryCondition boundary) {
            CartesianPICConfig config;
            config.backend            = PoissonSolverType::P3M;
            config.p3mCutoff          = 0.025;
            config.boundaryConditions = {boundary, boundary, boundary};
            return config;
        }

        TEST(SpaceChargeConfigTest, DerivesP3MParticleLayoutAndBoundaryMode) {
            SpaceChargeConfig open = p3mConfig(FieldBoundaryCondition::Open);
            validateSpaceChargeConfig(open);
            const CartesianDomainConfig3D openDomain = makeCartesianDomainConfig(open);
            EXPECT_EQ(openDomain.layoutType, ParticleLayoutType::SpatialOverlap);
            EXPECT_DOUBLE_EQ(openDomain.overlapCutoff, 0.025);
            EXPECT_FALSE(openDomain.periodicParticleBoundary);

            SpaceChargeConfig periodic = p3mConfig(FieldBoundaryCondition::Periodic);
            validateSpaceChargeConfig(periodic);
            EXPECT_TRUE(makeCartesianDomainConfig(periodic).periodicParticleBoundary);
        }

        TEST(SpaceChargeConfigTest, RejectsInvalidP3MCombinations) {
            auto config      = p3mConfig(FieldBoundaryCondition::Open);
            config.p3mCutoff = 0.0;
            EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(config)), OpalException);

            config                       = p3mConfig(FieldBoundaryCondition::Open);
            config.boundaryConditions[1] = FieldBoundaryCondition::Periodic;
            EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(config)), OpalException);

            config = p3mConfig(FieldBoundaryCondition::Open);
            config.binning.emplace();
            EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(config)), OpalException);

            config            = p3mConfig(FieldBoundaryCondition::Open);
            config.correction = {.kind = SpaceChargeCorrectionType::ImageCharge};
            EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(config)), OpalException);
        }

        TEST(SpaceChargeConfigTest, DerivesIndependentFFT2D5Domain) {
            FFT2D5Config config;
            config.grid.meshSize       = {16, 18, 20};
            config.grid.decomposition  = {false, true, false};
            config.pipeSizeX           = 0.1;
            config.pipeSizeY           = 0.2;
            config.beamRadius          = 0.01;
            config.referencePathFile   = "design-path.dat";
            SpaceChargeConfig selected = config;
            validateSpaceChargeConfig(selected);

            EXPECT_EQ(algorithmType(selected), SpaceChargeAlgorithmType::FFT2D5);
            const auto domain = makeCartesianDomainConfig(selected);
            EXPECT_EQ(domain.meshSize, config.grid.meshSize);
            EXPECT_EQ(domain.decomposition, config.grid.decomposition);
            EXPECT_EQ(domain.layoutType, ParticleLayoutType::Spatial);
            EXPECT_FALSE(domain.periodicParticleBoundary);
        }

        TEST(SpaceChargeConfigTest, RejectsUnsupportedCorrectionCombinations) {
            CartesianPICConfig shifted;
            shifted.backend         = PoissonSolverType::Open;
            shifted.correction.kind = SpaceChargeCorrectionType::ShiftedGreen;
            EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(shifted)), OpalException);

            shifted.binning.emplace();
            EXPECT_NO_THROW(validateSpaceChargeConfig(SpaceChargeConfig(shifted)));

            CartesianPICConfig binnedDump;
            binnedDump.backend = PoissonSolverType::Open;
            binnedDump.binning.emplace();
            binnedDump.correction.kind               = SpaceChargeCorrectionType::ImageCharge;
            binnedDump.correction.planeDumpFrequency = 1;
            EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(binnedDump)), OpalException);
        }

        TEST(SpaceChargeConfigTest, RejectsInvalidEnumValues) {
            CartesianPICConfig cartesian;
            cartesian.backend = static_cast<PoissonSolverType>(255);
            EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(cartesian)), OpalException);

            cartesian                 = {};
            cartesian.correction.kind = static_cast<SpaceChargeCorrectionType>(255);
            EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(cartesian)), OpalException);

            FFT2D5Config fft2d5;
            fft2d5.longitudinalFieldMode = static_cast<FFT2D5LongitudinalFieldMode>(255);
            EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(fft2d5)), OpalException);
        }

        TEST(SpaceChargeConfigBuilderTest, RejectsRecognizedCGBeforeRuntimeConstruction) {
            TestableFieldSolverCmd command;
            command.setType("CG");
            EXPECT_THROW(static_cast<void>(buildSpaceChargeConfig(command, {})), OpalException);
        }

    }  // namespace
}  // namespace opalx::spacecharge
