#include <gtest/gtest.h>
#include <limits>

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeConfigBuilder.h"
#include "Structure/BinningCmd.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/OpalException.h"

namespace opalx::spacecharge {
    namespace {

        class TestableFieldSolverCmd : public FieldSolverCmd {
        public:
            void setType(const std::string& type) {
                Attributes::setPredefinedString(itsAttr[FIELDSOLVER::TYPE], type);
            }
            void setBinsName(const std::string& name) {
                Attributes::setString(itsAttr[FIELDSOLVER::BINS], name);
            }
        };

        /** @brief Register a BINNING command @p name with MAXBINS = @p maximumBins. */
        void defineBinning(const std::string& name, int maximumBins) {
            BinningCmd prototype;
            BinningCmd* binning = prototype.clone(name);
            Attributes::setReal(*binning->findAttribute("MAXBINS"), maximumBins);
            OpalData::getInstance()->define(binning);
        }

        CartesianPIC3DConfig p3mConfig(FieldBoundaryCondition boundary) {
            CartesianPIC3DConfig config;
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

            config                = p3mConfig(FieldBoundaryCondition::Open);
            config.dirichletPlane = {.kind = DirichletPlaneType::ImageCharge};
            EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(config)), OpalException);
        }

        TEST(SpaceChargeConfigTest, DerivesIndependentFFT2D5Domain) {
            FFT2D5Config config;
            config.grid.meshSize       = {16, 18, 20};
            config.grid.decomposition  = {false, false, false};
            config.pipeSizeX           = 0.1;
            config.pipeSizeY           = 0.2;
            config.beamRadius          = 0.01;
            config.referencePathFile   = "design-path.dat";
            SpaceChargeConfig selected = config;
            validateSpaceChargeConfig(selected);

            EXPECT_TRUE(std::holds_alternative<FFT2D5Config>(selected));
            const auto domain = makeCartesianDomainConfig(selected);
            EXPECT_EQ(domain.meshSize, config.grid.meshSize);
            EXPECT_EQ(domain.decomposition, config.grid.decomposition);
            EXPECT_EQ(domain.layoutType, ParticleLayoutType::Spatial);
            EXPECT_FALSE(domain.periodicParticleBoundary);
        }

        TEST(SpaceChargeConfigTest, RejectsUnsupportedDirichletPlaneCombinations) {
            CartesianPIC3DConfig shifted;
            shifted.backend             = PoissonSolverType::Open;
            shifted.dirichletPlane.kind = DirichletPlaneType::ShiftedGreen;
            EXPECT_NO_THROW(validateSpaceChargeConfig(SpaceChargeConfig(shifted)));

            shifted.binning.emplace();
            EXPECT_NO_THROW(validateSpaceChargeConfig(SpaceChargeConfig(shifted)));

            CartesianPIC3DConfig binnedDump;
            binnedDump.backend = PoissonSolverType::Open;
            binnedDump.binning.emplace();
            binnedDump.dirichletPlane.kind               = DirichletPlaneType::ImageCharge;
            binnedDump.dirichletPlane.planeDumpFrequency = 1;
            EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(binnedDump)), OpalException);
        }

        TEST(SpaceChargeConfigTest, RejectsInvalidEnumValues) {
            CartesianPIC3DConfig cartesian;
            cartesian.backend = static_cast<PoissonSolverType>(255);
            EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(cartesian)), OpalException);

            cartesian                     = {};
            cartesian.dirichletPlane.kind = static_cast<DirichletPlaneType>(255);
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

        TEST(SpaceChargeConfigTest, ValidatesPoissonDomainBoundaryMatrix) {
            for (const auto backend :
                 {PoissonSolverType::None, PoissonSolverType::Open, PoissonSolverType::PeriodicFFT,
                  PoissonSolverType::P3M, PoissonSolverType::ConjugateGradient}) {
                for (const auto boundary :
                     {FieldBoundaryCondition::Open, FieldBoundaryCondition::Periodic,
                      FieldBoundaryCondition::Dirichlet}) {
                    PoissonSolverConfig config;
                    config.type      = backend;
                    config.p3mCutoff = backend == PoissonSolverType::P3M ? 0.1 : 0.0;
                    config.boundaryConditions.fill(boundary);
                    const bool accepted = boundary != FieldBoundaryCondition::Dirichlet
                                          && backend != PoissonSolverType::ConjugateGradient
                                          && (backend != PoissonSolverType::Open
                                              || boundary == FieldBoundaryCondition::Open)
                                          && (backend != PoissonSolverType::PeriodicFFT
                                              || boundary == FieldBoundaryCondition::Periodic);
                    if (accepted) {
                        EXPECT_NO_THROW(validatePoissonSolverConfig(config));
                    } else {
                        EXPECT_THROW(validatePoissonSolverConfig(config), OpalException);
                    }
                    config.boundaryConditions[0] = FieldBoundaryCondition::Dirichlet;
                    EXPECT_THROW(validatePoissonSolverConfig(config), OpalException);
                }
            }
        }

        TEST(SpaceChargeConfigTest, ShiftedGreenIsIndependentOfBinningAndKernelDiscretization) {
            CartesianPIC3DConfig config;
            config.backend             = PoissonSolverType::Open;
            config.dirichletPlane.kind = DirichletPlaneType::ShiftedGreen;
            for (auto green : {GreenFunctionType::Standard, GreenFunctionType::Integrated}) {
                config.greenFunction = green;
                config.binning.reset();
                EXPECT_NO_THROW(validateSpaceChargeConfig(config));
                config.binning.emplace();
                EXPECT_NO_THROW(validateSpaceChargeConfig(config));
            }
            config.boundaryConditions.fill(FieldBoundaryCondition::Periodic);
            EXPECT_THROW(validateSpaceChargeConfig(config), OpalException);
        }

        TEST(SpaceChargeConfigBuilderTest, ParserDefersCompatibilityAndReadsCurrentAttributes) {
            TestableFieldSolverCmd command;
            command.setType("P3M");
            command.setNX(8);
            command.setNY(8);
            command.setNZ(8);
            EXPECT_NO_THROW(command.execute());
            EXPECT_THROW((void)buildSpaceChargeConfig(command, {}), OpalException);
            command.setType("OPEN");
            EXPECT_EQ(command.getFieldSolverCmdType(), FieldSolverCmdType::OPEN);
            const auto snapshot =
                    std::get<CartesianPIC3DConfig>(buildSpaceChargeConfig(command, {}));
            command.setNX(16);
            EXPECT_EQ(snapshot.grid.meshSize[0], 8u);
            EXPECT_EQ(
                    std::get<CartesianPIC3DConfig>(buildSpaceChargeConfig(command, {}))
                            .grid.meshSize[0],
                    16u);
        }

        TEST(SpaceChargeConfigBuilderTest, RejectsInvalidMeshValuesBeforeIntegerConversion) {
            TestableFieldSolverCmd command;
            command.setType("OPEN");
            command.setNY(8);
            command.setNZ(8);
            for (double value :
                 {0.0, -1.0, 8.5, std::numeric_limits<double>::quiet_NaN(),
                  std::numeric_limits<double>::infinity(),
                  double(std::numeric_limits<int>::max()) + 1.0}) {
                command.setNX(value);
                EXPECT_THROW((void)buildSpaceChargeConfig(command, {}), OpalException);
            }
        }

        TEST(SpaceChargeConfigBuilderTest, BuildsBarnesHutWithoutMeshAndRejectsUnsupportedModes) {
            TestableFieldSolverCmd command;
            command.setType("BH");
            EXPECT_EQ(command.getFieldSolverCmdType(), FieldSolverCmdType::BH);
            // BHSOFTENING defaults to zero, which the self pair cannot evaluate.
            EXPECT_THROW((void)buildSpaceChargeConfig(command, {}), OpalException);
            command.setBHSoftening(1.0e-6);
            command.setBHTheta(0.3);
            const auto config = std::get<BarnesHutConfig>(buildSpaceChargeConfig(command, {}));
            EXPECT_DOUBLE_EQ(config.theta, 0.3);
            EXPECT_DOUBLE_EQ(config.softening, 1.0e-6);
            EXPECT_FALSE(config.leafBasedSoftening);
            EXPECT_FALSE(config.directSum);
            EXPECT_EQ(config.grid.meshSize, (std::array<std::size_t, 3>{8, 8, 8}));

            command.setBHDirectSum(true);
            EXPECT_TRUE(std::get<BarnesHutConfig>(buildSpaceChargeConfig(command, {})).directSum);

            command.setBHSoftening(0.0);
            command.setBHLeafBasedSoftening(true);
            // The direct sum has no octree leaves to soften with.
            EXPECT_THROW((void)buildSpaceChargeConfig(command, {}), OpalException);
            command.setBHDirectSum(false);
            EXPECT_NO_THROW((void)buildSpaceChargeConfig(command, {}));

            command.setBHTheta(1.5);
            EXPECT_THROW((void)buildSpaceChargeConfig(command, {}), OpalException);
            command.setBHTheta(0.5);
            command.setNX(8.5);
            EXPECT_THROW((void)buildSpaceChargeConfig(command, {}), OpalException);
        }

        TEST(SpaceChargeConfigBuilderTest, BuildsBarnesHutKinematicsFromBinsAndBHRELATIVISTIC) {
            // AUTO follows BINS: none is electrostatic, one bin is one rest-frame solve; more bins
            // would be energy binning, which TYPE=BH does not implement. FRAME overrides BINS.
            TestableFieldSolverCmd command;
            command.setType("BH");
            command.setBHSoftening(1.0e-6);
            const auto kinematics = [&command] {
                return std::get<BarnesHutConfig>(buildSpaceChargeConfig(command, {})).kinematics;
            };
            EXPECT_EQ(kinematics(), BarnesHutKinematics::Static);
            defineBinning("BH_ONE_BIN", 1);
            defineBinning("BH_TWO_BINS", 2);
            command.setBinsName("BH_ONE_BIN");
            EXPECT_EQ(kinematics(), BarnesHutKinematics::Frame);
            command.setBinsName("BH_TWO_BINS");
            EXPECT_THROW((void)buildSpaceChargeConfig(command, {}), OpalException);
            command.setBHRelativistic("FRAME");
            EXPECT_THROW((void)buildSpaceChargeConfig(command, {}), OpalException);
            command.setBinsName("NONE");
            EXPECT_EQ(kinematics(), BarnesHutKinematics::Frame);
            command.setBHRelativistic("AUTO");
            EXPECT_EQ(kinematics(), BarnesHutKinematics::Static);
        }

        TEST(SpaceChargeConfigBuilderTest, BuildsBarnesHutRelativeSoftening) {
            // BHSOFTENINGREL replaces BHSOFTENING (h in rest-frame spacings) and excludes BHLEAFH.
            TestableFieldSolverCmd command;
            command.setType("BH");
            command.setBHRelativeSoftening(0.5);
            const auto config = std::get<BarnesHutConfig>(buildSpaceChargeConfig(command, {}));
            EXPECT_DOUBLE_EQ(config.relativeSoftening, 0.5);
            EXPECT_DOUBLE_EQ(config.softening, 0.0);
            command.setBHDirectSum(true);
            EXPECT_NO_THROW((void)buildSpaceChargeConfig(command, {}));
            command.setBHDirectSum(false);
            command.setBHLeafBasedSoftening(true);
            EXPECT_THROW((void)buildSpaceChargeConfig(command, {}), OpalException);
            command.setBHLeafBasedSoftening(false);
            command.setBHRelativeSoftening(-1.0);
            EXPECT_THROW((void)buildSpaceChargeConfig(command, {}), OpalException);
        }

        TEST(SpaceChargeConfigBuilderTest, BuildsBarnesHutBoxMarginAndSlack) {
            TestableFieldSolverCmd command;
            command.setType("BH");
            command.setBHSoftening(1.0e-6);
            auto config = std::get<BarnesHutConfig>(buildSpaceChargeConfig(command, {}));
            EXPECT_DOUBLE_EQ(config.boxMargin, 0.05);  // defaults: the previous constants
            EXPECT_DOUBLE_EQ(config.boxSlack, 0.25);
            EXPECT_FALSE(config.traversalStatistics);
            command.setBHStatistics(true);
            EXPECT_TRUE(std::get<BarnesHutConfig>(buildSpaceChargeConfig(command, {})).traversalStatistics);
            EXPECT_EQ(config.precision, BarnesHutPrecision::Double);
            command.setBHPrecision("MIXED");
            EXPECT_EQ(std::get<BarnesHutConfig>(buildSpaceChargeConfig(command, {})).precision,
                      BarnesHutPrecision::Mixed);
            command.setBHPrecision("DOUBLE");
            command.setBHBoxMargin(0.2);
            command.setBHBoxSlack(1.0);
            config = std::get<BarnesHutConfig>(buildSpaceChargeConfig(command, {}));
            EXPECT_DOUBLE_EQ(config.boxMargin, 0.2);
            EXPECT_DOUBLE_EQ(config.boxSlack, 1.0);
            command.setBHBoxMargin(-0.1);
            EXPECT_THROW((void)buildSpaceChargeConfig(command, {}), OpalException);
            command.setBHBoxMargin(0.2);
            command.setBHBoxSlack(std::numeric_limits<double>::infinity());
            EXPECT_THROW((void)buildSpaceChargeConfig(command, {}), OpalException);
        }

        TEST(SpaceChargeConfigTest, ValidatesBarnesHutValues) {
            BarnesHutConfig config;
            config.softening = 1.0e-6;
            EXPECT_NO_THROW(validateSpaceChargeConfig(SpaceChargeConfig(config)));
            config.bucketSizeFocus = 0;
            EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(config)), OpalException);
            config.bucketSizeFocus = 64;
            config.softening       = -1.0;
            EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(config)), OpalException);
            config.softening = std::numeric_limits<double>::infinity();
            EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(config)), OpalException);
        }

    }  // namespace
}  // namespace opalx::spacecharge
