#include <gtest/gtest.h>

#include "AbstractObjects/OpalData.h"
#include "Algorithms/Quaternion.hpp"
#include "PartBunch/BunchStateHandler.h"
#include "SpaceCharge/CartesianPIC/CartesianPICAlgorithm.h"
#include "Structure/DataSink.h"
#include "Utilities/Options.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <memory>

namespace opalx::spacecharge {
    namespace {

        using Vector    = Vector_t<double, 3>;
        using Particles = ParticleContainer<double, 3>;

        class CartesianPICAlgorithmsTest : public ::testing::Test {
        protected:
            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
                OpalData::getInstance()->storeInputFn("cartesian_pic_algorithms.opal");
                gmsg                = new Inform(nullptr, -1);
                Options::enableHDF5 = false;
            }

            static void TearDownTestSuite() {
                delete gmsg;
                gmsg = nullptr;
                std::remove("cartesian_pic_algorithms.stat");
                std::remove("cartesian_pic_algorithms.lbal");
                ippl::finalize();
            }

            struct Run {
                CartesianDomain<double, 3> domain;
                std::shared_ptr<BunchStateHandler> state = std::make_shared<BunchStateHandler>();
                Particles particles;
                Particles secondary;
                DataSink sink;
                std::unique_ptr<CartesianPICAlgorithm> algorithm;
                std::array<std::uint8_t, 2> activity{1, 0};

                explicit Run(
                        const CartesianPICConfig& config, const CoordinateSystemTrafo& pose = {},
                        Vector momentum = Vector(0.0))
                    : domain(makeCartesianDomainConfig(config)),
                      particles(domain.mesh(), domain.layout()),
                      secondary(domain.mesh(), domain.layout()) {
                    particles.setBunchStateHandler(state);
                    secondary.setBunchStateHandler(state);
                    particles.createParticles(4);
                    particles.setQ(1.0e-12);
                    particles.setM(Physics::m_e);
                    auto positions = particles.R.getHostMirror();
                    auto momenta   = particles.P.getHostMirror();
                    auto dt        = particles.dt.getHostMirror();
                    const std::array<Vector, 4> initial{
                            Vector(-0.002, -0.001, 0.010), Vector(0.002, 0.001, 0.014),
                            Vector(-0.001, 0.002, 0.012), Vector(0.001, -0.002, 0.016)};
                    for (std::size_t i = 0; i < initial.size(); ++i) {
                        positions(i) = initial[i];
                        momenta(i)   = momentum;
                        dt(i)        = 1.0e-12;
                    }
                    Kokkos::deep_copy(particles.R.getView(), positions);
                    Kokkos::deep_copy(particles.P.getView(), momenta);
                    Kokkos::deep_copy(particles.dt.getView(), dt);
                    pose.transformBunchTo(particles.R.getView(), particles.getLocalNum());
                    pose.rotateBunchTo(particles.P.getView(), particles.getLocalNum());
                    particles.updateMoments();

                    secondary.createParticles(1);
                    secondary.setM(Physics::m_e);
                    Kokkos::deep_copy(secondary.R.getView(), Vector(0.0));
                    Kokkos::deep_copy(secondary.P.getView(), Vector(0.0));
                    Kokkos::deep_copy(secondary.E.getView(), Vector(7.0));
                    Kokkos::deep_copy(secondary.B.getView(), Vector(9.0));
                    const std::array containers{&particles, &secondary};
                    algorithm = std::make_unique<CartesianPICAlgorithm>(
                            config, containers,
                            std::make_unique<CartesianPICFieldStorage<double, 3>>(domain), &sink,
                            state);
                }

                SpaceChargeSolveResult solve(
                        std::size_t step = 0, CoordinateFrameTransforms frames = {}) {
                    SpaceChargeStepState state;
                    state.step     = step;
                    state.timeStep = 1.0e-12;
                    state.frames   = frames;
                    return algorithm->solve(SpaceChargeSolveContext(activity, state));
                }
            };

            static CartesianPICConfig config() {
                CartesianPICConfig result;
                result.backend            = PoissonSolverType::Open;
                result.grid.meshSize      = {16, 16, 16};
                result.grid.decomposition = {false, false, false};
                return result;
            }

            static void expectFieldsEqual(Particles& actual, Particles& expected) {
                auto actualE = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), actual.E.getView());
                auto expectedE = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), expected.E.getView());
                auto actualB = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), actual.B.getView());
                auto expectedB = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), expected.B.getView());
                // These serial tests retain particle order through migration and their single-bin
                // hash.
                for (std::size_t i = 0; i < actual.getLocalNum(); ++i) {
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_TRUE(std::isfinite(actualE(i)[d]));
                        EXPECT_NEAR(
                                actualE(i)[d], expectedE(i)[d],
                                1.0e-9 * std::max(1.0, std::abs(expectedE(i)[d])));
                        EXPECT_NEAR(actualB(i)[d], expectedB(i)[d], 1.0e-15);
                    }
                }
            }
        };

        TEST_F(CartesianPICAlgorithmsTest, ShiftedGreenWorksWithEitherKernelWithoutBinning) {
            for (auto green : {GreenFunctionType::Standard, GreenFunctionType::Integrated}) {
                auto wholeConfig          = config();
                wholeConfig.greenFunction = green;
                wholeConfig.correction    = {.kind = SpaceChargeCorrectionType::ShiftedGreen};
                auto binnedConfig         = wholeConfig;
                binnedConfig.binning.emplace();
                binnedConfig.binning->maximumBins = 1;
                binnedConfig.binning->adaptive    = false;
                Run whole(wholeConfig), binned(binnedConfig);
                EXPECT_EQ(whole.solve().backendSolves, 2u);
                EXPECT_EQ(binned.solve().backendSolves, 2u);
                expectFieldsEqual(whole.particles, binned.particles);
                auto directConfig       = wholeConfig;
                directConfig.correction = {};
                Run direct(directConfig);
                EXPECT_EQ(direct.solve().backendSolves, 1u);
                auto correctedE = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), whole.particles.E.getView());
                auto directE = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), direct.particles.E.getView());
                double change = 0.0;
                for (std::size_t i = 0; i < whole.particles.getLocalNum(); ++i) {
                    change += std::abs(correctedE(i)[2] - directE(i)[2]);
                }
                EXPECT_GT(change, 1.0e-6);
            }
        }

        TEST_F(CartesianPICAlgorithmsTest, CorrectionExpiryRestoresDirectFieldsAndMesh) {
            for (auto correction :
                 {SpaceChargeCorrectionType::ShiftedGreen,
                  SpaceChargeCorrectionType::ImageCharge}) {
                for (bool binned : {false, true}) {
                    auto correctedConfig       = config();
                    correctedConfig.correction = {.kind = correction, .maximumSteps = 2};
                    if (binned) {
                        correctedConfig.binning.emplace();
                        correctedConfig.binning->maximumBins = 1;
                        correctedConfig.binning->adaptive    = false;
                    }
                    auto directConfig       = correctedConfig;
                    directConfig.correction = {};
                    Run corrected(correctedConfig), direct(directConfig);
                    EXPECT_EQ(
                            corrected.solve(1).backendSolves,
                            correction == SpaceChargeCorrectionType::ShiftedGreen || binned ? 2u
                                                                                            : 1u);
                    EXPECT_EQ(
                            corrected.domain.layoutExtents()[2],
                            correction == SpaceChargeCorrectionType::ImageCharge ? 32u : 16u);
                    EXPECT_EQ(corrected.solve(2).backendSolves, 1u);
                    EXPECT_EQ(corrected.domain.layoutExtents()[2], 16u);
                    EXPECT_EQ(direct.solve(2).backendSolves, 1u);
                    expectFieldsEqual(corrected.particles, direct.particles);
                    EXPECT_EQ(corrected.solve(3).backendSolves, 1u);
                    expectFieldsEqual(corrected.particles, direct.particles);
                }
            }
        }

        TEST_F(CartesianPICAlgorithmsTest,
               ReplacesFieldsUsingExistingAllocationsAndPreservesOtherContainers) {
            for (auto backend : {PoissonSolverType::None, PoissonSolverType::Open}) {
                auto values    = config();
                values.backend = backend;
                Run run(values), clean(values);
                auto* originalE = run.particles.E.getView().data();
                auto* originalB = run.particles.B.getView().data();
                Kokkos::deep_copy(run.particles.E.getView(), Vector(123.0));
                Kokkos::deep_copy(run.particles.B.getView(), Vector(456.0));
                (void)run.solve();
                (void)clean.solve();
                expectFieldsEqual(run.particles, clean.particles);
                EXPECT_EQ(run.particles.E.getView().data(), originalE);
                EXPECT_EQ(run.particles.B.getView().data(), originalB);
                auto secondaryE = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), run.secondary.E.getView());
                auto secondaryB = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), run.secondary.B.getView());
                for (unsigned d = 0; d < 3; ++d) {
                    EXPECT_DOUBLE_EQ(secondaryE(0)[d], 7.0);
                    EXPECT_DOUBLE_EQ(secondaryB(0)[d], 9.0);
                }
            }
        }

        TEST_F(CartesianPICAlgorithmsTest,
               RotatedSolveRestoresMomentumAndRotatesFieldsConsistently) {
            auto values = config();
            values.binning.emplace();
            values.binning->maximumBins = 1;
            values.binning->adaptive    = false;
            const CoordinateSystemTrafo pose(
                    Vector(0.1, 0.2, 0.3), Quaternion(std::cos(0.3), 0.0, std::sin(0.3), 0.0));
            Run baseline(values, {}, Vector(0.0, 0.0, 1.0));
            Run rotated(values, pose, Vector(0.0, 0.0, 1.0));
            auto originalR = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), rotated.particles.R.getView());
            auto originalP = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), rotated.particles.P.getView());
            (void)baseline.solve();
            (void)rotated.solve(0, {pose.inverted(), pose});
            auto restoredR = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), rotated.particles.R.getView());
            auto restoredP = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), rotated.particles.P.getView());
            for (std::size_t i = 0; i < rotated.particles.getLocalNum(); ++i) {
                for (unsigned d = 0; d < 3; ++d) {
                    EXPECT_NEAR(restoredR(i)[d], originalR(i)[d], 1.0e-14);
                    EXPECT_NEAR(restoredP(i)[d], originalP(i)[d], 1.0e-14);
                }
            }
            pose.rotateBunchTo(baseline.particles.E.getView(), baseline.particles.getLocalNum());
            pose.rotateBunchTo(baseline.particles.B.getView(), baseline.particles.getLocalNum());
            expectFieldsEqual(rotated.particles, baseline.particles);
        }

    }  // namespace
}  // namespace opalx::spacecharge
