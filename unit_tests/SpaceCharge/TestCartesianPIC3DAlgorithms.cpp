#include <gtest/gtest.h>

#include "AbstractObjects/OpalData.h"
#include "Algorithms/Quaternion.hpp"
#include "PartBunch/BunchStateHandler.h"
#include "SpaceCharge/CartesianPIC3D/CartesianPIC3DAlgorithm.h"
#include "Structure/DataSink.h"
#include "Utilities/Options.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>

namespace opalx::spacecharge {
    namespace {

        using Vector    = Vector_t<double, 3>;
        using Particles = ParticleContainer<double, 3>;

        class CartesianPIC3DAlgorithmsTest : public ::testing::Test {
        protected:
            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
                OpalData::getInstance()->storeInputFn("cartesian_pic3d_algorithms.opal");
                gmsg                = new Inform(nullptr, -1);
                Options::enableHDF5 = false;
            }

            static void TearDownTestSuite() {
                delete gmsg;
                gmsg = nullptr;
                std::remove("cartesian_pic3d_algorithms.stat");
                std::remove("cartesian_pic3d_algorithms.lbal");
                ippl::finalize();
            }

            struct Run {
                CartesianDomain<double, 3> domain;
                std::shared_ptr<BunchStateHandler> state = std::make_shared<BunchStateHandler>();
                Particles particles;
                Particles secondary;
                DataSink sink;
                std::unique_ptr<CartesianPIC3DAlgorithm> algorithm;
                std::array<std::uint8_t, 2> activity{1, 0};

                explicit Run(
                        const CartesianPIC3DConfig& config, const CoordinateSystemTrafo& pose = {},
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
                    algorithm = std::make_unique<CartesianPIC3DAlgorithm>(
                            config, containers,
                            std::make_unique<CartesianPIC3DFieldStorage<double, 3>>(domain), &sink,
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

            static CartesianPIC3DConfig config() {
                CartesianPIC3DConfig result;
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

        TEST_F(CartesianPIC3DAlgorithmsTest, ShiftedGreenWorksWithEitherKernelWithoutBinning) {
            for (auto green : {GreenFunctionType::Standard, GreenFunctionType::Integrated}) {
                auto wholeConfig           = config();
                wholeConfig.greenFunction  = green;
                wholeConfig.dirichletPlane = {.kind = DirichletPlaneType::ShiftedGreen};
                auto binnedConfig          = wholeConfig;
                binnedConfig.binning.emplace();
                binnedConfig.binning->maximumBins = 1;
                binnedConfig.binning->adaptive    = false;
                Run whole(wholeConfig), binned(binnedConfig);
                EXPECT_EQ(whole.solve().backendSolves, 2u);
                EXPECT_EQ(binned.solve().backendSolves, 2u);
                expectFieldsEqual(whole.particles, binned.particles);
                auto directConfig           = wholeConfig;
                directConfig.dirichletPlane = {};
                Run direct(directConfig);
                EXPECT_EQ(direct.solve().backendSolves, 1u);
                auto planeE = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), whole.particles.E.getView());
                auto directE = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), direct.particles.E.getView());
                double change = 0.0;
                for (std::size_t i = 0; i < whole.particles.getLocalNum(); ++i) {
                    change += std::abs(planeE(i)[2] - directE(i)[2]);
                }
                EXPECT_GT(change, 1.0e-6);
            }
        }

        TEST_F(CartesianPIC3DAlgorithmsTest, DirichletPlaneExpiryRestoresDirectFieldsAndMesh) {
            for (auto dirichletPlane :
                 {DirichletPlaneType::ShiftedGreen, DirichletPlaneType::ImageCharge}) {
                for (bool binned : {false, true}) {
                    auto planeConfig           = config();
                    planeConfig.dirichletPlane = {.kind = dirichletPlane, .maximumSteps = 2};
                    if (binned) {
                        planeConfig.binning.emplace();
                        planeConfig.binning->maximumBins = 1;
                        planeConfig.binning->adaptive    = false;
                    }
                    auto directConfig           = planeConfig;
                    directConfig.dirichletPlane = {};
                    Run withPlane(planeConfig), direct(directConfig);
                    EXPECT_EQ(
                            withPlane.solve(1).backendSolves,
                            dirichletPlane == DirichletPlaneType::ShiftedGreen || binned ? 2u : 1u);
                    EXPECT_EQ(
                            withPlane.domain.layoutExtents()[2],
                            dirichletPlane == DirichletPlaneType::ImageCharge ? 32u : 16u);
                    EXPECT_EQ(withPlane.solve(2).backendSolves, 1u);
                    EXPECT_EQ(withPlane.domain.layoutExtents()[2], 16u);
                    EXPECT_EQ(direct.solve(2).backendSolves, 1u);
                    expectFieldsEqual(withPlane.particles, direct.particles);
                    EXPECT_EQ(withPlane.solve(3).backendSolves, 1u);
                    expectFieldsEqual(withPlane.particles, direct.particles);
                }
            }
        }

        TEST_F(CartesianPIC3DAlgorithmsTest,
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

        TEST_F(CartesianPIC3DAlgorithmsTest,
               RotatedSolveRestoresMomentumAndRotatesFieldsConsistently) {
            auto values = config();
            values.binning.emplace();
            values.binning->maximumBins = 1;
            values.binning->adaptive    = false;
            const CoordinateSystemTrafo pose(
                    Vector(0.1, 0.2, 0.3), Quaternion(std::cos(0.3), 0.0, std::sin(0.3), 0.0));
            Run baseline(values, {}, Vector(0.0, 0.0, 1.0));
            Run rotated(values, pose, Vector(0.0, 0.0, 1.0));
            auto originalR = Kokkos::create_mirror(rotated.particles.R.getView());
            auto originalP = Kokkos::create_mirror(rotated.particles.P.getView());
            Kokkos::deep_copy(originalR, rotated.particles.R.getView());
            Kokkos::deep_copy(originalP, rotated.particles.P.getView());
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

        TEST_F(CartesianPIC3DAlgorithmsTest,
               DynamicSolveFrameTranslationPreservesFieldsAndDoesNotRequirePathHistory) {
            for (auto backend : {PoissonSolverType::None, PoissonSolverType::Open}) {
                SCOPED_TRACE(static_cast<int>(backend));
                auto values = config();
                values.backend = backend;
                values.repartitionFrequency = 0;
                values.binning.emplace();
                values.binning->maximumBins = 1;
                values.binning->adaptive = false;
                values.grid.boundingBoxIncreasePercent = 50.;
                // Same physical particles, with no image plane or fixed domain.
                // Only the solver coordinate origin differs between these runs.
                // pz=3/4 gives gamma=5/4, so the bin's longitudinal mesh
                // stretch also preserves the exactly representable geometry.
                Run centred(values, {}, Vector(0.0, 0.0, 0.75));
                Run translated(values, {}, Vector(0.0, 0.0, 0.75));
                // Each span is 15/1024 m. Fifty-percent padding on both sides
                // yields a 30/1024 m span and exactly binary spacing 2/1024 m
                // on 16 nodes. Particle positions and mesh origins then remain
                // exactly representable at both solver origins, isolating
                // translation invariance from the conditioning test below.
                const std::array<Vector, 4> initial{
                        Vector(-8., -7., 8.) / 1024., Vector(7., 8., 23.) / 1024.,
                        Vector(-4., 4., 13.) / 1024., Vector(3., -3., 18.) / 1024.};
                for (auto* run : {&centred, &translated}) {
                    auto positions = run->particles.R.getHostMirror();
                    for (std::size_t i = 0; i < initial.size(); ++i) positions(i) = initial[i];
                    Kokkos::deep_copy(run->particles.R.getView(), positions);
                    run->particles.markMomentsDirty();
                    run->particles.updateMoments();
                }
                const double pathTranslation = 825.;
                const CoordinateSystemTrafo solveToTracker(
                        Vector(0.0, 0.0, pathTranslation), Quaternion(1.0, 0.0, 0.0, 0.0));
                const auto originalR = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), centred.particles.R.getView());
                const auto originalP = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), centred.particles.P.getView());
                for (std::size_t step = 0; step < 2; ++step) {
                    SCOPED_TRACE(step);
                    const auto expectedSolves = backend == PoissonSolverType::Open ? 1u : 0u;
                    EXPECT_EQ(centred.solve(step).backendSolves, expectedSolves);
                    EXPECT_EQ(translated.solve(step, {solveToTracker.inverted(), solveToTracker})
                                      .backendSolves, expectedSolves);
                    expectFieldsEqual(translated.particles, centred.particles);
                    const auto centredR = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), centred.particles.R.getView());
                    const auto restoredR = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), translated.particles.R.getView());
                    const auto restoredP = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), translated.particles.P.getView());
                    for (std::size_t i = 0; i < centred.particles.getLocalNum(); ++i) {
                        for (unsigned d = 0; d < 3; ++d) {
                            EXPECT_DOUBLE_EQ(centredR(i)[d], originalR(i)[d]);
                            EXPECT_DOUBLE_EQ(restoredP(i)[d], originalP(i)[d]);
                            EXPECT_DOUBLE_EQ(restoredR(i)[d], originalR(i)[d]);
                        }
                    }
                }
            }
        }

        TEST_F(CartesianPIC3DAlgorithmsTest,
               CentredNoOpSolveAvoidsQuantizationFromLargePathTranslation) {
            auto values = config();
            values.backend = PoissonSolverType::None;
            values.repartitionFrequency = 0;
            Run centred(values, {}, Vector(0.0, 0.0, 0.073));
            Run translated(values, {}, Vector(0.0, 0.0, 0.073));
            auto originalR = centred.particles.R.getHostMirror();
            Kokkos::deep_copy(originalR, centred.particles.R.getView());
            for (std::size_t i = 0; i < centred.particles.getLocalNum(); ++i)
                originalR(i)[2] = (static_cast<double>(i) - 1.5) * 1e-12;
            for (auto* run : {&centred, &translated}) {
                Kokkos::deep_copy(run->particles.R.getView(), originalR);
                run->particles.markMomentsDirty();
                run->particles.updateMoments();
            }
            const double pathTranslation = 825.;
            const double coordinateSpacing = std::nextafter(pathTranslation,
                    std::numeric_limits<double>::infinity()) - pathTranslation;
            const CoordinateSystemTrafo solveToTracker(
                    Vector(0.0, 0.0, pathTranslation), Quaternion(1.0, 0.0, 0.0, 0.0));
            EXPECT_EQ(centred.solve().backendSolves, 0u);
            EXPECT_EQ(translated.solve(0, {solveToTracker.inverted(), solveToTracker}).backendSolves,
                      0u);
            const auto centredR = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), centred.particles.R.getView());
            const auto restoredR = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), translated.particles.R.getView());
            double maximumError = 0;
            for (std::size_t i = 0; i < centred.particles.getLocalNum(); ++i) {
                for (unsigned d = 0; d < 3; ++d)
                    EXPECT_DOUBLE_EQ(centredR(i)[d], originalR(i)[d]);
                maximumError = std::max(maximumError, std::abs(restoredR(i)[2] - originalR(i)[2]));
            }
            // Adding a picometre offset to 825 m quantizes it on that coordinate's
            // binary grid. The subsequent subtraction is exact (Sterbenz lemma).
            // This checks the known conditioning loss directly, with no field
            // comparison or relaxed field tolerance at that ill-conditioned scale.
            EXPECT_LE(maximumError, 0.5 * coordinateSpacing);
            EXPECT_GT(maximumError, coordinateSpacing / 16);
        }

        TEST_F(CartesianPIC3DAlgorithmsTest, TrivialP3MPrimaryPreservesSecondaryAndResumesSolving) {
            for (std::size_t primaryCount : {0u, 1u}) {
                SCOPED_TRACE(primaryCount);
                auto values               = config();
                values.backend            = PoissonSolverType::P3M;
                values.p3mCutoff          = 0.1;
                values.grid.meshSize      = {8, 8, 8};
                values.grid.decomposition = {false, false, true};
                CartesianDomain<double, 3> domain(makeCartesianDomainConfig(values));
                auto state = std::make_shared<BunchStateHandler>();
                Particles primary(
                        domain.mesh(), domain.layout(), false,
                        Particles::LayoutType::SpatialOverlap, values.p3mCutoff, ippl::BC::NO);
                Particles secondary(
                        domain.mesh(), domain.layout(), false,
                        Particles::LayoutType::SpatialOverlap, values.p3mCutoff, ippl::BC::NO);
                primary.setBunchStateHandler(state);
                secondary.setBunchStateHandler(state);
                const bool root = ippl::Comm->rank() == 0;
                primary.createParticles(root ? primaryCount : 0);
                secondary.createParticles(root ? 16 : 0);
                primary.setQ(1.0e-12);
                primary.setM(Physics::m_e);
                secondary.setQ(1.0e-12);
                secondary.setM(Physics::m_e);
                const Vector initialPosition(0.1, -0.1, 0.2);
                const Vector initialMomentum(0.2, 0.3, 0.4);
                Kokkos::deep_copy(primary.R.getView(), initialPosition);
                Kokkos::deep_copy(primary.P.getView(), initialMomentum);
                Kokkos::deep_copy(primary.dt.getView(), 1.0e-12);
                Kokkos::deep_copy(primary.E.getView(), Vector(123.0));
                Kokkos::deep_copy(primary.B.getView(), Vector(456.0));
                auto positions = secondary.R.getHostMirror();
                for (std::size_t i = 0; i < secondary.getLocalNum(); ++i) {
                    positions(i) =
                            Vector(i & 1 ? 1.0 : -1.0, i & 2 ? 1.0 : -1.0, i & 4 ? 1.0 : -1.0);
                }
                Kokkos::deep_copy(secondary.R.getView(), positions);
                Kokkos::deep_copy(secondary.P.getView(), Vector(0.0));
                Kokkos::deep_copy(secondary.E.getView(), Vector(7.0));
                Kokkos::deep_copy(secondary.B.getView(), Vector(9.0));
                DataSink sink;
                const std::array containers{&primary, &secondary};
                CartesianPIC3DAlgorithm algorithm(
                        values, containers,
                        std::make_unique<CartesianPIC3DFieldStorage<double, 3>>(domain), &sink,
                        state);
                const std::array<std::uint8_t, 2> activity{1, 1};
                SpaceChargeStepState step;
                step.timeStep = 1.0e-12;
                step.mpiSize  = ippl::Comm->size();
                const CoordinateSystemTrafo pose(
                        Vector(0.1, 0.2, 0.3), Quaternion(std::cos(0.3), 0.0, std::sin(0.3), 0.0));
                step.frames = {pose.inverted(), pose};
                const SpaceChargeSolveContext context(activity, step);

                EXPECT_EQ(algorithm.solve(context).backendSolves, 0u);
                EXPECT_EQ(primary.getTotalNum(), primaryCount);
                EXPECT_EQ(secondary.getTotalNum(), 16u);
                EXPECT_FALSE(primary.isMomentsDirty());
                EXPECT_FALSE(secondary.isMomentsDirty());
                const auto restoredR = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), primary.R.getView());
                const auto restoredP = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), primary.P.getView());
                const auto electric = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), primary.E.getView());
                const auto magnetic = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), primary.B.getView());
                const auto secondaryE = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), secondary.E.getView());
                const auto secondaryB = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), secondary.B.getView());
                for (unsigned d = 0; d < 3; ++d) {
                    EXPECT_LT(domain.lower()[d], -1.0);
                    EXPECT_GT(domain.upper()[d], 1.0);
                    for (std::size_t i = 0; i < primary.getLocalNum(); ++i) {
                        EXPECT_DOUBLE_EQ(restoredR(i)[d], initialPosition[d]);
                        EXPECT_DOUBLE_EQ(restoredP(i)[d], initialMomentum[d]);
                        EXPECT_DOUBLE_EQ(electric(i)[d], 0.0);
                        EXPECT_DOUBLE_EQ(magnetic(i)[d], 0.0);
                    }
                    for (std::size_t i = 0; i < secondary.getLocalNum(); ++i) {
                        EXPECT_DOUBLE_EQ(secondaryE(i)[d], 7.0);
                        EXPECT_DOUBLE_EQ(secondaryB(i)[d], 9.0);
                    }
                }

                primary.createParticles(root ? 8 : 0);
                primary.setQ(1.0e-12);
                primary.setM(Physics::m_e);
                auto grownPositions = primary.R.getHostMirror();
                for (std::size_t i = 0; i < primary.getLocalNum(); ++i) {
                    grownPositions(i) =
                            Vector(i & 1 ? 0.5 : -0.5, i & 2 ? 0.5 : -0.5, i & 4 ? 0.5 : -0.5);
                }
                Kokkos::deep_copy(primary.R.getView(), grownPositions);
                Kokkos::deep_copy(primary.P.getView(), Vector(0.0));
                Kokkos::deep_copy(primary.dt.getView(), 1.0e-12);
                EXPECT_EQ(algorithm.solve(context).backendSolves, 1u);
                EXPECT_EQ(primary.getTotalNum(), primaryCount + 8);
                EXPECT_EQ(secondary.getTotalNum(), 16u);
            }
        }

    }  // namespace
}  // namespace opalx::spacecharge
