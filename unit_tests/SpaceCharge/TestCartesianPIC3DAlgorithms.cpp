/**
 * @file TestCartesianPIC3DAlgorithms.cpp
 * @brief End-to-end tests of `CartesianPIC3DAlgorithm::solve()` with the NONE, OPEN, and P3M
 * backends.
 *
 * Each test runs complete space-charge solves for a small bunch of four electrons on a 16^3 mesh
 * and checks the particle E and B fields that come out. Covers:
 * - Dirichlet planes: a shifted-Green plane gives the same fields with and without a single bin
 *   and with either Green function; after `ZEROFACE_MAXSTEPS`, both plane types return to the
 *   plain solve and mesh.
 * - P3M with a shifted-Green plane: the image field equals the STANDARD OPEN image, and the
 *   short-range correction is added once, to the real bunch only.
 * - Binned P3M: with one bin at rest it reproduces the unbinned P3M solve; a moving pair follows
 *   the Lorentz-boosted Coulomb field with B = v x E / c^2. With two bins, a close pair across
 *   the bins follows the Coulomb field of each source boosted with its own velocity. The
 *   shifted-Green image of a bunch with two energies equals the OPEN image with one and with two
 *   bins.
 * - A mesh-tied P3M cutoff (`RCUT_CELLS`) gives the same fields as the equivalent fixed cutoff,
 *   also after the mesh spacing changes.
 * - Solve frames: rotated and translated frames give the same fields and restore R and P; a large
 *   path translation changes positions only by the expected round-off.
 * - Bookkeeping: E and B keep their allocations, other containers are left untouched, and a P3M
 *   solve with fewer than two primary particles is skipped without blocking later solves.
 *
 * The fixture builds a primary and a secondary particle container, the field storage, and the
 * algorithm on one rank. P3M sorts particles by cell, so the binned P3M tests compare runs
 * particle by particle at matching positions.
 */

#include <gtest/gtest.h>

#include "AbstractObjects/OpalData.h"
#include "Algorithms/Quaternion.hpp"
#include "PartBunch/BunchStateHandler.h"
#include "SpaceCharge/CartesianPIC3D/CartesianPIC3DAlgorithm.h"
#include "Structure/DataSink.h"
#include "Utilities/Options.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

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
                      particles(
                              domain.mesh(), domain.layout(), false,
                              config.backend == PoissonSolverType::P3M
                                      ? Particles::LayoutType::SpatialOverlap
                                      : Particles::LayoutType::Spatial,
                              resolveOverlapCutoff<3>(
                                      config.p3mCutoff, config.p3mCutoffCells, domain.spacing()),
                              config.backend == PoissonSolverType::P3M ? ippl::BC::NO
                                                                       : ippl::BC::PERIODIC),
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
                // Spatial layouts on one rank keep the particle order through migration and the
                // single-bin hash. Compare overlap-layout (P3M) runs with
                // expectFieldsEqualByPosition() instead.
                for (std::size_t i = 0; i < actual.getLocalNum(); ++i) {
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_TRUE(std::isfinite(actualE(i)[d]));
                        EXPECT_NEAR(
                                actualE(i)[d], expectedE(i)[d],
                                1.0e-9 * std::max(1.0, std::abs(expectedE(i)[d])));
                        EXPECT_NEAR(
                                actualB(i)[d], expectedB(i)[d],
                                1.0e-9 * std::max(1.0, std::abs(expectedB(i)[d])));
                    }
                }
            }

            static double largestComponent(const Kokkos::View<Vector*, Kokkos::HostSpace>& view) {
                double largest = 0.0;
                for (std::size_t i = 0; i < view.extent(0); ++i) {
                    for (unsigned d = 0; d < 3; ++d) {
                        largest = std::max(largest, std::abs(view(i)[d]));
                    }
                }
                return largest;
            }

            template <typename Attribute>
            static Kokkos::View<Vector*, Kokkos::HostSpace> localHostCopy(
                    const Attribute& attribute, std::size_t count) {
                Kokkos::View<Vector*, Kokkos::HostSpace> result("host", count);
                Kokkos::deep_copy(
                        result,
                        Kokkos::subview(
                                attribute.getView(), std::make_pair(std::size_t{0}, count)));
                return result;
            }

            /**
             * @brief Index of the particle of @p expected at the position of each local particle
             * of @p actual.
             *
             * Overlap layouts sort particles by cell. A threaded or GPU sort can order particles
             * that share a cell differently in each run, so P3M runs are paired by position.
             * @return Empty, with a recorded failure, if the particle counts or positions differ.
             */
            static std::vector<std::size_t> matchByPosition(
                    Particles& actual, Particles& expected) {
                const std::size_t count = actual.getLocalNum();
                if (expected.getLocalNum() != count) {
                    ADD_FAILURE() << "The runs hold " << count << " and " << expected.getLocalNum()
                                  << " particles.";
                    return {};
                }
                const auto actualR   = localHostCopy(actual.R, count);
                const auto expectedR = localHostCopy(expected.R, count);
                std::vector<std::size_t> match(count, count);
                for (std::size_t i = 0; i < count; ++i) {
                    for (std::size_t j = 0; j < count; ++j) {
                        const Vector offset = actualR(i) - expectedR(j);
                        if (std::sqrt(offset.dot(offset)) < 1.0e-12) {
                            match[i] = j;
                        }
                    }
                    if (match[i] == count) {
                        ADD_FAILURE() << "No particle at the position of particle " << i;
                        return {};
                    }
                }
                return match;
            }

            static void expectFieldsEqualByPosition(Particles& actual, Particles& expected) {
                const auto match = matchByPosition(actual, expected);
                ASSERT_EQ(match.size(), actual.getLocalNum());
                const std::size_t count    = match.size();
                const auto actualE         = localHostCopy(actual.E, count),
                           expectedE       = localHostCopy(expected.E, count);
                const auto actualB         = localHostCopy(actual.B, count),
                           expectedB       = localHostCopy(expected.B, count);
                const double magneticScale = largestComponent(expectedB);
                for (std::size_t i = 0; i < count; ++i) {
                    const std::size_t j = match[i];
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_TRUE(std::isfinite(actualE(i)[d]));
                        EXPECT_NEAR(
                                actualE(i)[d], expectedE(j)[d],
                                1.0e-9 * std::max(1.0, std::abs(expectedE(j)[d])));
                        EXPECT_NEAR(
                                actualB(i)[d], expectedB(j)[d],
                                1.0e-9 * std::max(std::abs(expectedB(j)[d]), magneticScale));
                    }
                }
            }

            // Moves particle 1 within the P3M cutoff of particle 0, then shifts the bunch in z.
            static void placeClosePair(Particles& particles, double zShift) {
                auto positions = particles.R.getHostMirror();
                Kokkos::deep_copy(positions, particles.R.getView());
                positions(1) = positions(0) + Vector(0.0002, 0.0001, 0.0001);
                for (std::size_t i = 0; i < particles.getLocalNum(); ++i)
                    positions(i)[2] += zShift;
                Kokkos::deep_copy(particles.R.getView(), positions);
                particles.updateMoments();
            }

            static CartesianPIC3DConfig binnedP3MConfig(std::size_t maximumBins) {
                auto result      = config();
                result.backend   = PoissonSolverType::P3M;
                result.p3mCutoff = 0.001;
                result.binning.emplace();
                result.binning->maximumBins = maximumBins;
                result.binning->adaptive    = false;
                return result;
            }

            static void setLongitudinalMomenta(
                    Particles& particles, const std::array<double, 4>& momenta) {
                auto momentum = particles.P.getHostMirror();
                for (std::size_t i = 0; i < momenta.size(); ++i)
                    momentum(i) = Vector(0.0, 0.0, momenta[i]);
                Kokkos::deep_copy(particles.P.getView(), momentum);
                particles.markMomentsDirty();
                particles.updateMoments();
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

        TEST_F(CartesianPIC3DAlgorithmsTest, P3MAddsFullImageAndRealParticleCorrectionOnce) {
            auto p3mConfig           = config();
            p3mConfig.backend        = PoissonSolverType::P3M;
            p3mConfig.p3mCutoff      = 0.001;
            auto openConfig          = config();
            openConfig.greenFunction = GreenFunctionType::Standard;
            for (double planeZ : {0.0, 0.125}) {
                p3mConfig.dirichletPlane = {
                        .kind         = DirichletPlaneType::ShiftedGreen,
                        .planeZ       = planeZ,
                        .maximumSteps = 2};
                openConfig.dirichletPlane = p3mConfig.dirichletPlane;
                Run p3mImage(p3mConfig), openImage(openConfig);
                auto p3mDirectConfig            = p3mConfig;
                auto openDirectConfig           = openConfig;
                p3mDirectConfig.dirichletPlane  = {};
                openDirectConfig.dirichletPlane = {};
                Run p3mDirect(p3mDirectConfig), openDirect(openDirectConfig);
                for (Run* run : {&p3mImage, &p3mDirect, &openImage, &openDirect}) {
                    auto positions = run->particles.R.getHostMirror();
                    Kokkos::deep_copy(positions, run->particles.R.getView());
                    positions(1) = positions(0) + Vector(0.0002, 0.0001, 0.0001);
                    for (std::size_t i = 0; i < run->particles.getLocalNum(); ++i)
                        positions(i)[2] += planeZ;
                    Kokkos::deep_copy(run->particles.R.getView(), positions);
                    run->particles.updateMoments();
                }

                for (std::size_t step : {0u, 1u}) {
                    EXPECT_EQ(p3mImage.solve(step).backendSolves, 2u);
                    EXPECT_EQ(openImage.solve(step).backendSolves, 2u);
                    EXPECT_EQ(p3mDirect.solve(step).backendSolves, 1u);
                    EXPECT_EQ(openDirect.solve(step).backendSolves, 1u);
                    EXPECT_EQ(p3mImage.domain.layoutExtents(), p3mConfig.grid.meshSize);
                    auto p3mE = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), p3mImage.particles.E.getView());
                    auto p3mDirectE = Kokkos::create_mirror(p3mDirect.particles.E.getView());
                    Kokkos::deep_copy(p3mDirectE, p3mDirect.particles.E.getView());
                    auto openE = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), openImage.particles.E.getView());
                    auto openDirectE = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), openDirect.particles.E.getView());
                    const std::size_t count = p3mImage.particles.getLocalNum();
                    const auto toP3MDirect =
                            matchByPosition(p3mImage.particles, p3mDirect.particles);
                    const auto toOpen = matchByPosition(p3mImage.particles, openImage.particles);
                    const auto toOpenDirect =
                            matchByPosition(p3mImage.particles, openDirect.particles);
                    ASSERT_EQ(toP3MDirect.size(), count);
                    ASSERT_EQ(toOpen.size(), count);
                    ASSERT_EQ(toOpenDirect.size(), count);
                    double imageNorm = 0.0;
                    for (std::size_t i = 0; i < count; ++i) {
                        for (unsigned d = 0; d < 3; ++d) {
                            const double image =
                                    openE(toOpen[i])[d] - openDirectE(toOpenDirect[i])[d];
                            imageNorm += image * image;
                            EXPECT_NEAR(
                                    p3mE(i)[d] - p3mDirectE(toP3MDirect[i])[d], image,
                                    1.0e-9 * std::max(1.0, std::abs(image)));
                        }
                    }
                    EXPECT_GT(imageNorm, 1.0e-6);

                    // Ensure this fixture would detect a lost or repeated particle correction.
                    Kokkos::deep_copy(p3mDirect.particles.E.getView(), Vector(0.0));
                    P3MShortRangeInteraction().apply(p3mDirect.particles);
                    auto particleE = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), p3mDirect.particles.E.getView());
                    double particleNorm = 0.0;
                    for (std::size_t i = 0; i < p3mDirect.particles.getLocalNum(); ++i)
                        particleNorm += particleE(i).dot(particleE(i));
                    EXPECT_GT(particleNorm, 1.0e-6);
                    Kokkos::deep_copy(p3mDirect.particles.E.getView(), p3mDirectE);
                }
                for (std::size_t step : {2u, 3u}) {
                    EXPECT_EQ(p3mImage.solve(step).backendSolves, 1u);
                    EXPECT_EQ(p3mDirect.solve(step).backendSolves, 1u);
                    expectFieldsEqualByPosition(p3mImage.particles, p3mDirect.particles);
                }
            }
        }

        TEST_F(CartesianPIC3DAlgorithmsTest, P3MSingleBinAtRestMatchesUnbinnedP3M) {
            const auto binnedConfig = binnedP3MConfig(1);
            auto wholeConfig        = binnedConfig;
            wholeConfig.binning.reset();
            Run whole(wholeConfig), binned(binnedConfig);
            for (Run* run : {&whole, &binned})
                placeClosePair(run->particles, 0.0);

            for (std::size_t step : {0u, 1u}) {
                EXPECT_EQ(whole.solve(step).backendSolves, 1u);
                EXPECT_EQ(binned.solve(step).backendSolves, 1u);
                expectFieldsEqualByPosition(binned.particles, whole.particles);
            }

            // Ensure this fixture would detect a lost or repeated particle correction.
            Kokkos::deep_copy(whole.particles.E.getView(), Vector(0.0));
            P3MShortRangeInteraction().apply(whole.particles);
            auto particleE = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), whole.particles.E.getView());
            double particleNorm = 0.0;
            for (std::size_t i = 0; i < whole.particles.getLocalNum(); ++i)
                particleNorm += particleE(i).dot(particleE(i));
            EXPECT_GT(particleNorm, 1.0e-6);
        }

        TEST_F(CartesianPIC3DAlgorithmsTest, P3MSingleBinFollowsBoostedCoulomb) {
            auto values          = binnedP3MConfig(1);
            values.grid.meshSize = {32, 32, 32};
            // pz = 4/3 gives gamma = 5/3 and v = 0.8 c along z.
            const double gamma = 5.0 / 3.0;
            Run run(values, {}, Vector(0.0, 0.0, 4.0 / 3.0));

            // A close pair, dominated by its short-range term, and two distant particles that
            // interact through the mesh only. The rest-frame mesh keeps alpha * h near 0.5.
            const std::array<Vector, 4> placed{
                    Vector(0.0, 0.0, 0.010), Vector(0.0001, 0.00005, 0.0101),
                    Vector(-0.002, 0.0015, 0.008), Vector(0.0018, -0.002, 0.0125)};
            auto positions = run.particles.R.getHostMirror();
            for (std::size_t i = 0; i < placed.size(); ++i)
                positions(i) = placed[i];
            Kokkos::deep_copy(run.particles.R.getView(), positions);
            run.particles.markMomentsDirty();
            run.particles.updateMoments();

            EXPECT_EQ(run.solve().backendSolves, 1u);

            const auto finalR = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), run.particles.R.getView());
            const auto electric = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), run.particles.E.getView());
            const auto magnetic = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), run.particles.B.getView());
            const double coulomb =
                    run.particles.getChargePerParticle() / (4.0 * Physics::pi * Physics::epsilon_0);
            const Vector velocity(0.0, 0.0, 0.8 * Physics::c);
            ASSERT_EQ(run.particles.getLocalNum(), placed.size());
            for (std::size_t i = 0; i < placed.size(); ++i) {
                // Coulomb field of every other particle in the bunch rest frame, boosted to lab.
                Vector expected(0.0);
                for (const Vector& source : placed) {
                    Vector rest = finalR(i) - source;
                    rest[2] *= gamma;
                    const double r = std::sqrt(rest.dot(rest));
                    if (r < 1.0e-9) {
                        continue;
                    }
                    const Vector restField = (coulomb / (r * r * r)) * rest;
                    expected += Vector(gamma * restField[0], gamma * restField[1], restField[2]);
                }
                const double scale = std::sqrt(expected.dot(expected));
                const Vector expectedB =
                        (1.0 / (Physics::c * Physics::c)) * cross(velocity, Vector(electric(i)));
                const double scaleB = std::sqrt(expectedB.dot(expectedB));
                EXPECT_GT(scaleB, 0.0);
                for (unsigned d = 0; d < 3; ++d) {
                    EXPECT_NEAR(electric(i)[d], expected[d], 0.03 * scale)
                            << "particle " << i << ", component " << d;
                    EXPECT_NEAR(magnetic(i)[d], expectedB[d], 1.0e-9 * scaleB);
                }
            }
        }

        TEST_F(CartesianPIC3DAlgorithmsTest, P3MTwoBinsFollowPerSourceBoostedCoulomb) {
            auto values          = binnedP3MConfig(2);
            values.grid.meshSize = {32, 32, 32};
            Run run(values);

            // The close pair 0-1 spans both bins: pz = 3/4 and 4/3 give gamma = 5/4 and 5/3.
            // Each bin holds a single momentum, so the mesh has no velocity mismatch and the
            // near field tests that each source enters with its own gamma and velocity.
            const std::array<Vector, 4> placed{
                    Vector(0.0, 0.0, 0.010), Vector(0.0001, 0.00005, 0.0101),
                    Vector(-0.002, 0.0015, 0.008), Vector(0.0018, -0.002, 0.0125)};
            auto positions = run.particles.R.getHostMirror();
            for (std::size_t i = 0; i < placed.size(); ++i)
                positions(i) = placed[i];
            Kokkos::deep_copy(run.particles.R.getView(), positions);
            setLongitudinalMomenta(run.particles, {0.75, 4.0 / 3.0, 0.75, 4.0 / 3.0});

            EXPECT_EQ(run.solve().backendSolves, 2u);

            auto host = [](const auto& attribute) {
                return Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), attribute.getView());
            };
            const auto finalR = host(run.particles.R), finalP = host(run.particles.P);
            const auto electric = host(run.particles.E), magnetic = host(run.particles.B);
            const double coulomb =
                    run.particles.getChargePerParticle() / (4.0 * Physics::pi * Physics::epsilon_0);
            ASSERT_EQ(run.particles.getLocalNum(), placed.size());
            for (std::size_t i = 0; i < placed.size(); ++i) {
                // Coulomb field of every other particle in its own rest frame, boosted to lab.
                Vector expectedE(0.0), expectedB(0.0);
                for (std::size_t j = 0; j < placed.size(); ++j) {
                    if (j == i) {
                        continue;
                    }
                    const double gamma = std::sqrt(1.0 + finalP(j).dot(finalP(j)));
                    Vector rest        = finalR(i) - finalR(j);
                    rest[2] *= gamma;
                    const double r         = std::sqrt(rest.dot(rest));
                    const Vector restField = (coulomb / (r * r * r)) * rest;
                    const Vector field(gamma * restField[0], gamma * restField[1], restField[2]);
                    const Vector velocity = (Physics::c / gamma) * Vector(finalP(j));
                    expectedE += field;
                    expectedB += (1.0 / (Physics::c * Physics::c)) * cross(velocity, field);
                }
                const double scaleE = std::sqrt(expectedE.dot(expectedE));
                const double scaleB = std::sqrt(expectedB.dot(expectedB));
                EXPECT_GT(scaleB, 0.0);
                for (unsigned d = 0; d < 3; ++d) {
                    EXPECT_NEAR(electric(i)[d], expectedE[d], 0.03 * scaleE)
                            << "particle " << i << ", component " << d;
                    EXPECT_NEAR(magnetic(i)[d], expectedB[d], 0.03 * scaleB)
                            << "particle " << i << ", component " << d;
                }
            }
        }

        TEST_F(CartesianPIC3DAlgorithmsTest, P3MBinnedShiftedImageMatchesOpenImage) {
            // pz = 3/4 and 4/3 give gamma = 5/4 and 5/3, so the image shift and field conversion
            // are boosted; with two bins each bin gets its own boost and image pass.
            const std::array<double, 4> momenta{0.75, 4.0 / 3.0, 0.75, 4.0 / 3.0};
            for (std::size_t bins : {1u, 2u}) {
                auto p3mConfig           = binnedP3MConfig(bins);
                auto openConfig          = config();
                openConfig.greenFunction = GreenFunctionType::Standard;
                openConfig.binning       = p3mConfig.binning;
                for (double planeZ : {0.0, 0.125}) {
                    p3mConfig.dirichletPlane = {
                            .kind         = DirichletPlaneType::ShiftedGreen,
                            .planeZ       = planeZ,
                            .maximumSteps = 2};
                    openConfig.dirichletPlane       = p3mConfig.dirichletPlane;
                    auto p3mDirectConfig            = p3mConfig;
                    auto openDirectConfig           = openConfig;
                    p3mDirectConfig.dirichletPlane  = {};
                    openDirectConfig.dirichletPlane = {};
                    Run p3mImage(p3mConfig), openImage(openConfig);
                    Run p3mDirect(p3mDirectConfig), openDirect(openDirectConfig);
                    for (Run* run : {&p3mImage, &p3mDirect, &openImage, &openDirect}) {
                        setLongitudinalMomenta(run->particles, momenta);
                        placeClosePair(run->particles, planeZ);
                    }

                    for (std::size_t step : {0u, 1u}) {
                        EXPECT_EQ(p3mImage.solve(step).backendSolves, 2 * bins);
                        EXPECT_EQ(openImage.solve(step).backendSolves, 2 * bins);
                        EXPECT_EQ(p3mDirect.solve(step).backendSolves, bins);
                        EXPECT_EQ(openDirect.solve(step).backendSolves, bins);
                        auto host = [](const auto& attribute) {
                            return Kokkos::create_mirror_view_and_copy(
                                    Kokkos::HostSpace(), attribute.getView());
                        };
                        const auto p3mE         = host(p3mImage.particles.E),
                                   p3mB         = host(p3mImage.particles.B);
                        const auto p3mDirectE   = host(p3mDirect.particles.E),
                                   p3mDirectB   = host(p3mDirect.particles.B);
                        const auto openE        = host(openImage.particles.E),
                                   openB        = host(openImage.particles.B);
                        const auto openDirectE  = host(openDirect.particles.E),
                                   openDirectB  = host(openDirect.particles.B);
                        const std::size_t count = p3mImage.particles.getLocalNum();
                        const auto toP3MDirect =
                                matchByPosition(p3mImage.particles, p3mDirect.particles);
                        const auto toOpen =
                                matchByPosition(p3mImage.particles, openImage.particles);
                        const auto toOpenDirect =
                                matchByPosition(p3mImage.particles, openDirect.particles);
                        ASSERT_EQ(toP3MDirect.size(), count);
                        ASSERT_EQ(toOpen.size(), count);
                        ASSERT_EQ(toOpenDirect.size(), count);

                        // The image is mesh-resolved in both backends; P3M adds its particle term
                        // to the real bunch once, so it cancels in the image-minus-direct
                        // difference.
                        double electricNorm = 0.0, magneticScale = 0.0;
                        for (std::size_t i = 0; i < count; ++i) {
                            for (unsigned d = 0; d < 3; ++d) {
                                const double imageE =
                                        openE(toOpen[i])[d] - openDirectE(toOpenDirect[i])[d];
                                electricNorm += imageE * imageE;
                                magneticScale = std::max(
                                        magneticScale, std::abs(
                                                               openB(toOpen[i])[d]
                                                               - openDirectB(toOpenDirect[i])[d]));
                            }
                        }
                        EXPECT_GT(electricNorm, 1.0e-6);
                        EXPECT_GT(magneticScale, 0.0);
                        for (std::size_t i = 0; i < count; ++i) {
                            for (unsigned d = 0; d < 3; ++d) {
                                const double imageE =
                                        openE(toOpen[i])[d] - openDirectE(toOpenDirect[i])[d];
                                const double imageB =
                                        openB(toOpen[i])[d] - openDirectB(toOpenDirect[i])[d];
                                EXPECT_NEAR(
                                        p3mE(i)[d] - p3mDirectE(toP3MDirect[i])[d], imageE,
                                        1.0e-9 * std::max(1.0, std::abs(imageE)));
                                EXPECT_NEAR(
                                        p3mB(i)[d] - p3mDirectB(toP3MDirect[i])[d], imageB,
                                        1.0e-8 * magneticScale);
                            }
                        }
                    }
                    for (std::size_t step : {2u, 3u}) {
                        EXPECT_EQ(p3mImage.solve(step).backendSolves, bins);
                        EXPECT_EQ(p3mDirect.solve(step).backendSolves, bins);
                        expectFieldsEqualByPosition(p3mImage.particles, p3mDirect.particles);
                    }
                }
            }
        }

        TEST_F(CartesianPIC3DAlgorithmsTest, MeshTiedP3MCutoffMatchesEquivalentFixedCutoff) {
            auto meshConfig           = config();
            meshConfig.backend        = PoissonSolverType::P3M;
            meshConfig.p3mCutoffCells = 2.0;
            Run meshTied(meshConfig);

            // Bring two particles within the cutoff, then rescale the bunch so that the second
            // step changes the mesh spacing, the cutoff, and alpha on the live mesh solver. The P3M
            // solve sorts particles by cell, and a threaded sort returns particles that share a
            // cell in any order. The particles are identical, so write the whole position set
            // instead of moving one index relative to another.
            auto start = localHostCopy(meshTied.particles.R, meshTied.particles.getLocalNum());
            start(1)   = start(0) + Vector(0.0002, 0.0001, 0.0001);
            auto place = [&start](Particles& particles, double scale) {
                auto positions = particles.R.getHostMirror();
                Kokkos::deep_copy(positions, particles.R.getView());
                for (std::size_t i = 0; i < particles.getLocalNum(); ++i) {
                    positions(i) = start(i);
                    positions(i) *= scale;
                }
                Kokkos::deep_copy(particles.R.getView(), positions);
                particles.updateMoments();
            };
            // Beam-frame cutoff from the particle bounds, as CartesianDomainUpdater derives it.
            auto beamFrameCutoff = [&](Particles& particles) {
                auto positions = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), particles.R.getView());
                const double relativeIncrement = meshConfig.grid.boundingBoxIncreasePercent / 100.0;
                double largestSpacing          = 0.0;
                for (unsigned d = 0; d < 3; ++d) {
                    double lower = positions(0)[d];
                    double upper = positions(0)[d];
                    for (std::size_t i = 1; i < particles.getLocalNum(); ++i) {
                        lower = std::min(lower, positions(i)[d]);
                        upper = std::max(upper, positions(i)[d]);
                    }
                    const double span = (upper - lower) * (1.0 + 2.0 * relativeIncrement);
                    largestSpacing    = std::max(
                            largestSpacing,
                            span / static_cast<double>(meshConfig.grid.meshSize[d] - 1));
                }
                return meshConfig.p3mCutoffCells * largestSpacing;
            };

            for (const auto& [step, scale] : {std::pair{0u, 1.0}, std::pair{1u, 0.5}}) {
                place(meshTied.particles, scale);
                auto fixedConfig           = meshConfig;
                fixedConfig.p3mCutoffCells = 0.0;
                fixedConfig.p3mCutoff      = beamFrameCutoff(meshTied.particles);
                Run fixed(fixedConfig);
                place(fixed.particles, scale);
                meshTied.solve(step);
                fixed.solve(0);
                expectFieldsEqualByPosition(meshTied.particles, fixed.particles);
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
                auto values                 = config();
                values.backend              = backend;
                values.repartitionFrequency = 0;
                values.binning.emplace();
                values.binning->maximumBins            = 1;
                values.binning->adaptive               = false;
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
                    for (std::size_t i = 0; i < initial.size(); ++i)
                        positions(i) = initial[i];
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
                    EXPECT_EQ(
                            translated.solve(step, {solveToTracker.inverted(), solveToTracker})
                                    .backendSolves,
                            expectedSolves);
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
            auto values                 = config();
            values.backend              = PoissonSolverType::None;
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
            const double coordinateSpacing =
                    std::nextafter(pathTranslation, std::numeric_limits<double>::infinity())
                    - pathTranslation;
            const CoordinateSystemTrafo solveToTracker(
                    Vector(0.0, 0.0, pathTranslation), Quaternion(1.0, 0.0, 0.0, 0.0));
            EXPECT_EQ(centred.solve().backendSolves, 0u);
            EXPECT_EQ(
                    translated.solve(0, {solveToTracker.inverted(), solveToTracker}).backendSolves,
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
