#include <gtest/gtest.h>

#include "AbstractObjects/OpalData.h"
#include "Algorithms/Quaternion.hpp"
#include "PartBunch/BunchStateHandler.h"
#include "PartBunch/CartesianDomain.h"
#include "SpaceCharge/BarnesHut/BarnesHutAlgorithm.h"
#include "SpaceCharge/CartesianPIC3D/CartesianPIC3DAlgorithm.h"
#include "Structure/DataSink.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <vector>

// Tests for the Barnes-Hut space-charge algorithm. They run on any rank count: particle sets are
// defined by a global index, distributed round-robin over ranks, and compared by particle ID
// after gathering, because Barnes-Hut migrates and reorders particles.

namespace opalx::spacecharge {
    namespace {

        using Vector    = Vector_t<double, 3>;
        using Particles = ParticleContainer<double, 3>;

        constexpr double kTimeStep = 1.0e-12;
        constexpr double kCharge   = -1.0e-15;  // electron-like sign on purpose
        const double kCoulomb      = 1.0 / (4.0 * Physics::pi * Physics::epsilon_0);

        std::uint64_t splitmix(std::uint64_t z) {
            z += 0x9e3779b97f4a7c15ULL;
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            return z ^ (z >> 31);
        }

        double uniform(std::uint64_t index, std::uint64_t salt) {
            return static_cast<double>(splitmix(index * 0x100000001b3ULL ^ salt) >> 11)
                   * (1.0 / 9007199254740992.0);
        }

        double normal(std::uint64_t index, std::uint64_t salt) {
            const double u1 = std::max(uniform(index, salt), 1.0e-300);
            const double u2 = uniform(index, salt + 7919);
            return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * Physics::pi * u2);
        }

        enum class Shape { Gaussian, UniformSphere };

        struct BunchSpec {
            std::size_t globalCount = 2000;
            Shape shape             = Shape::Gaussian;
            Vector sigma{1.0e-3, 1.5e-3, 3.0e-3};  ///< rms sizes or sphere radius in x
            Vector center{0.0, 0.0, 0.02};
            bool variableDt = false;  ///< dt in (0.5, 1] * kTimeStep like emitted particles
        };

        Vector positionFor(const BunchSpec& spec, std::uint64_t k) {
            if (spec.shape == Shape::Gaussian) {
                return Vector(
                        spec.center[0] + spec.sigma[0] * normal(k, 1),
                        spec.center[1] + spec.sigma[1] * normal(k, 2),
                        spec.center[2] + spec.sigma[2] * normal(k, 3));
            }
            // Uniform ball of radius sigma[0]: radius ~ u^(1/3), isotropic direction.
            const double radius = spec.sigma[0] * std::cbrt(uniform(k, 11));
            const double cosT   = 2.0 * uniform(k, 12) - 1.0;
            const double sinT   = std::sqrt(std::max(0.0, 1.0 - cosT * cosT));
            const double phi    = 2.0 * Physics::pi * uniform(k, 13);
            return Vector(
                    spec.center[0] + radius * sinT * std::cos(phi),
                    spec.center[1] + radius * sinT * std::sin(phi),
                    spec.center[2] + radius * cosT);
        }

        /** @brief Host copy of one particle, gathered from every rank. */
        struct Record {
            std::int64_t id;
            double r[3], p[3], e[3], b[3];
            double dt;
        };

        std::vector<Record> gatherAll(Particles& particles) {
            const std::size_t n = particles.getLocalNum();
            auto R  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.R.getView());
            auto P  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.P.getView());
            auto E  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.E.getView());
            auto B  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.B.getView());
            auto dt = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.dt.getView());
            auto ID = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.ID.getView());
            std::vector<Record> local(n);
            for (std::size_t i = 0; i < n; ++i) {
                local[i].id = ID(i);
                local[i].dt = dt(i);
                for (unsigned d = 0; d < 3; ++d) {
                    local[i].r[d] = R(i)[d];
                    local[i].p[d] = P(i)[d];
                    local[i].e[d] = E(i)[d];
                    local[i].b[d] = B(i)[d];
                }
            }
            MPI_Comm comm  = ippl::Comm->getCommunicator();
            const int size = ippl::Comm->size();
            int bytes      = static_cast<int>(n * sizeof(Record));
            std::vector<int> counts(size), offsets(size);
            MPI_Allgather(&bytes, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
            std::exclusive_scan(counts.begin(), counts.end(), offsets.begin(), 0);
            std::vector<Record> all(
                    static_cast<std::size_t>(offsets.back() + counts.back()) / sizeof(Record));
            MPI_Allgatherv(
                    local.data(), bytes, MPI_BYTE, all.data(), counts.data(), offsets.data(),
                    MPI_BYTE, comm);
            std::sort(all.begin(), all.end(), [](const Record& a, const Record& b) {
                return a.id < b.id;
            });
            return all;
        }

        /** @brief Reference field of the softened Coulomb sum used by the ryoanji P2P kernel. */
        std::vector<std::array<double, 3>> directSum(
                const std::vector<Record>& particles, double charge, double softening) {
            const double h2 = 4.0 * softening * softening;
            std::vector<std::array<double, 3>> field(particles.size(), {0.0, 0.0, 0.0});
#pragma omp parallel for schedule(static)
            for (std::size_t i = 0; i < particles.size(); ++i) {
                double acc[3] = {0.0, 0.0, 0.0};
                for (std::size_t j = 0; j < particles.size(); ++j) {
                    const double dx[3] = {
                            particles[i].r[0] - particles[j].r[0],
                            particles[i].r[1] - particles[j].r[1],
                            particles[i].r[2] - particles[j].r[2]};
                    const double r2 = std::max(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2], h2);
                    const double q  = charge * particles[j].dt / kTimeStep;
                    const double w  = q / (r2 * std::sqrt(r2));
                    for (unsigned d = 0; d < 3; ++d) {
                        acc[d] += w * dx[d];
                    }
                }
                for (unsigned d = 0; d < 3; ++d) {
                    field[i][d] = kCoulomb * acc[d];
                }
            }
            return field;
        }

        /** @brief Small deterministic displacement of every local particle. */
        void drift(Particles& particles, double shift) {
            auto R = particles.R.getView();
            Kokkos::parallel_for(
                    "drift", particles.getLocalNum(), KOKKOS_LAMBDA(const std::size_t i) {
                        R(i)[0] += shift * (static_cast<double>(i % 3) - 1.0);
                        R(i)[2] += shift;
                    });
            Kokkos::fence();
        }

        double norm(const double* v) { return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }

        /** @brief sqrt(sum |E - E_ref|^2 / sum |E_ref|^2) over all particles. */
        double rmsRelativeError(
                const std::vector<Record>& actual, const std::vector<std::array<double, 3>>& reference) {
            double num = 0.0, den = 0.0;
            for (std::size_t i = 0; i < actual.size(); ++i) {
                for (unsigned d = 0; d < 3; ++d) {
                    const double diff = actual[i].e[d] - reference[i][d];
                    num += diff * diff;
                    den += reference[i][d] * reference[i][d];
                }
            }
            return std::sqrt(num / den);
        }

        class BarnesHutAlgorithmTest : public ::testing::Test {
        protected:
            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
                OpalData::getInstance()->storeInputFn("barnes_hut_algorithm.opal");
                gmsg                = new Inform(nullptr, -1);
                Options::enableHDF5 = false;
            }

            static void TearDownTestSuite() {
                delete gmsg;
                gmsg = nullptr;
                std::remove("barnes_hut_algorithm.stat");
                std::remove("barnes_hut_algorithm.lbal");
                ippl::finalize();
            }

            static BarnesHutConfig config(double softening = 1.0e-6, double theta = 0.5) {
                BarnesHutConfig result;
                result.theta     = theta;
                result.softening = softening;
                return result;
            }

            /** @brief Primary + secondary containers with an attached algorithm. */
            struct Run {
                CartesianDomain<double, 3> domain;
                std::shared_ptr<BunchStateHandler> state = std::make_shared<BunchStateHandler>();
                Particles particles;
                Particles secondary;
                std::unique_ptr<BarnesHutAlgorithm> algorithm;
                std::array<std::uint8_t, 2> activity{1, 0};

                Run(const BarnesHutConfig& values, const BunchSpec& spec, double charge = kCharge)
                    : domain(makeCartesianDomainConfig(values)),
                      particles(domain.mesh(), domain.layout()),
                      secondary(domain.mesh(), domain.layout()) {
                    particles.setBunchStateHandler(state);
                    secondary.setBunchStateHandler(state);
                    fill(spec);
                    particles.setQ(charge);
                    particles.setM(Physics::m_e);

                    secondary.createParticles(ippl::Comm->rank() == 0 ? 1 : 0);
                    secondary.setM(Physics::m_e);
                    Kokkos::deep_copy(secondary.R.getView(), Vector(0.0));
                    Kokkos::deep_copy(secondary.P.getView(), Vector(0.0));
                    Kokkos::deep_copy(secondary.E.getView(), Vector(7.0));
                    Kokkos::deep_copy(secondary.B.getView(), Vector(9.0));
                    const std::array containers{&particles, &secondary};
                    algorithm = std::make_unique<BarnesHutAlgorithm>(values, containers, state);
                }

                /// Round-robin: global particle k lives on rank k % size.
                void fill(const BunchSpec& spec) {
                    const std::size_t rank = ippl::Comm->rank();
                    const std::size_t size = ippl::Comm->size();
                    std::vector<std::uint64_t> mine;
                    for (std::size_t k = rank; k < spec.globalCount; k += size) {
                        mine.push_back(k);
                    }
                    const std::size_t offset = particles.getLocalNum();
                    particles.createParticles(mine.size());
                    auto R  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.R.getView());
                    auto P  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.P.getView());
                    auto dt = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.dt.getView());
                    for (std::size_t i = 0; i < mine.size(); ++i) {
                        const std::uint64_t k = mine[i];
                        R(offset + i)         = positionFor(spec, k);
                        P(offset + i) = Vector(0.01 * normal(k, 21), 0.01 * normal(k, 22), 1.0);
                        dt(offset + i) =
                                spec.variableDt ? kTimeStep * (0.5 + 0.5 * uniform(k, 31)) : kTimeStep;
                    }
                    Kokkos::deep_copy(particles.R.getView(), R);
                    Kokkos::deep_copy(particles.P.getView(), P);
                    Kokkos::deep_copy(particles.dt.getView(), dt);
                    particles.markMomentsDirty();
                }

                SpaceChargeSolveResult solve(CoordinateFrameTransforms frames = {}) {
                    SpaceChargeStepState step;
                    step.timeStep = kTimeStep;
                    step.mpiSize  = ippl::Comm->size();
                    step.frames   = frames;
                    return algorithm->solve(SpaceChargeSolveContext(activity, step));
                }
            };

            static void expectConserved(const std::vector<Record>& before, const std::vector<Record>& after) {
                ASSERT_EQ(before.size(), after.size());
                for (std::size_t i = 0; i < before.size(); ++i) {
                    ASSERT_EQ(before[i].id, after[i].id) << "particle set changed";
                    EXPECT_EQ(before[i].dt, after[i].dt);
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_EQ(before[i].r[d], after[i].r[d]);
                        EXPECT_EQ(before[i].p[d], after[i].p[d]);
                    }
                }
            }
        };

        TEST_F(BarnesHutAlgorithmTest, ConfigAndParserValidation) {
            EXPECT_NO_THROW(validateSpaceChargeConfig(SpaceChargeConfig(config())));
            for (double theta : {0.0, -0.1, 1.5, std::nan("")}) {
                auto values  = config();
                values.theta = theta;
                EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(values)), OpalException);
            }
            auto unsoftened      = config(0.0);
            EXPECT_THROW(validateSpaceChargeConfig(SpaceChargeConfig(unsoftened)), OpalException);
            unsoftened.leafBasedSoftening = true;
            EXPECT_NO_THROW(validateSpaceChargeConfig(SpaceChargeConfig(unsoftened)));
            const auto domain = makeCartesianDomainConfig(SpaceChargeConfig(config()));
            EXPECT_EQ(domain.layoutType, ParticleLayoutType::Spatial);
            EXPECT_FALSE(domain.periodicParticleBoundary);
        }

        TEST_F(BarnesHutAlgorithmTest, DirectSumTiny) {
            // Few particles end up in very few leaves, so nearly every interaction is P2P.
            BunchSpec spec;
            spec.globalCount = 48;
            Run run(config(2.0e-5), spec);
            const auto before = gatherAll(run.particles);
            EXPECT_EQ(run.solve().backendSolves, 1u);
            const auto after = gatherAll(run.particles);
            expectConserved(before, after);
            const auto reference = directSum(before, kCharge, 2.0e-5);
            const double error   = rmsRelativeError(after, reference);
            EXPECT_LT(error, ippl::Comm->size() == 1 ? 1.0e-10 : 5.0e-3);
            for (const Record& r : after) {
                EXPECT_EQ(r.b[0], 0.0);
                EXPECT_EQ(r.b[1], 0.0);
                EXPECT_EQ(r.b[2], 0.0);
            }
        }

        TEST_F(BarnesHutAlgorithmTest, DirectSumMediumConvergesWithTheta) {
            BunchSpec spec;
            spec.globalCount = 4000;
            spec.variableDt  = true;  // checks the Q dt_i / dt_step source weighting
            double previous  = 1.0;
            for (double theta : {0.7, 0.5, 0.3}) {
                SCOPED_TRACE(theta);
                Run run(config(1.0e-6, theta), spec);
                const auto before = gatherAll(run.particles);
                (void)run.solve();
                const auto after     = gatherAll(run.particles);
                const auto reference = directSum(before, kCharge, 1.0e-6);
                const double error   = rmsRelativeError(after, reference);
                if (ippl::Comm->rank() == 0) {
                    std::printf("[BH] theta=%.2f rms relative error=%.3e\n", theta, error);
                }
                EXPECT_LT(error, theta >= 0.7 ? 2.0e-2 : 5.0e-3);
                EXPECT_LT(error, previous);
                previous = error;
            }
        }

        TEST_F(BarnesHutAlgorithmTest, UniformSphereMatchesAnalyticField) {
            BunchSpec spec;
            spec.globalCount = 100000;
            spec.shape       = Shape::UniformSphere;
            spec.sigma       = Vector(1.0e-3, 0.0, 0.0);
            const double radius   = spec.sigma[0];
            const double spacing  = radius * std::cbrt(4.0 / 3.0 * Physics::pi / spec.globalCount);
            Run run(config(0.5 * spacing), spec);
            (void)run.solve();
            const auto after   = gatherAll(run.particles);
            const double total = kCharge * spec.globalCount;
            double num = 0.0, den = 0.0;
            for (const Record& rec : after) {
                const double rel[3] = {
                        rec.r[0] - spec.center[0], rec.r[1] - spec.center[1],
                        rec.r[2] - spec.center[2]};
                const double r = norm(rel);
                if (r < 0.2 * radius || r > 0.9 * radius) {
                    continue;  // discreteness dominates near the centre, the edge is not sharp
                }
                for (unsigned d = 0; d < 3; ++d) {
                    const double expected = kCoulomb * total * rel[d] / (radius * radius * radius);
                    num += (rec.e[d] - expected) * (rec.e[d] - expected);
                    den += expected * expected;
                }
            }
            const double error = std::sqrt(num / den);
            if (ippl::Comm->rank() == 0) {
                std::printf("[BH] uniform sphere rms relative error=%.3e\n", error);
            }
            EXPECT_LT(error, 5.0e-2);
        }

        TEST_F(BarnesHutAlgorithmTest, AgreesWithOpenFFTSolver) {
            BunchSpec spec;
            spec.globalCount = 200000;
            spec.sigma       = Vector(1.0e-3, 1.0e-3, 2.0e-3);

            CartesianPIC3DConfig pic;
            pic.backend            = PoissonSolverType::Open;
            pic.grid.meshSize      = {32, 32, 32};
            pic.grid.decomposition = {false, false, true};
            CartesianDomain<double, 3> domain(makeCartesianDomainConfig(pic));
            auto state = std::make_shared<BunchStateHandler>();
            Particles reference(domain.mesh(), domain.layout());
            reference.setBunchStateHandler(state);
            DataSink sink;
            {
                // Reuse the BH fixture's deterministic filler for an identical particle set.
                Run bh(config(2.0e-5), spec);
                const std::size_t rank = ippl::Comm->rank(), size = ippl::Comm->size();
                std::vector<std::uint64_t> mine;
                for (std::size_t k = rank; k < spec.globalCount; k += size) {
                    mine.push_back(k);
                }
                reference.createParticles(mine.size());
                reference.setQ(kCharge);
                reference.setM(Physics::m_e);
                auto R  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), reference.R.getView());
                for (std::size_t i = 0; i < mine.size(); ++i) {
                    R(i) = positionFor(spec, mine[i]);
                }
                Kokkos::deep_copy(reference.R.getView(), R);
                Kokkos::deep_copy(reference.P.getView(), Vector(0.0));
                Kokkos::deep_copy(reference.dt.getView(), kTimeStep);
                reference.updateMoments();
                const std::array containers{&reference};
                CartesianPIC3DAlgorithm open(
                        pic, containers,
                        std::make_unique<CartesianPIC3DFieldStorage<double, 3>>(domain), &sink,
                        state);
                const std::array<std::uint8_t, 1> activity{1};
                SpaceChargeStepState step;
                step.timeStep = kTimeStep;
                step.mpiSize  = ippl::Comm->size();
                (void)open.solve(SpaceChargeSolveContext(activity, step));
                (void)bh.solve();

                // Both containers assign IDs as rank + size * i in fill order: match by ID.
                const auto bhRecords  = gatherAll(bh.particles);
                const auto picRecords = gatherAll(reference);
                ASSERT_EQ(bhRecords.size(), picRecords.size());
                double num = 0.0, den = 0.0;
                for (std::size_t i = 0; i < bhRecords.size(); ++i) {
                    ASSERT_EQ(bhRecords[i].id, picRecords[i].id);
                    for (unsigned d = 0; d < 3; ++d) {
                        ASSERT_EQ(bhRecords[i].r[d], picRecords[i].r[d]);
                        const double diff = bhRecords[i].e[d] - picRecords[i].e[d];
                        num += diff * diff;
                        den += picRecords[i].e[d] * picRecords[i].e[d];
                    }
                }
                const double difference = std::sqrt(num / den);
                if (ippl::Comm->rank() == 0) {
                    std::printf("[BH] BH vs OPEN FFT (32^3) rms relative difference=%.3e\n",
                                difference);
                }
                EXPECT_LT(difference, 5.0e-2);
            }
        }

        TEST_F(BarnesHutAlgorithmTest, ReplacesFieldsAndPreservesOtherContainers) {
            BunchSpec spec;
            spec.globalCount = 500;
            Run run(config(), spec), clean(config(), spec);
            Kokkos::deep_copy(run.particles.E.getView(), Vector(123.0));
            Kokkos::deep_copy(run.particles.B.getView(), Vector(456.0));
            (void)run.solve();
            (void)clean.solve();
            const auto a = gatherAll(run.particles);
            const auto b = gatherAll(clean.particles);
            ASSERT_EQ(a.size(), b.size());
            for (std::size_t i = 0; i < a.size(); ++i) {
                for (unsigned d = 0; d < 3; ++d) {
                    EXPECT_EQ(a[i].e[d], b[i].e[d]);
                    EXPECT_EQ(a[i].b[d], 0.0);
                }
            }
            const auto secondary = gatherAll(run.secondary);
            ASSERT_EQ(secondary.size(), 1u);
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_EQ(secondary[0].e[d], 7.0);
                EXPECT_EQ(secondary[0].b[d], 9.0);
            }
        }

        TEST_F(BarnesHutAlgorithmTest, RotatedFrameRestoresStateAndRotatesFields) {
            BunchSpec spec;
            spec.globalCount = 1000;
            const CoordinateSystemTrafo pose(
                    Vector(0.1, 0.2, 0.3), Quaternion(std::cos(0.3), 0.0, std::sin(0.3), 0.0));
            Run baseline(config(), spec);
            Run rotated(config(), spec);
            pose.transformBunchTo(rotated.particles.R.getView(), rotated.particles.getLocalNum());
            pose.rotateBunchTo(rotated.particles.P.getView(), rotated.particles.getLocalNum());
            const auto original = gatherAll(rotated.particles);
            (void)baseline.solve();
            (void)rotated.solve({pose.inverted(), pose});
            const auto restored = gatherAll(rotated.particles);
            ASSERT_EQ(original.size(), restored.size());
            for (std::size_t i = 0; i < original.size(); ++i) {
                ASSERT_EQ(original[i].id, restored[i].id);
                for (unsigned d = 0; d < 3; ++d) {
                    EXPECT_NEAR(restored[i].r[d], original[i].r[d], 1.0e-14);
                    EXPECT_NEAR(restored[i].p[d], original[i].p[d], 1.0e-14);
                }
            }
            pose.rotateBunchTo(baseline.particles.E.getView(), baseline.particles.getLocalNum());
            const auto expected = gatherAll(baseline.particles);
            double scale        = 0.0;
            for (const Record& r : expected) {
                scale = std::max(scale, norm(r.e));
            }
            for (std::size_t i = 0; i < expected.size(); ++i) {
                ASSERT_EQ(expected[i].id, restored[i].id);
                for (unsigned d = 0; d < 3; ++d) {
                    // The rotated solve frame equals the baseline frame up to roundoff, which can
                    // flip SFC keys and therefore multipole groupings at the theta level.
                    EXPECT_NEAR(restored[i].e[d], expected[i].e[d], 1.0e-2 * scale);
                    EXPECT_EQ(restored[i].b[d], 0.0);
                }
            }
        }

        TEST_F(BarnesHutAlgorithmTest, TrivialPrimaryThenGrowFromOneRank) {
            for (std::size_t primaryCount : {0u, 1u}) {
                SCOPED_TRACE(primaryCount);
                BunchSpec none;
                none.globalCount = 0;
                Run run(config(), none);
                const bool root = ippl::Comm->rank() == 0;
                run.particles.createParticles(root ? primaryCount : 0);
                Kokkos::deep_copy(run.particles.R.getView(), Vector(0.1, -0.1, 0.2));
                Kokkos::deep_copy(run.particles.dt.getView(), kTimeStep);
                Kokkos::deep_copy(run.particles.E.getView(), Vector(123.0));
                Kokkos::deep_copy(run.particles.B.getView(), Vector(456.0));
                EXPECT_EQ(run.solve().backendSolves, 0u);
                const auto trivial = gatherAll(run.particles);
                ASSERT_EQ(trivial.size(), primaryCount);
                for (const Record& r : trivial) {
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_EQ(r.e[d], 0.0);
                        EXPECT_EQ(r.b[d], 0.0);
                    }
                }

                // All new particles start on rank 0 (like a cathode on one rank); the solve
                // spreads them over all ranks.
                BunchSpec grow;
                grow.globalCount = root ? 4096 : 0;
                if (root) {
                    grow.globalCount = 4096;
                }
                std::size_t created = 0;
                if (root) {
                    const std::size_t offset = run.particles.getLocalNum();
                    run.particles.createParticles(4096);
                    auto R  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), run.particles.R.getView());
                    auto dt = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), run.particles.dt.getView());
                    for (std::size_t i = 0; i < 4096; ++i) {
                        R(offset + i)  = positionFor(grow, i);
                        dt(offset + i) = kTimeStep;
                    }
                    Kokkos::deep_copy(run.particles.R.getView(), R);
                    Kokkos::deep_copy(run.particles.dt.getView(), dt);
                    created = 4096;
                } else {
                    run.particles.createParticles(0);
                }
                static_cast<void>(created);
                const auto before = gatherAll(run.particles);
                EXPECT_EQ(run.solve().backendSolves, 1u);
                EXPECT_EQ(run.particles.getTotalNum(), primaryCount + 4096);
                const auto after = gatherAll(run.particles);
                expectConserved(before, after);
                if (ippl::Comm->size() > 1) {
                    std::size_t local = run.particles.getLocalNum();
                    std::size_t minimum = 0;
                    ippl::Comm->allreduce(local, minimum, 1, std::less<std::size_t>());
                    EXPECT_GT(minimum, 0u) << "Barnes-Hut should spread particles over all ranks";
                }
                const auto reference = directSum(before, kCharge, 1.0e-6);
                EXPECT_LT(rmsRelativeError(after, reference), 5.0e-3);
            }
        }

        TEST_F(BarnesHutAlgorithmTest, CopyBackConservesParticlesAndDisablesLayoutUpdate) {
            BunchSpec spec;
            spec.globalCount = 20000;
            spec.variableDt  = true;
            Run run(config(), spec);
            const auto before = gatherAll(run.particles);
            (void)run.solve();
            const auto after = gatherAll(run.particles);
            expectConserved(before, after);
            std::size_t local = run.particles.getLocalNum(), total = 0;
            ippl::Comm->allreduce(local, total, 1, std::plus<std::size_t>());
            EXPECT_EQ(total, spec.globalCount);
            EXPECT_EQ(run.particles.getTotalNum(), spec.globalCount);
            std::set<std::int64_t> ids;
            for (const Record& r : after) {
                ids.insert(r.id);
            }
            EXPECT_EQ(ids.size(), spec.globalCount);
            EXPECT_TRUE(run.particles.isDecompositionOwnedExternally());
            EXPECT_THROW(run.particles.update(), OpalException);
        }

        TEST_F(BarnesHutAlgorithmTest, ShuffleAndSignFlipInvariance) {
            BunchSpec spec;
            spec.globalCount = 6000;
            Run run(config(), spec), negated(config(), spec, -kCharge), shuffled(config(), spec);
            // Reverse the local order of the shuffled run before its first solve.
            {
                auto R  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), shuffled.particles.R.getView());
                auto P  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), shuffled.particles.P.getView());
                auto dt = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), shuffled.particles.dt.getView());
                auto ID = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), shuffled.particles.ID.getView());
                const std::size_t n = shuffled.particles.getLocalNum();
                for (std::size_t i = 0; i < n / 2; ++i) {
                    std::swap(R(i), R(n - 1 - i));
                    std::swap(P(i), P(n - 1 - i));
                    std::swap(dt(i), dt(n - 1 - i));
                    std::swap(ID(i), ID(n - 1 - i));
                }
                Kokkos::deep_copy(shuffled.particles.R.getView(), R);
                Kokkos::deep_copy(shuffled.particles.P.getView(), P);
                Kokkos::deep_copy(shuffled.particles.dt.getView(), dt);
                Kokkos::deep_copy(shuffled.particles.ID.getView(), ID);
            }
            (void)run.solve();
            (void)negated.solve();
            (void)shuffled.solve();
            const auto a = gatherAll(run.particles);
            const auto n = gatherAll(negated.particles);
            const auto s = gatherAll(shuffled.particles);
            double scale = 0.0;
            for (const Record& r : a) {
                scale = std::max(scale, norm(r.e));
            }
            ASSERT_GT(scale, 0.0);
            for (std::size_t i = 0; i < a.size(); ++i) {
                ASSERT_EQ(a[i].id, n[i].id);
                ASSERT_EQ(a[i].id, s[i].id);
                for (unsigned d = 0; d < 3; ++d) {
                    EXPECT_NEAR(n[i].e[d], -a[i].e[d], 1.0e-12 * scale);
                    EXPECT_NEAR(s[i].e[d], a[i].e[d], 1.0e-12 * scale);
                }
            }
        }

        TEST_F(BarnesHutAlgorithmTest, IncrementalMigrationOverManySolves) {
            BunchSpec spec;
            spec.globalCount   = 50000;
            Run run(config(), spec);
            const std::size_t steps = 100;
            std::size_t maxMoved    = 0;
            for (std::size_t step = 0; step < steps; ++step) {
                // Record which IDs this rank owns, drift every particle a little, then solve.
                auto IDs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), run.particles.ID.getView());
                std::set<std::int64_t> owned(IDs.data(), IDs.data() + run.particles.getLocalNum());
                drift(run.particles, 2.0e-6);
                ASSERT_EQ(run.solve().backendSolves, 1u);
                auto after = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), run.particles.ID.getView());
                std::size_t arrived = 0;
                for (std::size_t i = 0; i < run.particles.getLocalNum(); ++i) {
                    arrived += owned.count(after(i)) == 0 ? 1 : 0;
                }
                std::size_t moved = 0;
                ippl::Comm->allreduce(arrived, moved, 1, std::plus<std::size_t>());
                if (step > 0) {
                    maxMoved = std::max(maxMoved, moved);
                }
                ASSERT_EQ(run.particles.getTotalNum(), spec.globalCount);
            }
            const auto records = gatherAll(run.particles);
            for (const Record& r : records) {
                ASSERT_TRUE(std::isfinite(r.e[0]) && std::isfinite(r.e[1]) && std::isfinite(r.e[2]));
            }
            const auto& stats = run.algorithm->statistics();
            EXPECT_EQ(stats.solves, steps);
            EXPECT_EQ(stats.containerBuilds, 1u);
            if (ippl::Comm->rank() == 0) {
                std::printf("[BH] max particles changing rank per step: %zu of %zu, in-place packs "
                            "%zu of %zu\n",
                            maxMoved, spec.globalCount, stats.inPlacePacks, stats.solves);
            }
            EXPECT_LT(maxMoved, spec.globalCount / 20);
        }

        TEST_F(BarnesHutAlgorithmTest, RebuildsContainerWhenBunchGrowsAndMatchesFreshSolve) {
            BunchSpec small;
            small.globalCount = 500;
            Run run(config(), small);
            (void)run.solve();
            EXPECT_EQ(run.algorithm->statistics().containerBuilds, 1u);

            BunchSpec large   = small;
            large.globalCount = 200000;
            // Add global indices [500, 200000) to the existing set.
            {
                const std::size_t rank = ippl::Comm->rank(), size = ippl::Comm->size();
                std::vector<std::uint64_t> mine;
                for (std::size_t k = small.globalCount + rank; k < large.globalCount; k += size) {
                    mine.push_back(k);
                }
                const std::size_t offset = run.particles.getLocalNum();
                run.particles.createParticles(mine.size());
                auto R  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), run.particles.R.getView());
                auto dt = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), run.particles.dt.getView());
                for (std::size_t i = 0; i < mine.size(); ++i) {
                    R(offset + i)  = positionFor(large, mine[i]);
                    dt(offset + i) = kTimeStep;
                }
                Kokkos::deep_copy(run.particles.R.getView(), R);
                Kokkos::deep_copy(run.particles.dt.getView(), dt);
            }
            (void)run.solve();
            EXPECT_EQ(run.algorithm->statistics().containerBuilds, 2u);
            EXPECT_EQ(run.particles.getTotalNum(), large.globalCount);

            Run fresh(config(), large);
            (void)fresh.solve();
            // Same physical particles, but IDs differ between the two runs: match by position.
            auto key = [](const Record& r) { return std::array<double, 3>{r.r[0], r.r[1], r.r[2]}; };
            std::map<std::array<double, 3>, std::array<double, 3>> freshField;
            for (const Record& r : gatherAll(fresh.particles)) {
                freshField[key(r)] = {r.e[0], r.e[1], r.e[2]};
            }
            double num = 0.0, den = 0.0;
            for (const Record& r : gatherAll(run.particles)) {
                const auto found = freshField.find(key(r));
                ASSERT_NE(found, freshField.end());
                for (unsigned d = 0; d < 3; ++d) {
                    num += (r.e[d] - found->second[d]) * (r.e[d] - found->second[d]);
                    den += found->second[d] * found->second[d];
                }
            }
            EXPECT_LT(std::sqrt(num / den), 1.0e-3);
        }

    }  // namespace
}  // namespace opalx::spacecharge
