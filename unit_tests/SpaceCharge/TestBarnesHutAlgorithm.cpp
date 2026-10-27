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
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <string>
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
            double transverseMomentum   = 0.01;  ///< rms px, py (beta gamma)
            double longitudinalMomentum = 1.0;   ///< pz (beta gamma)
            double momentumSpread       = 0.0;   ///< relative rms spread of pz
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

        Vector momentumFor(const BunchSpec& spec, std::uint64_t k) {
            return Vector(
                    spec.transverseMomentum * normal(k, 21), spec.transverseMomentum * normal(k, 22),
                    spec.longitudinalMomentum * (1.0 + spec.momentumSpread * normal(k, 23)));
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

        /**
         * @brief Reference field of the softened Coulomb sum used by the ryoanji P2P kernel, at
         * @p targets from @p sources (source charge = charge * dt / kTimeStep).
         */
        std::vector<std::array<double, 3>> softenedField(
                const std::vector<Record>& targets, const std::vector<Record>& sources,
                double charge, double softening) {
            const double h2 = 4.0 * softening * softening;
            std::vector<std::array<double, 3>> field(targets.size(), {0.0, 0.0, 0.0});
#pragma omp parallel for schedule(static)
            for (std::size_t i = 0; i < targets.size(); ++i) {
                double acc[3] = {0.0, 0.0, 0.0};
                for (const Record& src : sources) {
                    const double dx[3] = {
                            targets[i].r[0] - src.r[0], targets[i].r[1] - src.r[1],
                            targets[i].r[2] - src.r[2]};
                    const double r2 = std::max(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2], h2);
                    const double q  = charge * src.dt / kTimeStep;
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

        std::vector<std::array<double, 3>> directSum(
                const std::vector<Record>& particles, double charge, double softening) {
            return softenedField(particles, particles, charge, softening);
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

        /** @brief Displacement that depends only on the particle ID, not on its local slot. */
        void driftById(Particles& particles, double shift) {
            auto R        = particles.R.getView();
            const auto ID = particles.ID.getView();
            Kokkos::parallel_for(
                    "driftById", particles.getLocalNum(), KOKKOS_LAMBDA(const std::size_t i) {
                        R(i)[0] += shift * (static_cast<double>(ID(i) % 3) - 1.0);
                        R(i)[2] += shift;
                    });
            Kokkos::fence();
        }

        /** @brief Apply @p f to every local position on the host. */
        template <class F>
        void transformPositions(Particles& particles, F f) {
            auto R = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.R.getView());
            for (std::size_t i = 0; i < particles.getLocalNum(); ++i) {
                f(R(i));
            }
            Kokkos::deep_copy(particles.R.getView(), R);
            particles.markMomentsDirty();
        }

        /** @brief Reference field with a grounded plane: every particle has a mirror image
         * (z' = 2 planeZ - z, opposite charge) that acts as a source only. */
        std::vector<std::array<double, 3>> directSumWithImages(
                const std::vector<Record>& particles, double charge, double softening, double planeZ) {
            std::vector<Record> sources = particles;
            for (const Record& r : particles) {
                Record image = r;
                image.r[2]   = 2.0 * planeZ - r.r[2];
                image.dt     = -r.dt;  // the source charge is charge * dt / kTimeStep
                sources.push_back(image);
            }
            return softenedField(particles, sources, charge, softening);
        }

        using Fields = std::vector<std::array<double, 3>>;

        /**
         * @brief Reference of binned CartesianPIC3D with one bin: the softened field of the bunch
         * with z stretched by gamma = <sqrt(1 + p^2)> (and of its images mirrored at the stretched
         * plane), boosted with v = c <p> / gamma; image charges move the other way and enter B
         * with the opposite sign.
         */
        struct RestFrameReference {
            Fields e, b;
            double gamma = 1.0;
        };

        /** @brief Mean gamma = <sqrt(1 + p^2)> and v = c <p> / gamma of @p particles. */
        double meanRestFrame(const std::vector<Record>& particles, double v[3]) {
            const double n = static_cast<double>(particles.size());
            double meanP[3] = {0.0, 0.0, 0.0}, gamma = 0.0;
            for (const Record& r : particles) {
                for (unsigned d = 0; d < 3; ++d) {
                    meanP[d] += r.p[d] / n;
                }
                gamma += std::sqrt(1.0 + r.p[0] * r.p[0] + r.p[1] * r.p[1] + r.p[2] * r.p[2]) / n;
            }
            for (unsigned d = 0; d < 3; ++d) {
                v[d] = Physics::c * meanP[d] / gamma;
            }
            return gamma;
        }

        /** @brief Boost the rest-frame field of all sources @p all (and of the image charges
         * alone, @p image, which enter B with the opposite sign) to the lab frame. */
        RestFrameReference boostReference(
                const std::vector<Record>& particles, const Fields& all, const Fields& image) {
            double v[3];
            const double gamma = meanRestFrame(particles, v);
            const double speed = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            RestFrameReference result;
            result.gamma = gamma;
            result.e.resize(particles.size());
            result.b.resize(particles.size());
            for (std::size_t i = 0; i < particles.size(); ++i) {
                double along = 0.0;
                for (unsigned d = 0; d < 3; ++d) {
                    along += speed > 0.0 ? all[i][d] * v[d] / speed : 0.0;
                }
                double moving[3];
                for (unsigned d = 0; d < 3; ++d) {
                    const double unit = speed > 0.0 ? v[d] / speed : 0.0;
                    result.e[i][d]    = gamma * all[i][d] - (gamma - 1.0) * along * unit;
                    moving[d]         = all[i][d] - 2.0 * image[i][d];
                }
                const double scale = gamma / (Physics::c * Physics::c);
                result.b[i]        = {scale * (v[1] * moving[2] - v[2] * moving[1]),
                                      scale * (v[2] * moving[0] - v[0] * moving[2]),
                                      scale * (v[0] * moving[1] - v[1] * moving[0])};
            }
            return result;
        }

        RestFrameReference restFrameDirectSum(
                const std::vector<Record>& particles, double charge, double softening,
                std::optional<double> planeZ = std::nullopt) {
            double v[3];
            const double gamma = meanRestFrame(particles, v);
            std::vector<Record> stretched = particles, images;
            for (Record& r : stretched) {
                r.r[2] *= gamma;
            }
            if (planeZ.has_value()) {
                for (const Record& r : stretched) {
                    Record image = r;
                    image.r[2]   = 2.0 * gamma * *planeZ - r.r[2];
                    image.dt     = -r.dt;
                    images.push_back(image);
                }
            }
            std::vector<Record> sources = stretched;
            sources.insert(sources.end(), images.begin(), images.end());
            const Fields all = softenedField(stretched, sources, charge, softening);
            const Fields image = images.empty() ? Fields(particles.size(), {0.0, 0.0, 0.0})
                                                : softenedField(stretched, images, charge, softening);
            return boostReference(particles, all, image);
        }

        /**
         * @brief Continuum field of the Gaussian bunch of @p spec (total charge @p total) at the
         * particles, in its rest frame (z stretched by gamma) and boosted to the lab frame. With
         * @p planeZ the grounded plane adds the image Gaussian (opposite charge, centre mirrored
         * at the stretched plane), which enters B with the opposite sign.
         *
         * E_i(x) = k Q x_i / sqrt(2 pi) int_0^inf du exp(-sum_j x_j^2 / (2 (s_j^2 + u)))
         *          / ((s_i^2 + u) prod_j sqrt(s_j^2 + u)), with u = e^t and the trapezoidal rule
         * in t (the integrand decays exponentially at both ends).
         */
        RestFrameReference gaussianContinuum(
                const std::vector<Record>& particles, const BunchSpec& spec, double total,
                std::optional<double> planeZ = std::nullopt) {
            double v[3];
            const double gamma = meanRestFrame(particles, v);
            const double s2[3] = {
                    spec.sigma[0] * spec.sigma[0], spec.sigma[1] * spec.sigma[1],
                    gamma * gamma * spec.sigma[2] * spec.sigma[2]};
            const double smallest = std::min({s2[0], s2[1], s2[2]});
            const double largest  = std::max({s2[0], s2[1], s2[2]});
            const double first = std::log(smallest) - 25.0, last = std::log(largest) + 40.0;
            const int steps    = static_cast<int>((last - first) / 0.05);
            const double dt    = (last - first) / steps;
            const double scale = kCoulomb * total / std::sqrt(2.0 * Physics::pi);
            // Field of the Gaussian at offset x from its centre (rest frame).
            const auto gaussian = [&](const double x[3], double field[3]) {
                double integral[3] = {0.0, 0.0, 0.0};
                for (int k = 0; k <= steps; ++k) {
                    const double u      = std::exp(first + k * dt);
                    const double weight = (k == 0 || k == steps) ? 0.5 : 1.0;
                    double exponent = 0.0, product = 1.0;
                    for (unsigned d = 0; d < 3; ++d) {
                        exponent += x[d] * x[d] / (2.0 * (s2[d] + u));
                        product *= std::sqrt(s2[d] + u);
                    }
                    const double f = weight * u * std::exp(-exponent) / product;
                    for (unsigned d = 0; d < 3; ++d) {
                        integral[d] += f / (s2[d] + u);
                    }
                }
                for (unsigned d = 0; d < 3; ++d) {
                    field[d] = scale * x[d] * integral[d] * dt;
                }
            };
            Fields all(particles.size()), image(particles.size(), {0.0, 0.0, 0.0});
#pragma omp parallel for schedule(static)
            for (std::size_t i = 0; i < particles.size(); ++i) {
                const double z = gamma * particles[i].r[2], centre = gamma * spec.center[2];
                const double x[3] = {
                        particles[i].r[0] - spec.center[0], particles[i].r[1] - spec.center[1],
                        z - centre};
                double real[3];
                gaussian(x, real);
                for (unsigned d = 0; d < 3; ++d) {
                    all[i][d] = real[d];
                }
                if (planeZ.has_value()) {
                    const double mirrored[3] = {x[0], x[1], z - (2.0 * gamma * *planeZ - centre)};
                    double field[3];
                    gaussian(mirrored, field);
                    for (unsigned d = 0; d < 3; ++d) {
                        image[i][d] = -field[d];
                        all[i][d] -= field[d];
                    }
                }
            }
            return boostReference(particles, all, image);
        }

        Fields fieldsOf(const std::vector<Record>& records, double (Record::*member)[3]) {
            Fields result(records.size());
            for (std::size_t i = 0; i < records.size(); ++i) {
                for (unsigned d = 0; d < 3; ++d) {
                    result[i][d] = (records[i].*member)[d];
                }
            }
            return result;
        }

        /** @brief sqrt(sum |a - b|^2 / sum |b|^2). */
        double rmsRelativeDifference(const Fields& a, const Fields& b) {
            double num = 0.0, den = 0.0;
            for (std::size_t i = 0; i < a.size(); ++i) {
                for (unsigned d = 0; d < 3; ++d) {
                    num += (a[i][d] - b[i][d]) * (a[i][d] - b[i][d]);
                    den += b[i][d] * b[i][d];
                }
            }
            return std::sqrt(num / den);
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
                        P(offset + i)         = momentumFor(spec, k);
                        dt(offset + i) =
                                spec.variableDt ? kTimeStep * (0.5 + 0.5 * uniform(k, 31)) : kTimeStep;
                    }
                    Kokkos::deep_copy(particles.R.getView(), R);
                    Kokkos::deep_copy(particles.P.getView(), P);
                    Kokkos::deep_copy(particles.dt.getView(), dt);
                    particles.markMomentsDirty();
                }

                SpaceChargeSolveResult solve(
                        CoordinateFrameTransforms frames = {}, std::size_t stepIndex = 0) {
                    SpaceChargeStepState step;
                    step.step     = stepIndex;
                    step.timeStep = kTimeStep;
                    step.mpiSize  = ippl::Comm->size();
                    step.frames   = frames;
                    return algorithm->solve(SpaceChargeSolveContext(activity, step));
                }
            };

            static BarnesHutConfig restFrameConfig(double softening = 2.0e-5, double theta = 0.5) {
                BarnesHutConfig result = config(softening, theta);
                result.kinematics      = BarnesHutKinematics::Frame;
                return result;
            }

            /// Like expectConserved, but z may differ by the rounding of the rest-frame stretch.
            static void expectConservedUpToStretch(
                    const std::vector<Record>& before, const std::vector<Record>& after) {
                ASSERT_EQ(before.size(), after.size());
                for (std::size_t i = 0; i < before.size(); ++i) {
                    ASSERT_EQ(before[i].id, after[i].id) << "particle set changed";
                    EXPECT_EQ(before[i].dt, after[i].dt);
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_NEAR(before[i].r[d], after[i].r[d], 1.0e-15);
                        EXPECT_EQ(before[i].p[d], after[i].p[d]);
                    }
                }
            }

            /**
             * @brief E and B of binned CartesianPIC3D (OPEN, one bin) on a @p mesh^3 grid for the
             * particle set of @p spec, sorted by ID like gatherAll().
             */
            static std::vector<Record> binnedOpenFields(
                    const BunchSpec& spec, std::size_t mesh, DirichletPlaneConfig plane) {
                CartesianPIC3DConfig pic;
                pic.backend            = PoissonSolverType::Open;
                pic.grid.meshSize      = {mesh, mesh, mesh};
                pic.grid.decomposition = {false, false, true};
                pic.dirichletPlane     = plane;
                pic.binning.emplace();
                pic.binning->maximumBins = 1;
                pic.binning->adaptive    = false;
                CartesianDomain<double, 3> domain(makeCartesianDomainConfig(pic));
                auto state = std::make_shared<BunchStateHandler>();
                Particles reference(domain.mesh(), domain.layout());
                reference.setBunchStateHandler(state);
                DataSink sink;
                // Same IDs (rank + size * i) and particles as Run::fill.
                const std::size_t rank = ippl::Comm->rank(), size = ippl::Comm->size();
                std::vector<std::uint64_t> mine;
                for (std::size_t k = rank; k < spec.globalCount; k += size) {
                    mine.push_back(k);
                }
                reference.createParticles(mine.size());
                reference.setQ(kCharge);
                reference.setM(Physics::m_e);
                auto R = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), reference.R.getView());
                auto P = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), reference.P.getView());
                for (std::size_t i = 0; i < mine.size(); ++i) {
                    R(i) = positionFor(spec, mine[i]);
                    P(i) = momentumFor(spec, mine[i]);
                }
                Kokkos::deep_copy(reference.R.getView(), R);
                Kokkos::deep_copy(reference.P.getView(), P);
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
                return gatherAll(reference);
            }

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
            // Fewer particles than the direct-sum limit (64 per rank).
            BunchSpec spec;
            spec.globalCount = 48;
            Run run(config(2.0e-5), spec);
            const auto before = gatherAll(run.particles);
            EXPECT_EQ(run.solve().backendSolves, 1u);
            const auto after = gatherAll(run.particles);
            expectConserved(before, after);
            const auto reference = directSum(before, kCharge, 2.0e-5);
            const double error   = rmsRelativeError(after, reference);
            // The direct sum evaluates the reference kernel; only summation order differs.
            EXPECT_LT(error, 1.0e-10);
            EXPECT_EQ(run.algorithm->statistics().directSolves, 1u);
            for (const Record& r : after) {
                EXPECT_EQ(r.b[0], 0.0);
                EXPECT_EQ(r.b[1], 0.0);
                EXPECT_EQ(r.b[2], 0.0);
            }
        }

        TEST_F(BarnesHutAlgorithmTest, ForcedDirectSumMatchesReference) {
            // BHDIRECT: a bunch far above the fallback limit is still summed directly, with and
            // without image charges, and never enters the tree.
            BunchSpec spec;
            spec.globalCount = 4000;
            spec.variableDt  = true;
            for (bool images : {false, true}) {
                SCOPED_TRACE(images);
                auto values      = config(2.0e-5);
                values.directSum = true;
                const double planeZ = spec.center[2] - 2.0 * spec.sigma[2];
                if (images) {
                    values.dirichletPlane = {.kind   = DirichletPlaneType::ImageCharge,
                                             .planeZ = planeZ};
                }
                Run run(values, spec);
                const auto before = gatherAll(run.particles);
                EXPECT_EQ(run.solve().backendSolves, 1u);
                const auto after = gatherAll(run.particles);
                expectConserved(before, after);
                const auto reference = images ? directSumWithImages(before, kCharge, 2.0e-5, planeZ)
                                              : directSum(before, kCharge, 2.0e-5);
                EXPECT_LT(rmsRelativeError(after, reference), 1.0e-10);
                EXPECT_EQ(run.algorithm->statistics().directSolves, 1u);
                EXPECT_EQ(run.algorithm->statistics().solves, 0u);
            }
        }

        TEST_F(BarnesHutAlgorithmTest, SmallAndShrinkingBunchesMatchDirectSum) {
            // Below 64 particles per rank the algorithm falls back to ryoanji's direct sum, which
            // uses the reference kernel exactly; above it Barnes-Hut must cope with small trees.
            const std::size_t directLimit = 64 * static_cast<std::size_t>(ippl::Comm->size());
            auto check = [&](Run& run, std::size_t count, const std::vector<Record>& before) {
                const auto after = gatherAll(run.particles);
                expectConserved(before, after);
                const double error = rmsRelativeError(after, directSum(before, kCharge, 2.0e-5));
                if (ippl::Comm->rank() == 0) {
                    std::printf("[BH] N=%zu %s rms relative error=%.3e\n", count,
                                count < directLimit ? "direct" : "tree", error);
                }
                // Tree limit ~3x the largest measured error (1.49e-4, N=1000 on 4 ranks): an
                // accuracy regression must fail here.
                EXPECT_LT(error, count < directLimit ? 1.0e-10 : 5.0e-4);
            };

            for (std::size_t count : {2u, 9u, 48u, 65u, 300u, 1000u}) {
                SCOPED_TRACE(count);
                BunchSpec spec;
                spec.globalCount = count;
                Run run(config(2.0e-5), spec);
                const auto before = gatherAll(run.particles);
                EXPECT_EQ(run.solve().backendSolves, 1u);
                check(run, count, before);
                const auto& stats = run.algorithm->statistics();
                EXPECT_EQ(stats.directSolves, count < directLimit ? 1u : 0u);
                EXPECT_EQ(stats.solves, count < directLimit ? 0u : 1u);
            }

            // Shrink from a tree-sized bunch to 20 particles (every 20th ID survives).
            BunchSpec spec;
            spec.globalCount = std::max<std::size_t>(400, 2 * directLimit);
            Run run(config(2.0e-5), spec);
            (void)run.solve();
            EXPECT_EQ(run.algorithm->statistics().solves, 1u);
            const std::int64_t stride = static_cast<std::int64_t>(spec.globalCount / 20);
            {
                auto ID      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), run.particles.ID.getView());
                auto invalid = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), run.particles.InvalidMask.getView());
                for (std::size_t i = 0; i < run.particles.getLocalNum(); ++i) {
                    // IDs are 0..N-1 (rank + size * i over the round-robin fill); keep 20 of them.
                    invalid(i) = ID(i) % stride != 0 || ID(i) >= 20 * stride;
                }
                Kokkos::deep_copy(run.particles.InvalidMask.getView(), invalid);
                (void)run.particles.deleteInvalidParticles();
            }
            ASSERT_EQ(run.particles.getTotalNum(), 20u);
            const auto before = gatherAll(run.particles);
            EXPECT_EQ(run.solve().backendSolves, 1u);
            check(run, 20, before);
            EXPECT_EQ(run.algorithm->statistics().directSolves, 1u);
        }

        TEST_F(BarnesHutAlgorithmTest, ImageChargesMatchDirectSumWithMirrors) {
            // Grounded plane 2 sigma_z behind the bunch centre (image field ~8% of the self-field,
            // a few particles lie behind it, so mirrors mix into the bunch tail) and 5 sigma_z
            // behind it (no particle behind the plane; reals and mirrors share only cells that
            // straddle the plane). The reference treats both identically. The direct-sum path
            // (48 particles) uses the reference kernel exactly; the tree path must stay at the
            // usual BH accuracy. Only real particles may come back to OPALX.
            for (double planeOffset : {2.0, 5.0}) {
                for (std::size_t count : {48u, 4000u}) {
                    SCOPED_TRACE(count);
                    SCOPED_TRACE(planeOffset);
                    BunchSpec spec;
                    spec.globalCount = count;
                    spec.variableDt  = true;
                    auto values      = config(2.0e-5);
                    values.dirichletPlane = {
                            .kind   = DirichletPlaneType::ImageCharge,
                            .planeZ = spec.center[2] - planeOffset * spec.sigma[2]};
                    Run run(values, spec);
                    const auto before = gatherAll(run.particles);
                    EXPECT_EQ(run.solve().backendSolves, 1u);
                    const auto after = gatherAll(run.particles);
                    expectConserved(before, after);
                    const double error = rmsRelativeError(
                            after, directSumWithImages(before, kCharge, 2.0e-5,
                                                       values.dirichletPlane.planeZ));
                    // The image field must actually matter: without mirrors the field is different.
                    const double withoutImages =
                            rmsRelativeError(after, directSum(before, kCharge, 2.0e-5));
                    const bool direct = count < 64 * static_cast<std::size_t>(ippl::Comm->size());
                    if (ippl::Comm->rank() == 0) {
                        std::printf("[BH] images plane=%.0f sigma_z N=%zu %s rms relative error=%.3e"
                                    " (vs no images %.3e)\n",
                                    planeOffset, count, direct ? "direct" : "tree", error,
                                    withoutImages);
                    }
                    // Tree measured 2.21e-4 (2 sigma_z) and 2.70e-4 (5 sigma_z), near the
                    // image-free level (1.73e-4); without dipole multipoles 3.1e-3 at 2 sigma_z
                    // (cubic domain).
                    EXPECT_LT(error, direct ? 1.0e-10 : 7.0e-4);
                    EXPECT_GT(withoutImages, planeOffset < 3.0 ? 1.0e-2 : 2.0e-3);
                    EXPECT_EQ(run.algorithm->statistics().dipoleSolves, direct ? 0u : 1u);
                }
            }
        }

        TEST_F(BarnesHutAlgorithmTest, ImageChargeErrorScalesWithTheta) {
            // Cells mixing reals and mirrors converge like the image-free case once their dipole
            // is kept (measured 9.11e-4, 2.21e-4, 2.06e-5; DirectSumMediumConvergesWithTheta:
            // 7.28e-4, 1.73e-4, 1.64e-5), so the image-free limits apply. Without the dipole the
            // error fell only like theta^2 (5.99e-3, 3.10e-3, 1.01e-3, cubic domain).
            BunchSpec spec;
            spec.globalCount = 4000;
            spec.variableDt  = true;
            const double planeZ = spec.center[2] - 2.0 * spec.sigma[2];
            double previous     = 1.0;
            for (double theta : {0.7, 0.5, 0.3}) {
                SCOPED_TRACE(theta);
                auto values           = config(2.0e-5, theta);
                values.dirichletPlane = {.kind = DirichletPlaneType::ImageCharge, .planeZ = planeZ};
                Run run(values, spec);
                const auto before = gatherAll(run.particles);
                (void)run.solve();
                const auto after   = gatherAll(run.particles);
                const double error = rmsRelativeError(
                        after, directSumWithImages(before, kCharge, 2.0e-5, planeZ));
                if (ippl::Comm->rank() == 0) {
                    std::printf("[BH] images theta=%.2f rms relative error=%.3e\n", theta, error);
                }
                EXPECT_LT(error, theta >= 0.7 ? 2.0e-3 : (theta >= 0.5 ? 5.0e-4 : 5.0e-5));
                EXPECT_LT(error, previous);
                previous = error;
            }
        }

        TEST_F(BarnesHutAlgorithmTest, GroundedPlaneHasNoTangentialField) {
            // Particles placed exactly on the plane see no tangential field: every charge and its
            // mirror contribute equal and opposite E_x, E_y there. The normal field is non-zero.
            BunchSpec spec;
            spec.globalCount = 4000;
            const double planeZ = spec.center[2] - 4.0 * spec.sigma[2];
            auto values         = config(2.0e-5);
            values.dirichletPlane = {.kind = DirichletPlaneType::ShiftedGreen, .planeZ = planeZ};
            Run run(values, spec);
            // Move every 10th global particle onto the plane (by ID, any rank layout).
            {
                auto R  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), run.particles.R.getView());
                auto ID = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), run.particles.ID.getView());
                for (std::size_t i = 0; i < run.particles.getLocalNum(); ++i) {
                    if (ID(i) % 10 == 0) {
                        R(i)[2] = planeZ;
                    }
                }
                Kokkos::deep_copy(run.particles.R.getView(), R);
            }
            (void)run.solve();
            double tangential = 0.0, normal = 0.0;
            std::size_t onPlane = 0;
            for (const Record& r : gatherAll(run.particles)) {
                if (r.r[2] != planeZ) {
                    continue;
                }
                ++onPlane;
                tangential += r.e[0] * r.e[0] + r.e[1] * r.e[1];
                normal += r.e[2] * r.e[2];
            }
            ASSERT_EQ(onPlane, spec.globalCount / 10);
            ASSERT_GT(normal, 0.0);
            if (ippl::Comm->rank() == 0) {
                std::printf("[BH] grounded plane: rms E_t / rms E_n = %.3e\n",
                            std::sqrt(tangential / normal));
            }
            // The tangential field cancels exactly only where a source cell and its mirror get the
            // same multipole decision; otherwise the residual is at the tree-error level. Measured
            // 9.83e-6 (4 ranks) with compact target groups (cstone group split fix, 2026-10-16);
            // 3.66e-8 before, when unsplit groups in mixed-dimension boxes made almost every
            // near interaction P2P (1.85e-6 with the cubic domain).
            EXPECT_LT(std::sqrt(tangential / normal), 2.5e-5);
        }

        TEST_F(BarnesHutAlgorithmTest, ImageChargesExpireAfterMaximumSteps) {
            BunchSpec spec;
            spec.globalCount = 2000;
            auto values      = config(2.0e-5);
            values.dirichletPlane = {.kind         = DirichletPlaneType::ImageCharge,
                                     .planeZ       = spec.center[2] - 4.0 * spec.sigma[2],
                                     .maximumSteps = 2};
            Run run(values, spec);
            for (std::size_t step = 0; step < 4; ++step) {
                SCOPED_TRACE(step);
                const auto before = gatherAll(run.particles);
                ASSERT_EQ(run.solve({}, step).backendSolves, 1u);
                const auto after = gatherAll(run.particles);
                expectConserved(before, after);
                const auto reference =
                        step < 2 ? directSumWithImages(before, kCharge, 2.0e-5,
                                                       values.dirichletPlane.planeZ)
                                 : directSum(before, kCharge, 2.0e-5);
                EXPECT_LT(rmsRelativeError(after, reference), 1.0e-3);
            }
            // Dipole multipoles only while the images are active, quadrupoles after expiry.
            EXPECT_EQ(run.algorithm->statistics().solves, 4u);
            EXPECT_EQ(run.algorithm->statistics().dipoleSolves, 2u);
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
                // Limits ~2.5x the measured errors (7.28e-4, 1.73e-4, 1.64e-5 on 1 and 4 ranks;
                // 1.49e-3, 3.60e-4, 3.99e-5 with the cubic domain), so an accuracy regression of
                // BH against the direct sum fails this test.
                const double limit = theta >= 0.7 ? 2.0e-3 : (theta >= 0.5 ? 5.0e-4 : 5.0e-5);
                EXPECT_LT(error, limit);
                EXPECT_LT(error, previous);
                previous = error;
            }
        }

        TEST_F(BarnesHutAlgorithmTest, StretchedBunchesUseMixedDimensionKeys) {
            // The Barnes-Hut domain is the bunch's bounding box, not a cube: for a pancake
            // (drift4m-like, 1 mm x 1 mm x 10 um) and a needle the short axes get fewer SFC key
            // bits, so the octree cells stay close to cubic. Accuracy must stay at the usual BH
            // level against the direct sum.
            struct Case {
                const char* name;
                Vector sigma;
            };
            for (const Case& c : {Case{"pancake", Vector(1.0e-3, 1.0e-3, 1.0e-5)},
                                  Case{"needle", Vector(1.0e-5, 1.0e-5, 1.0e-3)}}) {
                SCOPED_TRACE(c.name);
                BunchSpec spec;
                spec.globalCount = 4000;
                spec.variableDt  = true;
                spec.sigma       = c.sigma;
                Run run(config(1.0e-6), spec);
                const auto before = gatherAll(run.particles);
                ASSERT_EQ(run.solve().backendSolves, 1u);
                const auto after = gatherAll(run.particles);
                expectConserved(before, after);
                const double error = rmsRelativeError(after, directSum(before, kCharge, 1.0e-6));
                const auto& bits   = run.algorithm->statistics().keyBits;
                if (ippl::Comm->rank() == 0) {
                    std::printf("[BH] %s key bits (%u, %u, %u) rms relative error=%.3e\n", c.name,
                                bits[0], bits[1], bits[2], error);
                }
                EXPECT_EQ(run.algorithm->statistics().solves, 1u);
                // sigma ratio 100: about log2(100) ~ 6-7 key bits fewer on each short axis; the
                // long axes keep (nearly) all of them.
                const unsigned most = std::max({bits[0], bits[1], bits[2]});
                for (unsigned d = 0; d < 3; ++d) {
                    if (c.sigma[d] < 1.0e-4) {
                        EXPECT_LE(bits[d] + 5, most) << "axis " << d;
                    } else {
                        EXPECT_GE(bits[d] + 1, most) << "axis " << d;
                    }
                }
                // Measured (1 and 4 ranks): pancake bits (21, 21, 14) 1.05e-4, needle (14, 14, 21)
                // 6.49e-5.
                EXPECT_LT(error, 3.0e-4);
            }
        }

        TEST_F(BarnesHutAlgorithmTest, StretchedBunchesOverManySolves) {
            // drift4m: a 1000-particle pancake solved step after step on several ranks. The first
            // sync converges cstone's trees from scratch, later ones update them incrementally on
            // the mixed-dimension keys.
            struct Case {
                const char* name;
                Vector sigma;
            };
            for (const Case& c : {Case{"pancake", Vector(1.0e-3, 1.0e-3, 1.0e-5)},
                                  Case{"needle", Vector(1.0e-5, 1.0e-5, 1.0e-3)}}) {
                SCOPED_TRACE(c.name);
                BunchSpec spec;
                spec.globalCount = 1000;
                spec.sigma       = c.sigma;
                Run run(config(1.0e-6), spec);
                const std::size_t steps = 20;
                for (std::size_t step = 0; step < steps; ++step) {
                    SCOPED_TRACE(step);
                    driftById(run.particles, 2.0e-6);
                    ASSERT_EQ(run.solve().backendSolves, 1u);
                    ASSERT_EQ(run.particles.getTotalNum(), spec.globalCount);
                }
                const auto& stats = run.algorithm->statistics();
                EXPECT_EQ(stats.solves, steps);
                const auto records = gatherAll(run.particles);
                const double error = rmsRelativeError(records, directSum(records, kCharge, 1.0e-6));
                if (ippl::Comm->rank() == 0) {
                    std::printf("[BH] %s over %zu solves: key bits (%u, %u, %u) rms relative "
                                "error=%.3e\n",
                                c.name, steps, stats.keyBits[0], stats.keyBits[1],
                                stats.keyBits[2], error);
                }
                // Measured (1 / 4 ranks): pancake 4.3e-5 / 5.9e-5, needle 3.4e-5 / 7.6e-5.
                EXPECT_LT(error, 3.0e-4);
            }
        }

        TEST_F(BarnesHutAlgorithmTest, DriftingPancakeKeepsItsKeyLayout) {
            // drift4m: the bunch moves 18 mm per step, ~1000 times its length, through a fixed
            // solve frame. The box moves along at its size and keeps the domain and key bits.
            BunchSpec spec;
            spec.globalCount = 1000;
            spec.sigma       = Vector(1.0e-3, 1.0e-3, 1.0e-5);
            spec.center      = Vector(0.0, 0.0, 0.009);
            Run run(config(1.0e-6), spec);
            const std::size_t steps = 10;
            std::array<unsigned, 3> first{};
            for (std::size_t step = 0; step < steps; ++step) {
                SCOPED_TRACE(step);
                if (step > 0) {
                    transformPositions(run.particles, [](Vector& r) { r[2] += 0.018; });
                }
                ASSERT_EQ(run.solve().backendSolves, 1u);
                if (step == 0) {
                    first = run.algorithm->statistics().keyBits;
                    EXPECT_GE(first[0], first[2] + 6);
                }
                EXPECT_EQ(run.algorithm->statistics().keyBits, first);
            }
            const auto& stats = run.algorithm->statistics();
            EXPECT_EQ(stats.solves, steps);
            EXPECT_EQ(stats.containerBuilds, 1u);
            EXPECT_EQ(stats.boxRebuilds, 0u);
            EXPECT_EQ(stats.boxMoves, steps - 1);
            const auto records = gatherAll(run.particles);
            const double error = rmsRelativeError(records, directSum(records, kCharge, 1.0e-6));
            if (ippl::Comm->rank() == 0) {
                std::printf("[BH] drifting pancake: rms relative error=%.3e\n", error);
            }
            EXPECT_LT(error, 3.0e-4);  // measured 8.5e-5 (1 rank) / 1.16e-4 (4 ranks)
        }

        TEST_F(BarnesHutAlgorithmTest, BoxFollowsADriftWithinItsMargin) {
            // A drift within the margin (SwissFEL booster, BHBOXMARGIN 0.5) must still move the box
            // at every solve; a resting box let keys slide and left overfull leaves.
            BunchSpec spec;
            spec.globalCount = 1000;
            spec.sigma       = Vector(1.0e-3, 1.0e-3, 1.0e-5);
            spec.center      = Vector(0.0, 0.0, 0.009);
            auto values      = config(1.0e-6);
            values.boxMargin = 0.5;  // ~3 sigma of room per side, ~15 drift steps
            Run run(values, spec);
            const std::size_t steps = 10;
            for (std::size_t step = 0; step < steps; ++step) {
                SCOPED_TRACE(step);
                if (step > 0) {
                    transformPositions(run.particles, [](Vector& r) { r[2] += 2.0e-6; });
                }
                ASSERT_EQ(run.solve().backendSolves, 1u);
            }
            const auto& stats = run.algorithm->statistics();
            EXPECT_EQ(stats.containerBuilds, 1u);
            EXPECT_EQ(stats.boxRebuilds, 0u);
            EXPECT_EQ(stats.boxMoves, steps - 1);
            const auto records = gatherAll(run.particles);
            const double error = rmsRelativeError(records, directSum(records, kCharge, 1.0e-6));
            if (ippl::Comm->rank() == 0) {
                std::printf("[BH] drift within the box margin: %zu moves, rms relative error=%.3e\n",
                            stats.boxMoves, error);
            }
            EXPECT_LT(error, 3.0e-4);  // as DriftingPancakeKeepsItsKeyLayout
        }

        TEST_F(BarnesHutAlgorithmTest, LengtheningPancakeRebuildsItsDomain) {
            // The bunch stretches from 10 um to 1 mm in z: the box grows, the key bits per axis
            // change, and every resize rebuilds the domain.
            BunchSpec spec;
            spec.globalCount = 1000;
            spec.sigma       = Vector(1.0e-3, 1.0e-3, 1.0e-5);
            Run run(config(1.0e-6), spec);
            const double centreZ   = spec.center[2];
            const std::size_t steps = 15;
            std::array<unsigned, 3> first{};
            for (std::size_t step = 0; step < steps; ++step) {
                SCOPED_TRACE(step);
                if (step > 0) {
                    // x1.39 per step: 100 times longer after 14 steps.
                    transformPositions(run.particles, [centreZ](Vector& r) {
                        r[2] = centreZ + 1.389 * (r[2] - centreZ);
                    });
                }
                ASSERT_EQ(run.solve().backendSolves, 1u);
                if (step == 0) {
                    first = run.algorithm->statistics().keyBits;
                }
            }
            const auto& stats = run.algorithm->statistics();
            EXPECT_EQ(stats.solves, steps);
            EXPECT_GE(first[0], first[2] + 6);
            EXPECT_GE(stats.keyBits[2] + 1, stats.keyBits[0]);
            EXPECT_EQ(stats.boxRebuilds, steps - 1);  // outgrows the box every step
            const auto records = gatherAll(run.particles);
            const double error = rmsRelativeError(records, directSum(records, kCharge, 1.0e-6));

            // A rebuilt domain must be as accurate as a new one for the final bunch.
            Run fresh(config(1.0e-6), spec);
            transformPositions(fresh.particles, [centreZ](Vector& r) {
                r[2] = centreZ + std::pow(1.389, 14) * (r[2] - centreZ);
            });
            (void)fresh.solve();
            const auto freshRecords = gatherAll(fresh.particles);
            const double freshError =
                    rmsRelativeError(freshRecords, directSum(freshRecords, kCharge, 1.0e-6));
            if (ippl::Comm->rank() == 0) {
                std::printf("[BH] lengthening pancake: key bits (%u, %u, %u) -> (%u, %u, %u), %zu "
                            "rebuilds, rms relative error=%.3e (fresh domain %.3e)\n",
                            first[0], first[1], first[2], stats.keyBits[0], stats.keyBits[1],
                            stats.keyBits[2], stats.boxRebuilds, error, freshError);
            }
            EXPECT_LT(std::abs(error - freshError), 0.25 * freshError);
            // Measured 3.28e-4 (1 rank) / 3.99e-4 (4 ranks): 1000 particles in a round bunch.
            EXPECT_LT(error, 8.0e-4);
        }

        TEST_F(BarnesHutAlgorithmTest, MixedPrecisionMatchesTheReference) {
            // BHPRECISION=MIXED: double coordinates (R, P, dt and ID come back exactly), float
            // charges, softening, multipoles and accelerations. The field error must stay at the
            // tree level, without and with the rest frame and near image charges (second pass).
            BunchSpec spec;
            spec.globalCount          = 4000;
            spec.variableDt           = true;
            spec.longitudinalMomentum = std::sqrt(3.0);
            const double planeZ       = spec.center[2] - 2.0 * spec.sigma[2];
            for (bool restFrame : {false, true}) {
                SCOPED_TRACE(restFrame);
                auto values      = restFrame ? restFrameConfig() : config(2.0e-5);
                values.precision = BarnesHutPrecision::Mixed;
                if (restFrame) {
                    values.dirichletPlane = {.kind = DirichletPlaneType::ShiftedGreen, .planeZ = planeZ};
                }
                Run run(values, spec);
                const auto before = gatherAll(run.particles);
                EXPECT_EQ(run.solve().backendSolves, 1u);
                const auto after = gatherAll(run.particles);
                double errorE = 0.0, errorB = 0.0;
                if (restFrame) {
                    expectConservedUpToStretch(before, after);
                    const auto reference = restFrameDirectSum(before, kCharge, 2.0e-5, planeZ);
                    errorE = rmsRelativeDifference(fieldsOf(after, &Record::e), reference.e);
                    errorB = rmsRelativeDifference(fieldsOf(after, &Record::b), reference.b);
                    EXPECT_EQ(run.algorithm->statistics().imagePasses, 1u);
                } else {
                    expectConserved(before, after);
                    errorE = rmsRelativeError(after, directSum(before, kCharge, 2.0e-5));
                }
                if (ippl::Comm->rank() == 0) {
                    std::printf("[BH] mixed precision%s: rms relative error E=%.3e B=%.3e\n",
                                restFrame ? " (rest frame, near images)" : "", errorE, errorB);
                }
                // Double precision measured 1.4e-4 / 1.7e-4 for these bunches; float adds ~1e-6.
                EXPECT_LT(errorE, 5.0e-4);
                EXPECT_LT(errorB, 5.0e-4);
            }
        }

        TEST_F(BarnesHutAlgorithmTest, BoxMarginAndSlackSetTheRebuildCount) {
            // Bunch length x1.389 (or /1.389) per solve, 14 solves. Growth rebuilds: every step for
            // margin 0.05, every third for 0.5; shrink rebuilds: every step for slack 0.25, every
            // third for 1.0. Accuracy as for a fresh domain.
            struct Case {
                double margin, slack, factor;
                std::size_t growth, shrink;
            };
            for (const Case& c : {Case{0.05, 0.25, 1.389, 14, 0}, Case{0.5, 0.25, 1.389, 4, 0},
                                  Case{0.05, 0.25, 1.0 / 1.389, 0, 14},
                                  Case{0.05, 1.0, 1.0 / 1.389, 0, 4}}) {
                SCOPED_TRACE(c.margin);
                SCOPED_TRACE(c.slack);
                SCOPED_TRACE(c.factor);
                BunchSpec spec;
                spec.globalCount = 1000;
                spec.sigma       = Vector(1.0e-3, 1.0e-3, c.factor > 1.0 ? 1.0e-5 : 1.0e-3);
                auto values      = config(1.0e-6);
                values.boxMargin = c.margin;
                values.boxSlack  = c.slack;
                Run run(values, spec);
                const double centreZ = spec.center[2];
                for (std::size_t step = 0; step < 15; ++step) {
                    if (step > 0) {
                        transformPositions(run.particles, [centreZ, &c](Vector& r) {
                            r[2] = centreZ + c.factor * (r[2] - centreZ);
                        });
                    }
                    ASSERT_EQ(run.solve().backendSolves, 1u);
                }
                const auto& stats = run.algorithm->statistics();
                EXPECT_EQ(stats.boxRebuilds - stats.boxShrinkRebuilds, c.growth);
                EXPECT_EQ(stats.boxShrinkRebuilds, c.shrink);
                EXPECT_EQ(stats.containerBuilds, 1 + c.growth + c.shrink);
                EXPECT_NE(run.algorithm->runSummary().find("BHBOXMARGIN"), std::string::npos);
                const auto records = gatherAll(run.particles);
                const double error = rmsRelativeError(records, directSum(records, kCharge, 1.0e-6));
                if (ippl::Comm->rank() == 0) {
                    std::printf("[BH] box margin %.2f slack %.2f, length x%.3f per solve: %zu "
                                "growth + %zu shrink rebuilds, rms relative error=%.3e\n",
                                c.margin, c.slack, c.factor,
                                stats.boxRebuilds - stats.boxShrinkRebuilds,
                                stats.boxShrinkRebuilds, error);
                }
                // Same accuracy as LengtheningPancakeRebuildsItsDomain (3.3e-4 / 4.0e-4).
                EXPECT_LT(error, 8.0e-4);
            }
        }

        TEST_F(BarnesHutAlgorithmTest, DriftingPancakeLosesAnOutlier) {
            // drift4m step 25: a far outlier is lost and the box shrinks to the remaining bunch,
            // which rebuilds the domain.
            BunchSpec spec;
            spec.globalCount = 1000;
            spec.sigma       = Vector(1.0e-3, 1.0e-3, 1.0e-5);
            spec.center      = Vector(0.0, 0.0, 0.009);
            Run run(config(1.0e-6), spec);
            {
                // Global particle 0 (rank 0) sits 8 sigma below the bunch in y.
                auto R  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), run.particles.R.getView());
                auto ID = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), run.particles.ID.getView());
                for (std::size_t i = 0; i < run.particles.getLocalNum(); ++i) {
                    if (ID(i) == 0) {
                        R(i)[1] = -8.0e-3;
                    }
                }
                Kokkos::deep_copy(run.particles.R.getView(), R);
                run.particles.markMomentsDirty();
            }
            const std::size_t steps = 10;
            for (std::size_t step = 0; step < steps; ++step) {
                SCOPED_TRACE(step);
                if (step > 0) {
                    transformPositions(run.particles, [](Vector& r) { r[2] += 0.018; });
                }
                if (step == 5) {
                    auto ID      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), run.particles.ID.getView());
                    auto invalid = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), run.particles.InvalidMask.getView());
                    for (std::size_t i = 0; i < run.particles.getLocalNum(); ++i) {
                        invalid(i) = ID(i) == 0;
                    }
                    Kokkos::deep_copy(run.particles.InvalidMask.getView(), invalid);
                    (void)run.particles.deleteInvalidParticles();
                    ASSERT_EQ(run.particles.getTotalNum(), spec.globalCount - 1);
                }
                ASSERT_EQ(run.solve().backendSolves, 1u);
            }
            const auto& stats = run.algorithm->statistics();
            EXPECT_EQ(stats.solves, steps);
            EXPECT_EQ(stats.boxRebuilds, 1u);  // after the outlier is gone
            const auto records = gatherAll(run.particles);
            const double error = rmsRelativeError(records, directSum(records, kCharge, 1.0e-6));
            if (ippl::Comm->rank() == 0) {
                std::printf("[BH] drifting pancake losing an outlier: %zu box rebuilds, %zu moves, "
                            "rms relative error=%.3e\n",
                            stats.boxRebuilds, stats.boxMoves, error);
            }
            EXPECT_LT(error, 3.0e-4);  // measured 8.5e-5 (1 rank) / 1.16e-4 (4 ranks)
        }

        TEST_F(BarnesHutAlgorithmTest, ZeroThicknessSheetMatchesDirectSum) {
            // All particles in one plane z = const: the zero-extent axis is widened to the
            // smallest extent the keys can carry (one key bit) instead of breaking the box.
            BunchSpec spec;
            spec.globalCount = 4000;
            spec.sigma       = Vector(1.0e-3, 1.5e-3, 0.0);
            Run run(config(2.0e-5), spec);
            const auto before = gatherAll(run.particles);
            for (const Record& r : before) {
                ASSERT_EQ(r.r[2], spec.center[2]);
            }
            ASSERT_EQ(run.solve().backendSolves, 1u);
            const auto after = gatherAll(run.particles);
            expectConserved(before, after);
            const double error = rmsRelativeError(after, directSum(before, kCharge, 2.0e-5));
            const auto& bits   = run.algorithm->statistics().keyBits;
            if (ippl::Comm->rank() == 0) {
                std::printf("[BH] sheet key bits (%u, %u, %u) rms relative error=%.3e\n", bits[0],
                            bits[1], bits[2], error);
            }
            EXPECT_EQ(run.algorithm->statistics().solves, 1u);
            EXPECT_EQ(bits[2], 1u);
            // In-plane sources only: the normal component vanishes up to rounding of the
            // expansion centres.
            double normal = 0.0, total = 0.0;
            for (const Record& r : after) {
                normal += r.e[2] * r.e[2];
                total += r.e[0] * r.e[0] + r.e[1] * r.e[1] + r.e[2] * r.e[2];
            }
            EXPECT_LT(std::sqrt(normal / total), 1.0e-8);
            // Measured 2.79e-4 (4 ranks) since cstone splits target groups in mixed-dimension
            // boxes (2026-10-16); before, groups were never split there (64 SFC-consecutive
            // particles, wide boxes, mostly P2P) and the error was 8.0e-5 / 7.9e-5.
            EXPECT_LT(error, 7.0e-4);
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

            // The CIC/FFT reference converges to the Coulomb field as the mesh is refined, so the
            // BH-vs-OPEN difference must shrink with resolution, not merely stay small.
            std::vector<double> differences;
            for (std::size_t mesh : {32u, 64u}) {
            CartesianPIC3DConfig pic;
            pic.backend            = PoissonSolverType::Open;
            pic.grid.meshSize      = {mesh, mesh, mesh};
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
                    std::printf("[BH] BH vs OPEN FFT (%zu^3) rms relative difference=%.3e\n",
                                mesh, difference);
                }
                differences.push_back(difference);
            }
            }
            EXPECT_LT(differences[1], differences[0]) << "the gap must shrink with the mesh size";
            EXPECT_LT(differences[1], 5.0e-2);
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
                const double error = rmsRelativeError(after, reference);
                if (ippl::Comm->rank() == 0) {
                    std::printf("[BH] grow from one rank: rms relative error=%.3e\n", error);
                }
                EXPECT_LT(error, 8.0e-4);  // measured 1.73e-4 / 3.25e-4
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

        TEST_F(BarnesHutAlgorithmTest, RepackEveryStepMatchesInPlace) {
            // The destroy + create path on a container that already has halos (taken whenever the
            // local count changes, e.g. during emission) must give the same fields as the in-place
            // path that IncrementalMigrationOverManySolves exercises.
            BunchSpec spec;
            spec.globalCount   = 20000;
            auto repackConfig  = config();
            repackConfig.repackEverySolve = true;
            Run inPlace(config(), spec), repacked(repackConfig, spec);
            for (std::size_t step = 0; step < 30; ++step) {
                SCOPED_TRACE(step);
                driftById(inPlace.particles, 2.0e-6);
                driftById(repacked.particles, 2.0e-6);
                ASSERT_EQ(inPlace.solve().backendSolves, 1u);
                ASSERT_EQ(repacked.solve().backendSolves, 1u);
                const auto a = gatherAll(inPlace.particles);
                const auto b = gatherAll(repacked.particles);
                ASSERT_EQ(a.size(), spec.globalCount);
                ASSERT_EQ(b.size(), spec.globalCount);
                double scale = 0.0, worst = 0.0;
                for (std::size_t i = 0; i < a.size(); ++i) {
                    ASSERT_EQ(a[i].id, b[i].id) << "particle sets differ";
                    scale = std::max(scale, norm(a[i].e));
                    for (unsigned d = 0; d < 3; ++d) {
                        ASSERT_EQ(a[i].r[d], b[i].r[d]);
                        ASSERT_TRUE(std::isfinite(b[i].e[d]));
                        worst = std::max(worst, std::abs(a[i].e[d] - b[i].e[d]));
                    }
                }
                // Not bit-identical: destroyed slots still enter cstone's bounding box, so the trees
                // differ slightly. Deviations must stay far below the BH error (~3.6e-4 rms).
                ASSERT_LT(worst, 1.0e-5 * scale) << "repacked fields deviate from in-place fields";
            }
            EXPECT_EQ(repacked.algorithm->statistics().inPlacePacks, 0u);
            EXPECT_GT(inPlace.algorithm->statistics().inPlacePacks, 0u);
        }

        TEST_F(BarnesHutAlgorithmTest, GrowsEveryStepLikeEmission) {
            // Emission-like: every rank adds particles at a thin "cathode" layer each step while the
            // bunch drifts away from it, so every solve takes the destroy + create path on a
            // container with halos and the bounding box keeps growing.
            BunchSpec spec;
            spec.globalCount = 2000;
            Run run(config(), spec);
            BunchSpec layer   = spec;
            layer.sigma       = Vector(1.0e-3, 1.5e-3, 1.0e-6);
            layer.center      = Vector(0.0, 0.0, spec.center[2] - 4.0 * spec.sigma[2]);
            const std::size_t rank = ippl::Comm->rank(), size = ippl::Comm->size();
            const std::size_t perStepGlobal = 40;
            std::uint64_t nextGlobal = spec.globalCount;
            for (std::size_t step = 0; step < 40; ++step) {
                SCOPED_TRACE(step);
                driftById(run.particles, 2.0e-5);
                // Round-robin share of this step's new particles, with fractional dt like emission.
                std::vector<std::uint64_t> mine;
                for (std::uint64_t k = nextGlobal + rank; k < nextGlobal + perStepGlobal; k += size) {
                    mine.push_back(k);
                }
                nextGlobal += perStepGlobal;
                const std::size_t offset = run.particles.getLocalNum();
                run.particles.createParticles(mine.size());
                auto R  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), run.particles.R.getView());
                auto dt = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), run.particles.dt.getView());
                for (std::size_t i = 0; i < mine.size(); ++i) {
                    R(offset + i)  = positionFor(layer, mine[i]);
                    dt(offset + i) = kTimeStep * (0.05 + 0.95 * uniform(mine[i], 41));
                }
                Kokkos::deep_copy(run.particles.R.getView(), R);
                Kokkos::deep_copy(run.particles.dt.getView(), dt);

                const auto before = gatherAll(run.particles);
                ASSERT_EQ(run.solve().backendSolves, 1u);
                ASSERT_EQ(run.particles.getTotalNum(), nextGlobal);
                const auto after = gatherAll(run.particles);
                expectConserved(before, after);
                for (const Record& r : after) {
                    ASSERT_TRUE(std::isfinite(r.e[0]) && std::isfinite(r.e[1]) && std::isfinite(r.e[2]));
                }
                if (step % 10 == 9) {
                    const double error = rmsRelativeError(after, directSum(before, kCharge, 1.0e-6));
                    if (ippl::Comm->rank() == 0) {
                        std::printf("[BH] emission-like step %zu N=%zu rms relative error=%.3e\n",
                                    step, static_cast<std::size_t>(nextGlobal), error);
                    }
                    ASSERT_LT(error, 4.0e-4);  // measured 5.4e-6 - 1.26e-4
                }
            }
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

        TEST_F(BarnesHutAlgorithmTest, RestFrameAtRestMatchesElectrostatic) {
            // A bunch at rest has gamma = 1 and v = 0: the rest-frame solve must reproduce the
            // electrostatic solve bit for bit, with B = 0, on the direct-sum and the tree path.
            for (std::size_t count : {48u, 2000u}) {
                SCOPED_TRACE(count);
                BunchSpec spec;
                spec.globalCount          = count;
                spec.transverseMomentum   = 0.0;
                spec.longitudinalMomentum = 0.0;
                Run plain(config(2.0e-5), spec), rest(restFrameConfig(), spec);
                EXPECT_EQ(plain.solve().backendSolves, 1u);
                EXPECT_EQ(rest.solve().backendSolves, 1u);
                const auto a = gatherAll(plain.particles);
                const auto b = gatherAll(rest.particles);
                ASSERT_EQ(a.size(), b.size());
                for (std::size_t i = 0; i < a.size(); ++i) {
                    ASSERT_EQ(a[i].id, b[i].id);
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_EQ(a[i].r[d], b[i].r[d]);
                        EXPECT_EQ(a[i].e[d], b[i].e[d]);
                        EXPECT_EQ(b[i].b[d], 0.0);
                    }
                }
                EXPECT_EQ(rest.algorithm->statistics().restFrameSolves, 1u);
                EXPECT_EQ(rest.algorithm->statistics().gamma, 1.0);
                EXPECT_EQ(plain.algorithm->statistics().restFrameSolves, 0u);
            }
        }

        TEST_F(BarnesHutAlgorithmTest, RestFrameMatchesBoostedDirectSum) {
            // gamma ~ 5 with momentum spread: E and B against the boosted rest-frame direct sum,
            // on the direct-sum path (small bunch and BHDIRECT, exact kernel) and the tree path.
            BunchSpec spec;
            spec.transverseMomentum   = 0.05;
            spec.longitudinalMomentum = std::sqrt(24.0);
            spec.momentumSpread       = 0.05;
            spec.variableDt           = true;
            for (std::size_t count : {48u, 4000u}) {
                for (bool forcedDirect : {false, true}) {
                    if (count == 48u && forcedDirect) {
                        continue;
                    }
                    SCOPED_TRACE(count);
                    SCOPED_TRACE(forcedDirect);
                    spec.globalCount = count;
                    auto values      = restFrameConfig();
                    values.directSum = forcedDirect;
                    Run run(values, spec);
                    const auto before = gatherAll(run.particles);
                    EXPECT_EQ(run.solve().backendSolves, 1u);
                    const auto after = gatherAll(run.particles);
                    expectConservedUpToStretch(before, after);
                    const auto reference = restFrameDirectSum(before, kCharge, 2.0e-5);
                    const double errorE  = rmsRelativeDifference(fieldsOf(after, &Record::e), reference.e);
                    const double errorB  = rmsRelativeDifference(fieldsOf(after, &Record::b), reference.b);
                    const bool direct =
                            forcedDirect || count < 64 * static_cast<std::size_t>(ippl::Comm->size());
                    if (ippl::Comm->rank() == 0) {
                        std::printf("[BH] rest frame gamma=%.4f N=%zu %s rms relative error "
                                    "E=%.3e B=%.3e\n",
                                    reference.gamma, count, direct ? "direct" : "tree", errorE,
                                    errorB);
                    }
                    EXPECT_NEAR(run.algorithm->statistics().gamma, reference.gamma,
                                1.0e-12 * reference.gamma);
                    EXPECT_EQ(run.algorithm->statistics().restFrameSolves, 1u);
                    // Tree measured E 1.41e-4, B 1.39e-4 (1 and 4 ranks).
                    EXPECT_LT(errorE, direct ? 1.0e-10 : 5.0e-4);
                    EXPECT_LT(errorB, direct ? 1.0e-10 : 5.0e-4);
                }
            }
        }

        TEST_F(BarnesHutAlgorithmTest, RestFrameImagesNearAndFar) {
            // gamma ~ 2. Near plane (2 sigma_z behind): images are tree sources, a second pass gives
            // their field. Far plane (50 mm): one quadrupole, checked against BHDIRECT. Measured:
            // near E 1.67e-4, B 1.58e-4; far E 1.21e-4, B 1.17e-4; quadrupole error 1.4e-5 (E),
            // 1.4e-4 (B) of the image field.
            BunchSpec spec;
            spec.globalCount          = 4000;
            spec.variableDt           = true;
            spec.longitudinalMomentum = std::sqrt(3.0);
            for (double planeZ : {spec.center[2] - 2.0 * spec.sigma[2], spec.center[2] - 0.05}) {
                const bool far = planeZ < spec.center[2] - 0.01;
                for (bool forcedDirect : {true, false}) {
                    SCOPED_TRACE(planeZ);
                    SCOPED_TRACE(forcedDirect);
                    auto values           = restFrameConfig();
                    values.directSum      = forcedDirect;
                    values.dirichletPlane = {.kind = DirichletPlaneType::ShiftedGreen, .planeZ = planeZ};
                    Run run(values, spec);
                    const auto before = gatherAll(run.particles);
                    EXPECT_EQ(run.solve().backendSolves, 1u);
                    const auto after = gatherAll(run.particles);
                    expectConservedUpToStretch(before, after);
                    const auto reference = restFrameDirectSum(before, kCharge, 2.0e-5, planeZ);
                    const auto unmirrored = restFrameDirectSum(before, kCharge, 2.0e-5);
                    const Fields e = fieldsOf(after, &Record::e), b = fieldsOf(after, &Record::b);
                    const double errorE = rmsRelativeDifference(e, reference.e);
                    const double errorB = rmsRelativeDifference(b, reference.b);
                    // Share of the image charges in the reference fields.
                    const double imageE = rmsRelativeDifference(unmirrored.e, reference.e);
                    const double imageB = rmsRelativeDifference(unmirrored.b, reference.b);
                    if (ippl::Comm->rank() == 0) {
                        std::printf("[BH] rest-frame images %s plane %s: rms relative error "
                                    "E=%.3e B=%.3e (image share E=%.3e B=%.3e)\n",
                                    far ? "far" : "near", forcedDirect ? "direct" : "tree", errorE,
                                    errorB, imageE, imageB);
                    }
                    const auto& stats = run.algorithm->statistics();
                    EXPECT_EQ(stats.farImageSolves, far ? 1u : 0u);
                    EXPECT_EQ(stats.imagePasses, far ? 0u : 1u);
                    EXPECT_EQ(stats.dipoleSolves, (far || forcedDirect) ? 0u : 1u);
                    EXPECT_GT(imageE, far ? 1.0e-4 : 5.0e-2);  // the plane must matter
                    if (forcedDirect && !far) {
                        EXPECT_LT(errorE, 1.0e-10);
                        EXPECT_LT(errorB, 1.0e-10);
                    } else if (forcedDirect) {
                        EXPECT_LT(errorE, 1.0e-3 * imageE);
                        EXPECT_LT(errorB, 1.0e-3 * imageB);
                    } else {
                        EXPECT_LT(errorE, 5.0e-4);
                        EXPECT_LT(errorB, 5.0e-4);
                    }
                }
            }
        }

        TEST_F(BarnesHutAlgorithmTest, RestFrameAgreesWithBinnedOpenFFT) {
            // BH and binned OPEN (one bin) against the boosted continuum field of the Gaussian (and
            // its image behind a shifted-Green plane): both within 2-3 % (graininess, sampling); a
            // wrong gamma, velocity or B sign would be O(1).
            BunchSpec spec;
            spec.globalCount          = 200000;
            spec.sigma                = Vector(1.0e-3, 1.0e-3, 2.0e-3);
            spec.transverseMomentum   = 0.05;
            spec.longitudinalMomentum = std::sqrt(24.0);
            spec.momentumSpread       = 0.02;
            const double planeZ       = spec.center[2] - 2.0 * spec.sigma[2];
            for (bool withPlane : {false, true}) {
                SCOPED_TRACE(withPlane);
                const DirichletPlaneConfig plane =
                        withPlane ? DirichletPlaneConfig{.kind   = DirichletPlaneType::ShiftedGreen,
                                                         .planeZ = planeZ}
                                  : DirichletPlaneConfig{};
                // h = 1e-4 m is about one rest-frame mean spacing (9.2e-5 m) of this bunch.
                auto values           = restFrameConfig(1.0e-4);
                values.dirichletPlane = plane;
                Run bh(values, spec);
                const auto particles = gatherAll(bh.particles);
                (void)bh.solve();
                const auto continuum = gaussianContinuum(
                        particles, spec, kCharge * spec.globalCount,
                        withPlane ? std::optional<double>(planeZ) : std::nullopt);
                const auto bhRecords  = gatherAll(bh.particles);
                const auto picRecords = binnedOpenFields(spec, 64, plane);
                ASSERT_EQ(bhRecords.size(), particles.size());
                ASSERT_EQ(picRecords.size(), particles.size());
                double errors[2][2];  // [BH, OPEN][E, B]
                const std::vector<Record>* records[2] = {&bhRecords, &picRecords};
                for (unsigned m = 0; m < 2; ++m) {
                    for (std::size_t i = 0; i < particles.size(); ++i) {
                        ASSERT_EQ((*records[m])[i].id, particles[i].id);
                    }
                    errors[m][0] = rmsRelativeDifference(fieldsOf(*records[m], &Record::e), continuum.e);
                    errors[m][1] = rmsRelativeDifference(fieldsOf(*records[m], &Record::b), continuum.b);
                }
                if (ippl::Comm->rank() == 0) {
                    std::printf("[BH] rest frame vs continuum (%s): BH h=1e-4 E=%.3e B=%.3e, binned "
                                "OPEN 64^3 E=%.3e B=%.3e\n",
                                withPlane ? "shifted Green" : "no plane", errors[0][0], errors[0][1],
                                errors[1][0], errors[1][1]);
                }
                // Without plane measured BH 3.17e-2, OPEN 1.91e-2 (diagnostic test).
                for (unsigned m = 0; m < 2; ++m) {
                    EXPECT_LT(errors[m][0], 5.0e-2);
                    EXPECT_LT(errors[m][1], 5.0e-2);
                }
            }
        }

        TEST_F(BarnesHutAlgorithmTest, RelativeSofteningFollowsTheRestFrameSpacing) {
            // BHSOFTENINGREL must equal BHSOFTENING with h = factor x rest-frame mean spacing: a
            // Gaussian at rest, at gamma ~ 5 (spacing x gamma^(1/3)), and a sheet (2D spacing).
            struct Case {
                const char* name;
                double pz, thickness;
            };
            for (const Case& c : {Case{"at rest", 0.0, 3.0e-3}, Case{"gamma 5", std::sqrt(24.0), 3.0e-3},
                                  Case{"sheet", 0.0, 0.0}}) {
                SCOPED_TRACE(c.name);
                BunchSpec spec;
                spec.globalCount          = 4000;
                spec.sigma                = Vector(1.0e-3, 1.5e-3, c.thickness);
                spec.transverseMomentum   = 0.0;
                spec.longitudinalMomentum = c.pz;
                auto relative              = c.pz > 0.0 ? restFrameConfig() : config(2.0e-5);
                relative.softening         = 0.0;
                relative.relativeSoftening = 0.5;
                Run run(relative, spec);
                const auto before = gatherAll(run.particles);
                (void)run.solve();
                const double h = run.algorithm->statistics().softening;

                // Host reference of the spacing.
                double v[3];
                const double gamma = meanRestFrame(before, v);
                const double n     = static_cast<double>(before.size());
                std::array<double, 3> edge{};
                for (unsigned d = 0; d < 3; ++d) {
                    const double scale = d == 2 ? gamma : 1.0;
                    double sum = 0.0, sum2 = 0.0;
                    for (const Record& r : before) {
                        sum += scale * r.r[d];
                        sum2 += scale * r.r[d] * scale * r.r[d];
                    }
                    const double mean = sum / n;
                    edge[d] = std::sqrt(2.0 * Physics::pi * std::max(0.0, sum2 / n - mean * mean));
                }
                std::sort(edge.begin(), edge.end());
                double spacing = 0.0;
                for (unsigned thin = 0; thin < 3; ++thin) {
                    double volume = 1.0;
                    for (unsigned d = thin; d < 3; ++d) {
                        volume *= edge[d];
                    }
                    spacing = std::pow(volume / n, 1.0 / (3.0 - thin));
                    if (edge[thin] > 0.0 && edge[thin] >= spacing) {
                        break;
                    }
                }
                if (ippl::Comm->rank() == 0) {
                    std::printf("[BH] relative softening %s: spacing %.3e m, h = %.3e m\n", c.name,
                                spacing, h);
                }
                EXPECT_NEAR(h, 0.5 * spacing, 1.0e-9 * spacing);

                auto fixed      = relative;
                fixed.softening = h;
                fixed.relativeSoftening = 0.0;
                Run reference(fixed, spec);
                (void)reference.solve();
                const auto a = gatherAll(run.particles);
                const auto b = gatherAll(reference.particles);
                ASSERT_EQ(a.size(), b.size());
                for (std::size_t i = 0; i < a.size(); ++i) {
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_EQ(a[i].e[d], b[i].e[d]);
                        EXPECT_EQ(a[i].b[d], b[i].b[d]);
                    }
                }
            }
        }

        TEST_F(BarnesHutAlgorithmTest, DISABLED_SolveTimeVsSoftening) {
            // Diagnostic: 262144 particles, h from 1/4 to 4 spacings, leaf-based h, both
            // precisions; the tree does not depend on h. BHSTATS prints the counts per solve.
            BunchSpec spec;
            spec.globalCount   = 262144;
            spec.sigma         = Vector(1.0e-3, 1.0e-3, 3.0e-3);
            const double volume = std::pow(2.0 * Physics::pi, 1.5) * spec.sigma[0] * spec.sigma[1]
                                  * spec.sigma[2];
            const double spacing = std::cbrt(volume / spec.globalCount);
            struct Case {
                double factor;  // h in mean spacings; 0 = leaf-based
                BarnesHutPrecision precision;
            };
            const auto D = BarnesHutPrecision::Double, M = BarnesHutPrecision::Mixed;
            for (const Case& c : {Case{0.25, D}, Case{0.5, D}, Case{1.0, D}, Case{2.0, D},
                                  Case{4.0, D}, Case{0.0, D}, Case{0.5, M}, Case{1.0, M}}) {
                auto values                = config(c.factor > 0.0 ? c.factor * spacing : 0.0);
                values.leafBasedSoftening  = c.factor == 0.0;
                values.precision           = c.precision;
                values.traversalStatistics = true;
                Run run(values, spec);
                (void)run.solve();  // builds the domain
                constexpr int repeats = 5;
                const auto start      = std::chrono::steady_clock::now();
                for (int r = 0; r < repeats; ++r) {
                    (void)run.solve();
                }
                const double seconds =
                        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
                                .count()
                        / repeats;
                if (ippl::Comm->rank() == 0) {
                    std::printf("[BH] solve time %s, h = %s: %.2f ms per solve\n",
                                c.precision == M ? "MIXED " : "DOUBLE",
                                c.factor > 0.0 ? (std::to_string(c.factor) + " spacings").c_str()
                                               : "leaf-based",
                                1e3 * seconds);
                }
            }
        }

        TEST_F(BarnesHutAlgorithmTest, DISABLED_RestFrameAgainstGaussianContinuum) {
            // Diagnostic: BH (several h) and binned OPEN (several meshes) against the boosted
            // continuum field. Measured, needle 1 x 1 x 10 mm (spacing 9.2e-5 m): BH h = 1e-5 ..
            // 2e-4 -> 10.6, 7.5, 4.6, 3.2, 2.8 % (~h^-1/2); OPEN 32^3, 64^3, 128^3 -> 3.3, 1.9,
            // 2.8 %. 1 x 1 x 2 mm: BH 5.7, 4.0, 2.4, 1.7, 2.5 %; OPEN 3.4, 1.7, 2.3 %. B as E.
            for (double sigmaZ : {2.0e-3, 0.4e-3}) {
                BunchSpec spec;
                spec.globalCount          = 200000;
                spec.sigma                = Vector(1.0e-3, 1.0e-3, sigmaZ);
                spec.transverseMomentum   = 0.05;
                spec.longitudinalMomentum = std::sqrt(24.0);
                spec.momentumSpread       = 0.02;
                Run sample(restFrameConfig(), spec);
                const auto particles = gatherAll(sample.particles);
                const auto continuum =
                        gaussianContinuum(particles, spec, kCharge * spec.globalCount);
                const double volume = std::pow(2.0 * Physics::pi, 1.5) * spec.sigma[0]
                                      * spec.sigma[1] * continuum.gamma * spec.sigma[2];
                const double spacing = std::cbrt(volume / spec.globalCount);
                if (ippl::Comm->rank() == 0) {
                    std::printf("[BH] continuum: lab sigma_z=%.1e m, gamma=%.3f, rest-frame mean "
                                "spacing %.2e m\n",
                                sigmaZ, continuum.gamma, spacing);
                }
                const auto report = [&](const char* label, const std::vector<Record>& records) {
                    ASSERT_EQ(records.size(), particles.size());
                    const double errorE =
                            rmsRelativeDifference(fieldsOf(records, &Record::e), continuum.e);
                    const double errorB =
                            rmsRelativeDifference(fieldsOf(records, &Record::b), continuum.b);
                    if (ippl::Comm->rank() == 0) {
                        std::printf("[BH] continuum: %-22s rms relative error E=%.3e B=%.3e\n",
                                    label, errorE, errorB);
                    }
                };
                for (double softening : {1.0e-5, 2.0e-5, 5.0e-5, 1.0e-4, 2.0e-4}) {
                    Run bh(restFrameConfig(softening), spec);
                    (void)bh.solve();
                    char label[64];
                    std::snprintf(label, sizeof(label), "BH h=%.0e (%.2f sp)", softening,
                                  softening / spacing);
                    report(label, gatherAll(bh.particles));
                }
                for (std::size_t mesh : {32u, 64u, 128u}) {
                    char label[64];
                    std::snprintf(label, sizeof(label), "OPEN %zu^3", mesh);
                    report(label, binnedOpenFields(spec, mesh, DirichletPlaneConfig{}));
                }
            }
        }

    }  // namespace
}  // namespace opalx::spacecharge
