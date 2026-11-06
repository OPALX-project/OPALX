/** Distributed FFT2D5 against direct discrete free-space convolution (no FFT oracle).
 * Run with 1, 2 and 4 MPI ranks. Sources initially reside only on rank zero;
 * they straddle every slab boundary and the periodic seam. All three fields,
 * Cartesian positions, momenta, IDs, charge and repeated migration are checked.
 */
#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include "PartBunch/BunchStateHandler.h"
#include "Physics/Physics.h"
#include "SpaceCharge/FFT2D5/FFT2D5Algorithm.h"
#include "Utilities/Options.h"

namespace opalx::spacecharge {
    namespace {
        using V          = Vector_t<double, 3>;
        constexpr int nx = 8, nz = 16;
        constexpr double dx = 1. / nx, dz = 1. / nz, q = 1.e-12, dt = 1.e-12;
        using Grid = std::array<double, nx * nx * nz>;
        int index(int i, int j, int k) { return (k * nx + j) * nx + i; }
        class FFT2D5ParallelTest : public ::testing::Test {
        protected:
            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
            }
            static void TearDownTestSuite() { ippl::finalize(); }
        };

        TEST_F(FFT2D5ParallelTest, MigrationHalosAndFieldsMatchDirectConvolution) {
            const auto path = std::filesystem::current_path() / "fft2d5-mpi-path.dat";
            if (ippl::Comm->rank() == 0) {
                std::ofstream out(path);
                out << "s rx ry rz px py pz ex ey ez bx by bz ke t\n"
                    << "0 0 0 0 0 0 1 0 0 0 0 0 0 0 0\n"
                    << "1 0 0 1 0 0 1 0 0 0 0 0 0 0 0\n";
            }
            MPI_Barrier(MPI_COMM_WORLD);
            const int previousBatch  = Options::fft2d5BatchSize;
            Options::fft2d5BatchSize = 8;
            for (bool closed : {false, true})
                for (bool cicZ : {false, true})
                    for (bool sparse : {false, true}) {
                        std::vector<V> sources =
                                sparse ? std::vector<V>{V(-.0625, .0625, .5001)}
                                       : std::vector<V>{
                                                 V(-.0625, .0625, .0001),  V(.0625, -.0625, .2499),
                                                 V(.1875, .0625, .2501),   V(-.1875, -.0625, .4999),
                                                 V(-.0625, .1875, .5001),  V(.0625, -.1875, .5002),
                                                 V(-.0625, -.0625, .7499), V(.0625, .0625, .7501),
                                                 V(.1875, -.0625, .9999),  V(-.1875, .1875, 0.),
                                                 V(.1875, -.1875, 1.)};
                        ippl::NDIndex<3> domain{nx, nx, nz};
                        Mesh_t<3> mesh(domain, V(dx, dx, dz), V(-.5, -.5, 0.));
                        FieldLayout_t<3> layout(
                                MPI_COMM_WORLD, domain, {false, false, true}, closed);
                        auto state = std::make_shared<BunchStateHandler>();
                        ::ParticleContainer<double> pc(mesh, layout);
                        pc.setBunchStateHandler(state);
                        pc.setQ(q);
                        pc.setM(Physics::m_e);
                        pc.createParticles(ippl::Comm->rank() == 0 ? sources.size() : 0);
                        if (ippl::Comm->rank() == 0) {
                            auto r  = pc.R.getHostMirror();
                            auto p  = pc.P.getHostMirror();
                            auto id = pc.ID.getHostMirror();
                            for (size_t n = 0; n < sources.size(); ++n) {
                                r(n)  = sources[n];
                                p(n)  = V(0., 0., .4);
                                id(n) = n;
                            }
                            Kokkos::deep_copy(pc.R.getView(), r);
                            Kokkos::deep_copy(pc.P.getView(), p);
                            Kokkos::deep_copy(pc.ID.getView(), id);
                        }
                        Kokkos::deep_copy(pc.dt.getView(), dt);
                        Kokkos::deep_copy(pc.InvalidMask.getView(), false);
                        pc.markMomentsDirty();
                        pc.updateMoments();
                        FFT2D5Config cfg;
                        cfg.grid.meshSize      = {nx, nx, nz};
                        cfg.grid.decomposition = {false, false, true};
                        cfg.pipeSizeX = cfg.pipeSizeY = 1.;
                        cfg.beamRadius                = .1;
                        cfg.closedRing                = closed;
                        cfg.scatterLongitudinally     = cicZ;
                        cfg.longitudinalFieldMode     = FFT2D5LongitudinalFieldMode::Open;
                        cfg.referencePathFile         = path.string();
                        std::array containers{&pc};
                        FFT2D5Algorithm solver(cfg, containers, state);
                        SpaceChargeStepState step;
                        step.timeStep = dt;
                        step.mpiSize  = ippl::Comm->size();
                        const std::array<std::uint8_t, 1> activity{1};
                        Grid rho{}, ex{}, ey{};
                        for (const auto& r : sources) {
                            const int i       = int(std::lround((r[0] + .5) / dx - .5)),
                                      j       = int(std::lround((r[1] + .5) / dx - .5));
                            const double cell = r[2] / dz - .5;
                            const int left    = int(std::floor(cell));
                            const double f    = cell - left;
                            for (int side = 0; side < (cicZ ? 2 : 1); ++side) {
                                int k = cicZ ? left + side : int(std::floor(r[2] / dz + .5));
                                const double w = cicZ ? (side ? f : 1. - f) : 1.;
                                if (closed) k = (k % nz + nz) % nz;
                                if (k >= 0 && k < nz) rho[index(i, j, k)] += q * w / (dx * dx * dz);
                            }
                        }
                        auto phi = [&](int i, int j, int k) {
                            double sum = 0.;
                            for (int a = 0; a < nx; ++a)
                                for (int b = 0; b < nx; ++b)
                                    sum += rho[index(a, b, k)] * dx * dx / Physics::epsilon_0
                                           * FFT2D5Poisson::cellAveragedGreensFunction(
                                                   (i - a) * dx, (j - b) * dx, dx, dx);
                            return sum;
                        };
                        for (int k = 0; k < nz; ++k)
                            for (int i = 0; i < nx; ++i)
                                for (int j = 0; j < nx; ++j) {
                                    ex[index(i, j, k)] =
                                            -(phi(i + 1, j, k) - phi(i - 1, j, k)) / (2 * dx);
                                    ey[index(i, j, k)] =
                                            -(phi(i, j + 1, k) - phi(i, j - 1, k)) / (2 * dx);
                                }
                        std::array<double, nz> lambda{}, gradient{};
                        for (int k = 0; k < nz; ++k)
                            for (int i = 0; i < nx; ++i)
                                for (int j = 0; j < nx; ++j)
                                    lambda[k] += rho[index(i, j, k)] * dx * dx;
                        for (int k = 0; k < nz; ++k)
                            gradient[k] =
                                    ((k + 1 < nz ? lambda[k + 1] : (closed ? lambda[0] : 0.))
                                     - (k > 0 ? lambda[k - 1] : (closed ? lambda[nz - 1] : 0.)))
                                    / (2 * dz);
                        for (int repeat = 0; repeat < 2; ++repeat) {
                            SCOPED_TRACE(
                                    ::testing::Message()
                                    << "closed=" << closed << " cic=" << cicZ
                                    << " sparse=" << sparse << " repeat=" << repeat
                                    << " rank=" << ippl::Comm->rank());
                            ASSERT_NO_THROW(static_cast<void>(solver.solve({activity, step})));
                            const auto local   = solver.getRho()->getLayout().getLocalNDIndex();
                            const auto density = Kokkos::create_mirror_view_and_copy(
                                    Kokkos::HostSpace(), solver.getRho()->getView());
                            const auto field = Kokkos::create_mirror_view_and_copy(
                                    Kokkos::HostSpace(), solver.getEField()->getView());
                            for (int k = 0; k < local[2].length(); ++k)
                                for (int i = 0; i < nx; ++i)
                                    for (int j = 0; j < nx; ++j) {
                                        const int n = index(i, j, k + local[2].first());
                                        EXPECT_NEAR(density(i + 1, j + 1, k + 1), rho[n], 1.e-21);
                                        EXPECT_NEAR(field(i + 1, j + 1, k + 1)[0], ex[n], 1.e-11);
                                        EXPECT_NEAR(field(i + 1, j + 1, k + 1)[1], ey[n], 1.e-11);
                                    }
                            const auto g = Kokkos::create_mirror_view_and_copy(
                                    Kokkos::HostSpace(), solver.getLineDensityGradient());
                            for (int k = 0; k < nz; ++k)
                                EXPECT_NEAR(g(k), gradient[k], 1.e-22);
                            const auto positions = Kokkos::create_mirror_view_and_copy(
                                    Kokkos::HostSpace(), pc.R.getView());
                            const auto momenta = Kokkos::create_mirror_view_and_copy(
                                    Kokkos::HostSpace(), pc.P.getView());
                            const auto ids = Kokkos::create_mirror_view_and_copy(
                                    Kokkos::HostSpace(), pc.ID.getView());
                            const auto electric = Kokkos::create_mirror_view_and_copy(
                                    Kokkos::HostSpace(), pc.E.getView());
                            const auto magnetic = Kokkos::create_mirror_view_and_copy(
                                    Kokkos::HostSpace(), pc.B.getView());
                            size_t count = pc.getLocalNum(), total = 0;
                            ippl::Comm->allreduce(count, total, 1, std::plus<size_t>());
                            EXPECT_EQ(total, sources.size());
                            for (size_t n = 0; n < count; ++n) {
                                ASSERT_GE(ids(n), 0);
                                ASSERT_LT(size_t(ids(n)), sources.size());
                                const auto r          = sources[ids(n)];
                                const int i           = int(std::lround((r[0] + .5) / dx - .5));
                                const int j           = int(std::lround((r[1] + .5) / dx - .5));
                                const double cell     = r[2] / dz - .5;
                                const int left        = int(std::floor(cell));
                                const double fraction = cell - left;
                                auto sample           = [&](const Grid& values, int k) {
                                    if (closed) k = (k % nz + nz) % nz;
                                    return k >= 0 && k < nz ? values[index(i, j, k)] : 0.;
                                };
                                double expectedX, expectedY;
                                if (cicZ) {
                                    expectedX = (1. - fraction) * sample(ex, left)
                                                + fraction * sample(ex, left + 1);
                                    expectedY = (1. - fraction) * sample(ey, left)
                                                + fraction * sample(ey, left + 1);
                                } else {
                                    expectedX = sample(ex, int(std::floor(r[2] / dz + .5)));
                                    expectedY = sample(ey, int(std::floor(r[2] / dz + .5)));
                                }
                                auto slope = [&](int k) {
                                    k = closed ? (k % nz + nz) % nz : std::clamp(k, 0, nz - 1);
                                    return gradient[k];
                                };
                                const double derivative =
                                        cicZ ? (1. - fraction) * slope(left)
                                                        + fraction * slope(left + 1)
                                             : slope(int(std::floor(r[2] / dz)));
                                const double expectedZ = -6.36
                                                         / (4 * Physics::pi * Physics::epsilon_0)
                                                         * derivative / 1.16;
                                EXPECT_NEAR(electric(n)[0], expectedX, 1.e-11);
                                EXPECT_NEAR(electric(n)[1], expectedY, 1.e-11);
                                EXPECT_NEAR(electric(n)[2], expectedZ, 1.e-11);
                                const double beta = .4 / std::sqrt(1.16);
                                EXPECT_NEAR(magnetic(n)[0], -beta * expectedY / Physics::c, 1.e-19);
                                EXPECT_NEAR(magnetic(n)[1], beta * expectedX / Physics::c, 1.e-19);
                                EXPECT_DOUBLE_EQ(magnetic(n)[2], 0.);
                                for (int d = 0; d < 3; ++d) {
                                    EXPECT_DOUBLE_EQ(positions(n)[d], sources[ids(n)][d]);
                                    EXPECT_DOUBLE_EQ(momenta(n)[d], d == 2 ? .4 : 0.);
                                }
                            }
                        }
                    }
            Options::fft2d5BatchSize = previousBatch;
            MPI_Barrier(MPI_COMM_WORLD);
            if (ippl::Comm->rank() == 0) std::filesystem::remove(path);
        }
    }  // namespace
}  // namespace opalx::spacecharge
