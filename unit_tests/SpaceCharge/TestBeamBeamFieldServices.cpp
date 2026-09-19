#include <gtest/gtest.h>

#include "AbstractObjects/OpalData.h"
#include "PartBunch/BunchStateHandler.h"
#include "Physics/Physics.h"
#include "SpaceCharge/CartesianPIC3D/CartesianPIC3DAlgorithm.h"
#include "SpaceCharge/CartesianPIC3D/FieldMirror.hpp"
#include "SpaceCharge/SpaceChargeSolver.h"
#include "Structure/DataSink.h"
#include "Utilities/Options.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

namespace opalx::spacecharge {
    namespace {
        template <typename View>
        auto snapshot(const View& view) {
            auto host = Kokkos::create_mirror(view);
            Kokkos::deep_copy(host, view);
            return host;
        }

        class BeamBeamFieldServicesTest : public ::testing::Test {
        protected:
            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
                OpalData::getInstance()->storeInputFn("beambeam_field_services.opal");
                gmsg                = new Inform(nullptr, -1);
                Options::enableHDF5 = false;
            }
            static void TearDownTestSuite() {
                delete gmsg;
                gmsg = nullptr;
                std::remove("beambeam_field_services.stat");
                std::remove("beambeam_field_services.lbal");
                ippl::finalize();
            }
        };

        TEST_F(BeamBeamFieldServicesTest, CopyDoublesEAndCancelsBAtExactOverlap) {
            CartesianDomainConfig3D domainConfig;
            domainConfig.meshSize      = {8, 8, 8};
            domainConfig.decomposition = {true, true, true};
            CartesianDomain<double, 3> domain(domainConfig);
            auto fields     = std::make_unique<CartesianPIC3DFieldStorage<double, 3>>(domain);
            auto state      = std::make_shared<BunchStateHandler>();
            using Container = ::ParticleContainer<double, 3>;
            Container primary(domain.mesh(), domain.layout());
            primary.setBunchStateHandler(state);
            // A finite field scale keeps the absolute roundoff floor well below physical B.
            primary.setQ(1.0e-9);
            primary.setM(Physics::m_e);
            constexpr std::size_t count = 64;
            primary.createParticles(count);
            auto positions = primary.R.getHostMirror();
            auto momenta   = primary.P.getHostMirror();
            auto timeSteps = primary.dt.getHostMirror();
            for (std::size_t pair = 0; pair < count / 2; ++pair) {
                const double phase = static_cast<double>(pair + 1);
                const double x     = 0.36 * std::sin(0.7 * phase);
                const double y     = 0.32 * std::cos(1.1 * phase);
                const double z     = 0.16 * (0.25 + static_cast<double>(pair % 7) / 7.0);
                for (std::size_t side = 0; side < 2; ++side) {
                    const auto index = 2 * pair + side;
                    positions(index) = Vector_t<double, 3>(x, y, side == 0 ? -z : z);
                    momenta(index)   = Vector_t<double, 3>(0.0, 0.0, 2.0);
                    timeSteps(index) = 1.0e-12;
                }
            }
            Kokkos::deep_copy(primary.R.getView(), positions);
            Kokkos::deep_copy(primary.P.getView(), momenta);
            Kokkos::deep_copy(primary.dt.getView(), timeSteps);
            primary.markMomentsDirty();
            primary.updateMoments();

            CartesianPIC3DConfig config;
            config.backend                      = PoissonSolverType::Open;
            config.grid.meshSize                = domainConfig.meshSize;
            config.grid.decomposition           = domainConfig.decomposition;
            config.binning                      = BinningConfig{};
            config.binning->maximumBins         = 1;
            config.binning->adaptive            = false;
            config.binning->tablePrintFrequency = 0;
            // BeamBeam must freeze ORB even when the run requests redistribution every step.
            config.repartitionFrequency   = 1;
            config.loadBalancingThreshold = 0.0;
            std::vector<Container*> containers{&primary};
            DataSink sink;
            SpaceChargeSolver solver(
                    std::make_unique<CartesianPIC3DAlgorithm>(
                            config, containers, std::move(fields), &sink, state),
                    1);
            auto& services = solver.beamBeamFields();
            state->setFixedCartesianDomain({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0});
            SpaceChargeStepState step;
            step.timeStep = 1.0e-12;
            step.mpiSize  = ippl::Comm->size();
            const std::array<std::uint8_t, 1> activity{1};
            SpaceChargeSolveContext context(activity, step);
            // Ordinary fixed-domain solves retain the existing configuration validation.
            EXPECT_THROW(solver.solve(context), OpalException);
            services.configure(BeamBeamSolvePolicy{true, false});
            solver.solve(context);
            EXPECT_EQ(solver.redistributionCount(), 0u);
            ASSERT_TRUE(services.depositedCharge().has_value());
            const double physicalCharge = *services.depositedCharge();
            EXPECT_NEAR(
                    physicalCharge, primary.getTotalCharge(),
                    1.0e-12 * std::abs(primary.getTotalCharge()));
            const auto physicalE        = snapshot(primary.E.getView());
            const auto physicalB        = snapshot(primary.B.getView());
            const auto originalR        = snapshot(primary.R.getView());
            const auto solvesBeforeCopy = solver.backendSolveCount();
            services.configure(BeamBeamSolvePolicy{true, true});
            solver.solve(context);
            EXPECT_EQ(solver.backendSolveCount(), solvesBeforeCopy + 1);
            EXPECT_EQ(solver.redistributionCount(), 0u);
            ASSERT_TRUE(services.depositedCharge().has_value());
            EXPECT_NEAR(
                    *services.depositedCharge(), 2.0 * physicalCharge,
                    1.0e-12 * std::abs(physicalCharge));
            const auto copiedE   = snapshot(primary.E.getView());
            const auto copiedB   = snapshot(primary.B.getView());
            const auto restoredR = snapshot(primary.R.getView());
            double electricScale = 0.0;
            double magneticScale = 0.0;
            for (std::size_t index = 0; index < primary.getLocalNum(); ++index) {
                for (unsigned dimension = 0; dimension < 3; ++dimension) {
                    electricScale = std::max(electricScale, std::abs(physicalE(index)[dimension]));
                    magneticScale = std::max(magneticScale, std::abs(physicalB(index)[dimension]));
                }
            }
            if (primary.getLocalNum() > 0) {
                EXPECT_GT(electricScale, 0.0);
                EXPECT_GT(magneticScale, 1.0e-10);
            }
            const double electricTolerance = std::max(1.0e-12, 1.0e-10 * electricScale);
            const double magneticTolerance = std::max(1.0e-12, 1.0e-10 * magneticScale);
            for (std::size_t index = 0; index < primary.getLocalNum(); ++index) {
                for (unsigned dimension = 0; dimension < 3; ++dimension) {
                    EXPECT_NEAR(
                            copiedE(index)[dimension], 2.0 * physicalE(index)[dimension],
                            electricTolerance);
                    EXPECT_NEAR(copiedB(index)[dimension], 0.0, magneticTolerance);
                    EXPECT_DOUBLE_EQ(restoredR(index)[dimension], originalR(index)[dimension]);
                }
            }

            // Copy-only fields reverse B while retaining E at exact overlap.
            services.configure(BeamBeamSolvePolicy{false, true});
            solver.solve(context);
            const auto onlyCopyE = snapshot(primary.E.getView());
            const auto onlyCopyB = snapshot(primary.B.getView());
            for (std::size_t index = 0; index < primary.getLocalNum(); ++index) {
                for (unsigned dimension = 0; dimension < 3; ++dimension) {
                    EXPECT_NEAR(
                            onlyCopyE(index)[dimension], physicalE(index)[dimension],
                            electricTolerance);
                    EXPECT_NEAR(
                            onlyCopyB(index)[dimension], -physicalB(index)[dimension],
                            magneticTolerance);
                }
            }

            services.configure(BeamBeamSolvePolicy{false, false});
            const auto solvesBeforeEmpty = solver.backendSolveCount();
            solver.solve(context);
            EXPECT_EQ(solver.backendSolveCount(), solvesBeforeEmpty);
            EXPECT_EQ(services.depositedCharge(), std::optional<double>(0.0));
            services.gatherFields(primary);
            const auto emptyE = snapshot(primary.E.getView());
            for (std::size_t index = 0; index < primary.getLocalNum(); ++index) {
                for (unsigned dimension = 0; dimension < 3; ++dimension) {
                    EXPECT_DOUBLE_EQ(emptyE(index)[dimension], 0.0);
                }
            }
            EXPECT_EQ(solver.redistributionCount(), 0u);
            EXPECT_EQ(services.configuration().repartitionFrequency, 1u);
            services.configure(std::nullopt);
            EXPECT_THROW(solver.solve(context), OpalException);
        }

        TEST_F(BeamBeamFieldServicesTest, MirrorFieldZHandlesNonSlab3DDecomposition) {
            if (ippl::Comm->size() != 4) {
                GTEST_SKIP() << "This mirror-field decomposition check is defined for 4 MPI ranks.";
            }
            constexpr unsigned Dim = 3;
            ippl::NDIndex<Dim> domain;
            for (unsigned axis = 0; axis < Dim; ++axis) {
                domain[axis] = ippl::Index(8);
            }
            std::array<bool, Dim> decomposition{true, true, true};
            Vector_t<double, Dim> spacing(1.0), origin(0.0);
            FieldLayout_t<Dim> layout(MPI_COMM_WORLD, domain, decomposition, false);
            Mesh_t<Dim> mesh(domain, spacing, origin);
            std::vector<ippl::NDIndex<Dim>> domains(4);
            domains[0][0] = ippl::Index(0, 3);
            domains[0][1] = ippl::Index(0, 7);
            domains[0][2] = ippl::Index(0, 3);
            domains[1][0] = ippl::Index(4, 7);
            domains[1][1] = ippl::Index(0, 7);
            domains[1][2] = ippl::Index(0, 3);
            domains[2][0] = ippl::Index(0, 7);
            domains[2][1] = ippl::Index(0, 3);
            domains[2][2] = ippl::Index(4, 7);
            domains[3][0] = ippl::Index(0, 7);
            domains[3][1] = ippl::Index(4, 7);
            domains[3][2] = ippl::Index(4, 7);
            layout.updateLayout(domains);
            VField_t<double, Dim> source, destination;
            source.initialize(mesh, layout);
            destination.initialize(mesh, layout);
            source            = 0.0;
            const auto& local = layout.getLocalNDIndex();
            const int ghost   = source.getNghost();
            auto values       = source.getHostMirror();
            Kokkos::deep_copy(values, source.getView());
            for (int i = local[0].first(); i <= local[0].last(); ++i) {
                for (int j = local[1].first(); j <= local[1].last(); ++j) {
                    for (int k = local[2].first(); k <= local[2].last(); ++k) {
                        const double base = i + 10 * j + 100 * k;
                        values(i - local[0].first() + ghost, j - local[1].first() + ghost,
                               k - local[2].first() + ghost) =
                                Vector_t<double, 3>(base, 1000.0 + base, 2000.0 + base);
                    }
                }
            }
            Kokkos::deep_copy(source.getView(), values);
            opalx::detail::mirrorField(source, destination, 2);
            const auto mirrored =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), destination.getView());
            for (int i = local[0].first(); i <= local[0].last(); ++i) {
                for (int j = local[1].first(); j <= local[1].last(); ++j) {
                    for (int k = local[2].first(); k <= local[2].last(); ++k) {
                        const double base = i + 10 * j + 100 * (7 - k);
                        for (unsigned dimension = 0; dimension < 3; ++dimension) {
                            EXPECT_DOUBLE_EQ(
                                    mirrored(
                                            i - local[0].first() + ghost,
                                            j - local[1].first() + ghost,
                                            k - local[2].first() + ghost)[dimension],
                                    1000.0 * dimension + base);
                        }
                    }
                }
            }
        }
    }  // namespace
}  // namespace opalx::spacecharge
