/**
 * @file TestCartesianDomainLoadBalancing.cpp
 * @brief Four-rank tests of pair-weighted ORB load balancing for P3M
 * (`OPTION, LOADBALANCINGPAIRWEIGHT`).
 *
 * The P3M short-range work grows with the number of particle pairs inside the cutoff, not with
 * the particle count, so a dense cluster can overload one rank even when every rank holds the
 * same number of particles. The test bunch is a uniform lattice plus a dense cluster with 11 % of
 * the particles but most of the pairs. Covers:
 * - Weighting ORB by pair counts spreads the pairs more evenly over the ranks than a cut by
 *   particle count, at the cost of less even particle counts; no particle is lost.
 * - If the pair-weighted cut would leave a rank region narrower than twice the P3M cutoff, the
 *   updater falls back to the cut by particle count and keeps the full cutoff.
 *
 * The tests need exactly four ranks (two ORB levels) and are skipped otherwise, so CTest only
 * runs them on one rank. Run them by hand from the build directory:
 * `mpirun -n 4 ./unit_tests/SpaceCharge/TestCartesianDomainLoadBalancing`
 */

#include <gtest/gtest.h>

#include "AbstractObjects/OpalData.h"
#include "PartBunch/BunchStateHandler.h"
#include "PartBunch/CartesianDomain.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"
#include "SpaceCharge/CartesianPIC3D/CartesianDomainUpdater.h"
#include "SpaceCharge/CartesianPIC3D/CartesianPIC3DFieldStorage.h"
#include "SpaceCharge/Poisson/PoissonSolver.h"
#include "SpaceCharge/SpaceChargeSolveContext.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

#include <array>
#include <cstdio>
#include <functional>
#include <memory>
#include <vector>

extern Inform* gmsg;

namespace opalx::spacecharge {
    namespace {

        using Container = ::ParticleContainer<double, 3>;

        // A uniform 16^3 lattice over [-0.94, 0.94]^3 plus a dense 8^3 cluster, about three mesh
        // cells wide, in the +x+y+z octant. The cluster holds 11 % of the particles but most of
        // the pairs.
        constexpr std::size_t latticeSide = 16;
        constexpr std::size_t clusterSide = 8;
        constexpr std::size_t totalParticles =
                latticeSide * latticeSide * latticeSide + clusterSide * clusterSide * clusterSide;
        const ippl::Vector<double, 3> clusterCenter(0.5, 0.5, 0.5);

        class CartesianDomainLoadBalancingTest : public ::testing::Test {
        protected:
            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
                OpalData::getInstance()->storeInputFn("cartesian_domain_load_balancing.opal");
                gmsg                = new Inform(nullptr, -1);
                Options::enableHDF5 = false;
            }

            static void TearDownTestSuite() {
                delete gmsg;
                gmsg = nullptr;
                std::remove("cartesian_domain_load_balancing.stat");
                std::remove("cartesian_domain_load_balancing.lbal");
                ippl::finalize();
            }

            void SetUp() override {
                if (ippl::Comm->size() != 4) {
                    GTEST_SKIP() << "Needs 4 MPI ranks for a two-level ORB.";
                }
            }
        };

        /** @brief One P3M domain with the lattice-plus-cluster bunch created on rank 0. */
        class Scenario {
        public:
            Scenario(double cutoffCells, double pairWeight)
                : domain_m(domainConfig(cutoffCells)), workspace_m(domain_m) {
                workspace_m.initializeFields(PoissonSolverType::None);
                particles_m = std::make_unique<Container>(
                        domain_m.mesh(), domain_m.layout(), false,
                        Container::LayoutType::SpatialOverlap,
                        resolveOverlapCutoff<3>(0.0, cutoffCells, domain_m.spacing()),
                        ippl::BC::NO);
                particles_m->setBunchStateHandler(bunchState_m);
                createBunch();
                containers_m = {particles_m.get()};

                PoissonSolverConfig poissonConfig;
                poissonConfig.type = PoissonSolverType::None;
                poisson_m          = makePoissonSolver(
                        poissonConfig,
                        {&workspace_m.chargeDensity(), &workspace_m.electricField()});

                CartesianPIC3DConfig values;
                values.backend                 = PoissonSolverType::P3M;
                values.grid.meshSize           = {latticeSide, latticeSide, latticeSide};
                values.grid.decomposition      = {true, true, true};
                values.p3mCutoffCells          = cutoffCells;
                values.repartitionFrequency    = 1;
                values.loadBalancingThreshold  = 0.0;
                values.loadBalancingPairWeight = pairWeight;
                updater_m = std::make_unique<CartesianDomainUpdater>(values, containers_m);
            }

            bool updateBeamFrame() {
                SpaceChargeStepState step;
                step.mpiSize = ippl::Comm->size();
                const std::array<std::uint8_t, 1> activity{1};
                SpaceChargeSolveContext context(activity, step);
                return updater_m->updateForSolve(
                        DomainCoordinateFrame::Beam, context, {}, nullptr, workspace_m, *poisson_m);
            }

            [[nodiscard]] Container& particles() { return *particles_m; }
            [[nodiscard]] CartesianDomain<double, 3>& domain() { return domain_m; }

            [[nodiscard]] double requestedCutoff(double cutoffCells) const {
                return resolveOverlapCutoff<3>(0.0, cutoffCells, domain_m.mesh().getMeshSpacing());
            }

            /**
             * @brief Largest over mean rank count of local particle pairs closer than radius.
             *
             * Counts only pairs whose particles share a rank, which is enough to compare two
             * decompositions of the same bunch. Collective.
             */
            [[nodiscard]] double pairMaxOverMean(double radius) const {
                const auto position = particles_m->R.getHostMirror();
                Kokkos::deep_copy(position, particles_m->R.getView());
                const std::size_t n = particles_m->getLocalNum();
                double pairs        = 0.0;
                for (std::size_t i = 0; i < n; ++i) {
                    for (std::size_t j = i + 1; j < n; ++j) {
                        const Vector_t<double, 3> distance = position(i) - position(j);
                        if (distance.dot(distance) < radius * radius) {
                            pairs += 1.0;
                        }
                    }
                }
                double extremes[2] = {pairs, pairs};
                ippl::Comm->allreduce(&extremes[0], 1, std::greater<double>());
                ippl::Comm->allreduce(&extremes[1], 1, std::plus<double>());
                return extremes[0] * ippl::Comm->size() / extremes[1];
            }

            /** @brief Largest local particle count over the mean count. Collective. */
            [[nodiscard]] double countMaxOverMean() const {
                double largest = static_cast<double>(particles_m->getLocalNum());
                ippl::Comm->allreduce(&largest, 1, std::greater<double>());
                return largest * ippl::Comm->size() / static_cast<double>(totalParticles);
            }

        private:
            static CartesianDomainConfig3D domainConfig(double cutoffCells) {
                CartesianDomainConfig3D config;
                config.meshSize           = {latticeSide, latticeSide, latticeSide};
                config.decomposition      = {true, true, true};
                config.layoutType         = ParticleLayoutType::SpatialOverlap;
                config.overlapCutoffCells = cutoffCells;
                return config;
            }

            void createBunch() {
                const std::size_t local = ippl::Comm->rank() == 0 ? totalParticles : 0;
                particles_m->createParticles(local);
                particles_m->setM(Physics::m_e);
                Kokkos::deep_copy(particles_m->P.getView(), Vector_t<double, 3>(0.0));
                if (local == 0) {
                    return;
                }
                auto positions    = particles_m->R.getHostMirror();
                std::size_t index = 0;
                const double step = 1.875 / latticeSide;
                for (std::size_t i = 0; i < latticeSide; ++i) {
                    for (std::size_t j = 0; j < latticeSide; ++j) {
                        for (std::size_t k = 0; k < latticeSide; ++k) {
                            positions(index++) = Vector_t<double, 3>(
                                    -0.9375 + (i + 0.5) * step, -0.9375 + (j + 0.5) * step,
                                    -0.9375 + (k + 0.5) * step);
                        }
                    }
                }
                const double clusterStep = 0.05;
                const double clusterLow  = clusterCenter[0] - 0.5 * clusterStep * clusterSide;
                for (std::size_t i = 0; i < clusterSide; ++i) {
                    for (std::size_t j = 0; j < clusterSide; ++j) {
                        for (std::size_t k = 0; k < clusterSide; ++k) {
                            positions(index++) = Vector_t<double, 3>(
                                    clusterLow + (i + 0.5) * clusterStep,
                                    clusterLow + (j + 0.5) * clusterStep,
                                    clusterLow + (k + 0.5) * clusterStep);
                        }
                    }
                }
                Kokkos::deep_copy(particles_m->R.getView(), positions);
            }

            CartesianDomain<double, 3> domain_m;
            CartesianPIC3DFieldStorage<double, 3> workspace_m;
            std::shared_ptr<BunchStateHandler> bunchState_m = std::make_shared<BunchStateHandler>();
            std::unique_ptr<Container> particles_m;
            std::vector<Container*> containers_m;
            std::unique_ptr<PoissonSolver> poisson_m;
            std::unique_ptr<CartesianDomainUpdater> updater_m;
        };

        std::size_t globalCount(const Container& particles) {
            std::size_t count = particles.getLocalNum();
            ippl::Comm->allreduce(&count, 1, std::plus<std::size_t>());
            return count;
        }

        TEST_F(CartesianDomainLoadBalancingTest, PairWeightEvensOutPairCounts) {
            Scenario byCount(0.5, 0.0);
            Scenario byPairs(0.5, 1.0);
            ASSERT_TRUE(byCount.updateBeamFrame());
            ASSERT_TRUE(byPairs.updateBeamFrame());
            EXPECT_EQ(globalCount(byCount.particles()), totalParticles);
            EXPECT_EQ(globalCount(byPairs.particles()), totalParticles);

            // Balancing pairs splits the cluster over more ranks and moves lattice particles onto
            // the others, so pair counts even out while particle counts spread apart.
            const double radius = 0.1;
            EXPECT_LT(byPairs.pairMaxOverMean(radius), byCount.pairMaxOverMean(radius));
            EXPECT_GT(byPairs.countMaxOverMean(), byCount.countMaxOverMean());
            EXPECT_DOUBLE_EQ(
                    byPairs.particles().getP3MLayout().getCutoff(), byPairs.requestedCutoff(0.5));
        }

        TEST_F(CartesianDomainLoadBalancingTest, PairWeightFallsBackWhenTheCutoffDoesNotFit) {
            Scenario byCount(2.5, 0.0);
            Scenario byPairs(2.5, 1.0);
            ASSERT_TRUE(byCount.updateBeamFrame());
            ASSERT_TRUE(byPairs.updateBeamFrame());

            // The pair-weighted cut leaves a 4-cell region next to the cluster, narrower than two
            // cutoffs of 2.5 cells, while the particle-count cut keeps every region at 6 cells or
            // more. The updater restores the layout and cuts by particle count, keeping the full
            // cutoff.
            EXPECT_DOUBLE_EQ(
                    byPairs.particles().getP3MLayout().getCutoff(), byPairs.requestedCutoff(2.5));
            const auto& countDomain = byCount.domain().layout().getLocalNDIndex();
            const auto& pairDomain  = byPairs.domain().layout().getLocalNDIndex();
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_EQ(pairDomain[d].first(), countDomain[d].first());
                EXPECT_EQ(pairDomain[d].last(), countDomain[d].last());
            }
            EXPECT_EQ(globalCount(byPairs.particles()), totalParticles);
        }

    }  // namespace
}  // namespace opalx::spacecharge
