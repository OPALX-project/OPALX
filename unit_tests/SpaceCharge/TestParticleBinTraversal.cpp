#include <gtest/gtest.h>

#include "PartBunch/BunchStateHandler.h"
#include "SpaceCharge/CartesianPIC3D/ParticleBinTraversal.h"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <vector>

namespace opalx::spacecharge {
    namespace {

        class ParticleBinTraversalTest : public ::testing::Test {
        protected:
            using Container = ::ParticleContainer<double, 3>;

            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
            }

            static void TearDownTestSuite() { ippl::finalize(); }

            void SetUp() override {
                ippl::NDIndex<3> indexDomain;
                for (unsigned dimension = 0; dimension < 3; ++dimension) {
                    indexDomain[dimension] = ippl::Index(8);
                }
                mesh_m = std::make_unique<Mesh_t<3>>(
                        indexDomain, Vector_t<double, 3>(0.125), Vector_t<double, 3>(0.0));
                layout_m = std::make_unique<FieldLayout_t<3>>(
                        MPI_COMM_WORLD, indexDomain,
                        std::array<bool, 3>{false, false, ippl::Comm->size() > 1}, true);
                particles_m = std::make_unique<Container>(*mesh_m, *layout_m);
                particles_m->setBunchStateHandler(std::make_shared<BunchStateHandler>());
                particles_m->createParticles(8);

                auto position = particles_m->R.getHostMirror();
                auto momentum = particles_m->P.getHostMirror();
                for (std::size_t index = 0; index < 8; ++index) {
                    position(index) = Vector_t<double, 3>(0.25, 0.25, 0.1 * (index + 1));
                    momentum(index) = Vector_t<double, 3>(0.0, 0.0, 0.1 * (index + 1));
                }
                Kokkos::deep_copy(particles_m->R.getView(), position);
                Kokkos::deep_copy(particles_m->P.getView(), momentum);
                particles_m->update();
            }

            void TearDown() override {
                particles_m.reset();
                layout_m.reset();
                mesh_m.reset();
            }

            [[nodiscard]] BinningConfig config(bool adaptive) const {
                BinningConfig values;
                values.name          = "PLAN_TEST";
                values.maximumBins   = 4;
                values.desiredWidth  = 0.15;
                values.alpha         = 1.0;
                values.beta          = 1.5;
                values.parameter     = BinningVariable::VelocityZ;
                values.adaptive      = adaptive;
                values.dumpFrequency = 1;
                return values;
            }

            std::unique_ptr<Mesh_t<3>> mesh_m;
            std::unique_ptr<FieldLayout_t<3>> layout_m;
            std::unique_ptr<Container> particles_m;
        };

        TEST_F(ParticleBinTraversalTest, RequiresPreparationBeforeLazyTraversal) {
            ParticleBinTraversal plan(*particles_m, config(false));
            EXPECT_THROW(static_cast<void>(plan.nextNonemptyBin()), OpalException);
        }

        TEST_F(ParticleBinTraversalTest, TraversesOnlyNonemptyBinsAndReturnsOptionalSnapshot) {
            ParticleBinTraversal plan(*particles_m, config(false));
            const BinPreparationResult prepared = plan.prepareBins(true);
            EXPECT_GT(prepared.mergedBinCount, 0u);
            EXPECT_TRUE(prepared.beforeMerge.has_value());
            EXPECT_FALSE(prepared.afterMerge.has_value());

            std::size_t visited       = 0;
            std::size_t particleCount = 0;
            while (const auto unit = plan.nextNonemptyBin()) {
                ++visited;
                particleCount += unit->globalParticleCount;
                EXPECT_GT(unit->globalParticleCount, 0u);
                EXPECT_GT(unit->gamma, 0.0);
            }
            EXPECT_LE(visited, prepared.mergedBinCount);
            EXPECT_EQ(particleCount, particles_m->getTotalNum());
            EXPECT_FALSE(plan.nextNonemptyBin().has_value());
        }

        TEST_F(ParticleBinTraversalTest, CapturesPostMergeSnapshotOnlyForAdaptiveBinning) {
            ParticleBinTraversal plan(*particles_m, config(true));
            const BinPreparationResult prepared = plan.prepareBins(true);
            EXPECT_TRUE(prepared.beforeMerge.has_value());
            EXPECT_TRUE(prepared.afterMerge.has_value());
            EXPECT_EQ(prepared.afterMerge->particleCounts.size(), prepared.mergedBinCount);
        }

        TEST_F(ParticleBinTraversalTest, SelectionsPartitionParticlesAndPreserveMoments) {
            auto momentum = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), particles_m->P.getView());
            for (std::size_t i = 0; i < particles_m->getLocalNum(); ++i) {
                momentum(i)[2] = (i % 2 == 0 ? 0.8 : 0.1) + 0.1 * ippl::Comm->rank();
                momentum(i)[0] = 2.0 * momentum(i)[2];
                momentum(i)[1] = -momentum(i)[2];
            }
            Kokkos::deep_copy(particles_m->P.getView(), momentum);

            ParticleBinTraversal plan(*particles_m, config(false));
            (void)plan.prepareBins(false);
            std::vector<bool> visited(particles_m->getLocalNum(), false);
            while (const auto unit = plan.nextNonemptyBin()) {
                const auto selection = unit->indexedSelection;
                const auto hash =
                        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), selection.hash());
                std::array<double, 4> expected{};
                for (auto i = selection.policy().begin(); i < selection.policy().end(); ++i) {
                    const auto index = static_cast<std::size_t>(hash(i));
                    EXPECT_LT(index, visited.size());
                    if (index >= visited.size()) continue;
                    EXPECT_FALSE(visited[index]);
                    visited[index] = true;
                    const auto p   = momentum(index);
                    for (unsigned d = 0; d < 3; ++d)
                        expected[d] += p[d];
                    expected[3] += std::sqrt(1.0 + p.dot(p));
                }
                for (double& component : expected) {
                    ippl::Comm->allreduce(component, 1, std::plus<double>());
                    component /= static_cast<double>(unit->globalParticleCount);
                }
                for (unsigned d = 0; d < 3; ++d) {
                    EXPECT_NEAR(unit->meanMomentum[d], expected[d], 1.0e-13);
                }
                EXPECT_NEAR(unit->gamma, expected[3], 1.0e-13);
            }
            EXPECT_TRUE(std::all_of(visited.begin(), visited.end(), [](bool seen) {
                return seen;
            }));
        }

    }  // namespace
}  // namespace opalx::spacecharge
