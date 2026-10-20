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
#include <memory>
#include <vector>

extern Inform* gmsg;

namespace opalx::spacecharge {
    namespace {

        /**
         * ORB load balancing must cut only the axes the layout distributes (PARFFTX/Y/Z).
         *
         * ORB runs only on two or more ranks. CI runs this test on one rank (NUM_PROCS 1), where
         * it checks that a z-only layout survives a full updateForSolve with ORB due but skipped;
         * the ORB cut itself is not executed there. Run it with `mpirun -n 2` to check the cut.
         */
        class CartesianDomainUpdaterORBTest : public ::testing::Test {
        protected:
            using Container = ::ParticleContainer<double, 3>;

            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
                OpalData::getInstance()->storeInputFn("cartesian_domain_updater_orb.opal");
                gmsg                = new Inform(nullptr, -1);
                Options::enableHDF5 = false;
            }

            static void TearDownTestSuite() {
                delete gmsg;
                gmsg = nullptr;
                std::remove("cartesian_domain_updater_orb.stat");
                std::remove("cartesian_domain_updater_orb.lbal");
                ippl::finalize();
            }

            /**
             * Every rank creates the same particles: uniform in x and y, and z = -1 + 2 u^3 with u
             * uniform, so most particles sit at low z and the initial z split is imbalanced.
             */
            static void fillSkewedBeam(Container& particles, std::size_t count) {
                particles.createParticles(count);
                auto positions = particles.R.getHostMirror();
                auto momenta   = particles.P.getHostMirror();
                for (std::size_t i = 0; i < count; ++i) {
                    const double n  = static_cast<double>(count);
                    const double ux = (static_cast<double>(i) + 0.5) / n;
                    const double uy = (static_cast<double>((7 * i) % count) + 0.5) / n;
                    const double uz = (static_cast<double>((13 * i) % count) + 0.5) / n;
                    positions(i)    = Vector_t<double, 3>(
                            -1.0 + 2.0 * ux, -0.25 + 0.5 * uy, -1.0 + 2.0 * uz * uz * uz);
                    momenta(i) = Vector_t<double, 3>(0.0);
                }
                Kokkos::deep_copy(particles.R.getView(), positions);
                Kokkos::deep_copy(particles.P.getView(), momenta);
                particles.setM(Physics::m_e);
            }

            [[nodiscard]] static std::vector<ippl::NDIndex<3>> localDomains(
                    const ippl::FieldLayout<3>& layout) {
                const auto domains = layout.getHostLocalDomains();
                return {domains.data(), domains.data() + domains.extent(0)};
            }

            /** @brief Run one beam-frame update that is due for redistribution at any imbalance. */
            [[nodiscard]] static bool updateOnce(
                    const CartesianDomainConfig3D& domainConfig,
                    CartesianPIC3DFieldStorage<double, 3>& workspace, Container& particles) {
                PoissonSolverConfig poissonConfig;
                poissonConfig.type = PoissonSolverType::None;
                auto poisson       = makePoissonSolver(
                        poissonConfig, {&workspace.chargeDensity(), &workspace.electricField()});

                CartesianPIC3DConfig values;
                values.backend                    = PoissonSolverType::None;
                values.grid.meshSize              = domainConfig.meshSize;
                values.grid.decomposition         = domainConfig.decomposition;
                values.layoutRebuildDecomposition = domainConfig.decomposition;
                values.repartitionFrequency       = 1;
                values.loadBalancingThreshold     = 0.0;
                std::vector<Container*> containers{&particles};
                CartesianDomainUpdater updater(values, containers);

                SpaceChargeStepState step;
                step.mpiSize = ippl::Comm->size();
                const std::array<std::uint8_t, 1> activity{1};
                SpaceChargeSolveContext context(activity, step);
                return updater.updateForSolve(
                        DomainCoordinateFrame::Beam, context, {}, nullptr, workspace, *poisson);
            }
        };

        TEST_F(CartesianDomainUpdaterORBTest, CutsOnlyDistributedAxes) {
            // x is the longest axis, so an unrestricted ORB would cut it first.
            CartesianDomainConfig3D domainConfig;
            domainConfig.meshSize      = {32, 8, 16};
            domainConfig.decomposition = {false, false, true};
            CartesianDomain<double, 3> domain(domainConfig);
            CartesianPIC3DFieldStorage<double, 3> workspace(domain);
            workspace.initializeFields(PoissonSolverType::None);

            auto bunchState = std::make_shared<BunchStateHandler>();
            Container particles(domain.mesh(), domain.layout());
            particles.setBunchStateHandler(bunchState);
            constexpr std::size_t perRank = 256;
            fillSkewedBeam(particles, perRank);

            const bool redistributed = updateOnce(domainConfig, workspace, particles);
            if (ippl::Comm->size() < 2) {
                EXPECT_FALSE(redistributed) << "ORB must not run on one rank.";
            } else {
                ASSERT_TRUE(redistributed);
            }

            const auto& layout = domain.layout();
            const auto& global = layout.getDomain();
            for (const auto& local : localDomains(layout)) {
                EXPECT_EQ(local[0], global[0]);
                EXPECT_EQ(local[1], global[1]);
            }
            EXPECT_EQ(layout.isParallel(), domainConfig.decomposition);
            EXPECT_EQ(
                    particles.getTotalNum(),
                    perRank * static_cast<std::size_t>(ippl::Comm->size()));
        }

    }  // namespace
}  // namespace opalx::spacecharge
