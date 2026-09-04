#include <gtest/gtest.h>

#include "PartBunch/BunchStateHandler.h"
#include "PartBunch/CartesianDomain.h"
#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/CartesianPIC/CartesianPICFieldStorage.h"
#include "SpaceCharge/CartesianPIC/P3MShortRangeInteraction.h"
#include "SpaceCharge/CartesianPIC/ParticleMeshFieldTransfer.h"
#include "SpaceCharge/CartesianPIC/RelativisticFieldComposer.h"
#include "Utilities/Options.h"

#include <cmath>
#include <memory>

namespace opalx::spacecharge {
    namespace {

        class CartesianPICFieldsTest : public ::testing::Test {
        protected:
            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
            }

            static void TearDownTestSuite() { ippl::finalize(); }

            CartesianDomainConfig3D storage() const {
                CartesianDomainConfig3D result;
                result.meshSize      = {8, 8, 8};
                result.decomposition = {false, false, false};
                return result;
            }
        };

        TEST_F(CartesianPICFieldsTest, TransferRestoresImageStateAndGathersVectorField) {
            Options::useQMAttributes = false;
            CartesianDomain<double, 3> domain(storage());
            CartesianPICFieldStorage<double, 3> fields(domain);
            fields.initializeFields(PoissonSolverType::None);

            using Container = ::ParticleContainer<double, 3>;
            Container particles(domain.mesh(), domain.layout());
            particles.setBunchStateHandler(std::make_shared<BunchStateHandler>());
            particles.createParticles(2);
            particles.setQ(2.0e-15);
            auto positions = particles.R.getHostMirror();
            auto timeSteps = particles.dt.getHostMirror();
            positions(0)   = Vector_t<double, 3>(-0.2, 0.1, 0.3);
            positions(1)   = Vector_t<double, 3>(0.25, -0.15, 0.45);
            timeSteps(0)   = 1.0e-12;
            timeSteps(1)   = 1.0e-12;
            Kokkos::deep_copy(particles.R.getView(), positions);
            Kokkos::deep_copy(particles.dt.getView(), timeSteps);
            const auto originalCharge =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.getQView());

            ParticleMeshFieldTransfer transfer;
            ParticleMeshFieldTransfer::ChargeNormalization normalization;
            normalization.timeStep       = 1.0e-12;
            normalization.gamma          = 1.0;
            normalization.selectedCharge = 4.0e-15;
            transfer.depositCharge(
                    particles, fields, ParticleMeshFieldTransfer::DepositKind::PrimaryAndImage,
                    ParticleMeshFieldTransfer::Selection::direct(0, 2), normalization,
                    {.enabled = true, .planeZ = 0.0});

            const auto restoredPositions =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.R.getView());
            const auto restoredTimeSteps = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), particles.dt.getView());
            const auto restoredCharge =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.getQView());
            for (std::size_t index = 0; index < 2; ++index) {
                for (unsigned dimension = 0; dimension < 3; ++dimension) {
                    EXPECT_DOUBLE_EQ(
                            restoredPositions(index)[dimension], positions(index)[dimension]);
                }
                EXPECT_DOUBLE_EQ(restoredTimeSteps(index), timeSteps(index));
            }
            EXPECT_DOUBLE_EQ(restoredCharge(0), originalCharge(0));

            const auto chargeDensity = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), fields.chargeDensity().getView());
            double absoluteChargeDensity = 0.0;
            for (std::size_t i = 0; i < chargeDensity.extent(0); ++i) {
                for (std::size_t j = 0; j < chargeDensity.extent(1); ++j) {
                    for (std::size_t k = 0; k < chargeDensity.extent(2); ++k) {
                        absoluteChargeDensity += std::abs(chargeDensity(i, j, k));
                    }
                }
            }
            EXPECT_GT(absoluteChargeDensity, 0.0);

            fields.electricField() = Vector_t<double, 3>(1.0, 2.0, 3.0);
            transfer.gatherVector(
                    particles.E, fields.electricField(), particles.R,
                    ParticleMeshFieldTransfer::GatherMode::Replace);
            const auto gathered =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.E.getView());
            for (std::size_t index = 0; index < 2; ++index) {
                for (unsigned dimension = 0; dimension < 3; ++dimension) {
                    EXPECT_NEAR(
                            gathered(index)[dimension], static_cast<double>(dimension + 1),
                            1.0e-14);
                }
            }
        }

        TEST_F(CartesianPICFieldsTest, ComposerHandlesDirectAndMirroredSources) {
            CartesianDomain<double, 3> domain(storage());
            CartesianPICFieldStorage<double, 3> fields(domain);
            fields.initializeFields(PoissonSolverType::Open);

            auto source = fields.electricField().getHostMirror();
            for (std::size_t i = 0; i < source.extent(0); ++i) {
                for (std::size_t j = 0; j < source.extent(1); ++j) {
                    for (std::size_t k = 0; k < source.extent(2); ++k) {
                        source(i, j, k) = Vector_t<double, 3>(
                                static_cast<double>(i + 1), static_cast<double>(j + 2),
                                static_cast<double>(k + 3));
                    }
                }
            }
            Kokkos::deep_copy(fields.electricField().getView(), source);

            RelativisticFieldComposer composer;
            RelativisticFieldComposer::Policy policy;
            policy.gamma = 1.0;
            composer.clearAccumulation(fields);
            composer.accumulate(fields, policy);
            Kokkos::fence();

            const int ghost           = fields.electricField().getNghost();
            const auto directElectric = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), fields.accumulatedElectricField().getView());
            const auto directMagnetic = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), fields.accumulatedMagneticField().getView());
            for (unsigned dimension = 0; dimension < 3; ++dimension) {
                EXPECT_DOUBLE_EQ(
                        directElectric(ghost, ghost, ghost)[dimension],
                        source(ghost, ghost, ghost)[dimension]);
                EXPECT_DOUBLE_EQ(directMagnetic(ghost, ghost, ghost)[dimension], 0.0);
            }

            composer.clearAccumulation(fields);
            policy.sourceRule = FieldSourceRule::ShiftedGreenImageZ;
            composer.accumulate(fields, policy);
            Kokkos::fence();
            const auto mirrored = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), fields.accumulatedElectricField().getView());
            const std::size_t mirroredK =
                    static_cast<std::size_t>(ghost) + fields.layoutExtents()[2] - 1;
            EXPECT_DOUBLE_EQ(mirrored(ghost, ghost, ghost)[0], -source(ghost, ghost, mirroredK)[0]);
            EXPECT_DOUBLE_EQ(mirrored(ghost, ghost, ghost)[1], -source(ghost, ghost, mirroredK)[1]);
            EXPECT_DOUBLE_EQ(mirrored(ghost, ghost, ghost)[2], source(ghost, ghost, mirroredK)[2]);
        }

        TEST_F(CartesianPICFieldsTest, P3MShortRangeProducesFinitePairField) {
            auto setup          = storage();
            setup.layoutType    = ParticleLayoutType::SpatialOverlap;
            setup.overlapCutoff = 0.5;
            CartesianDomain<double, 3> domain(setup);
            using Container = ::ParticleContainer<double, 3>;
            Container particles(
                    domain.mesh(), domain.layout(), false, Container::LayoutType::SpatialOverlap,
                    setup.overlapCutoff, ippl::BC::NO);
            particles.setBunchStateHandler(std::make_shared<BunchStateHandler>());
            particles.createParticles(2);
            particles.setQ(1.0e-15);
            auto positions = particles.R.getHostMirror();
            positions(0)   = Vector_t<double, 3>(-0.1, 0.0, 0.0);
            positions(1)   = Vector_t<double, 3>(0.1, 0.0, 0.0);
            Kokkos::deep_copy(particles.R.getView(), positions);
            particles.E = Vector_t<double, 3>(0.0);
            particles.update();

            P3MShortRangeInteraction interaction(setup.overlapCutoff);
            interaction.apply(particles);
            const auto electric =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.E.getView());
            double magnitude = 0.0;
            for (std::size_t index = 0; index < particles.getLocalNum(); ++index) {
                for (unsigned dimension = 0; dimension < 3; ++dimension) {
                    EXPECT_TRUE(std::isfinite(electric(index)[dimension]));
                    magnitude += std::abs(electric(index)[dimension]);
                }
            }
            EXPECT_GT(magnitude, 0.0);
        }

    }  // namespace
}  // namespace opalx::spacecharge
