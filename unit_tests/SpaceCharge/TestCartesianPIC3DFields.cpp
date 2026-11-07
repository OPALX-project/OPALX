/**
 * @file TestCartesianPIC3DFields.cpp
 * @brief Tests of the CartesianPIC3D parts that move charges and fields between particles and
 * the mesh, and of the P3M short-range field.
 *
 * Covers:
 * - `ParticleMeshFieldTransfer`: depositing a bunch together with its image charges restores R,
 *   dt, and Q afterwards, and gathering a uniform mesh field returns it unchanged.
 * - `RelativisticFieldComposer`: at gamma = 1 a direct source passes through unchanged with
 *   B = 0, and the shifted-Green source is mirrored in z with the image signs.
 * - `P3MShortRangeInteraction`: the lab-frame pair field is finite and nonzero. The moving-source
 *   field of each pair matches the erfc pair field evaluated in the source rest frame and
 *   Lorentz-transformed to E and B, with each source's own momentum and per-particle charges. A
 *   pair that only one partner reaches within the rest-frame cutoff acts in one direction only,
 *   and sources at rest reproduce the electrostatic field with B = 0.
 *
 * Each test builds only the objects it needs on an 8^3 mesh on one rank; no Poisson solve runs.
 */

#include <gtest/gtest.h>

#include "PartBunch/BunchStateHandler.h"
#include "PartBunch/CartesianDomain.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"
#include "SpaceCharge/CartesianPIC3D/CartesianPIC3DFieldStorage.h"
#include "SpaceCharge/CartesianPIC3D/ParticleMeshFieldTransfer.h"
#include "SpaceCharge/CartesianPIC3D/RelativisticFieldComposer.h"
#include "SpaceCharge/Poisson/P3MAdapters.h"
#include "Utilities/Options.h"

#include <array>
#include <cmath>
#include <memory>

namespace opalx::spacecharge {
    namespace {

        class CartesianPIC3DFieldsTest : public ::testing::Test {
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

        TEST_F(CartesianPIC3DFieldsTest, TransferRestoresImageStateAndGathersVectorField) {
            Options::useQMAttributes = false;
            CartesianDomain<double, 3> domain(storage());
            CartesianPIC3DFieldStorage<double, 3> fields(domain);
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
            const auto originalCharge = Kokkos::create_mirror(particles.getQView());
            Kokkos::deep_copy(originalCharge, particles.getQView());

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

        TEST_F(CartesianPIC3DFieldsTest, ComposerHandlesDirectAndMirroredSources) {
            CartesianDomain<double, 3> domain(storage());
            CartesianPIC3DFieldStorage<double, 3> fields(domain);
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

        TEST_F(CartesianPIC3DFieldsTest, P3MShortRangeProducesFinitePairField) {
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

            P3MShortRangeInteraction interaction;
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

        // Field of a source charge at a target separated by @p separation in the frame of the
        // separation, from the P3M short-range kernel erfc(alpha r)/r with alpha = 2 / cutoff.
        Vector_t<double, 3> shortRangePairField(
                const Vector_t<double, 3>& separation, double charge, double cutoff) {
            const double r = std::sqrt(separation.dot(separation));
            if (r >= cutoff) {
                return Vector_t<double, 3>(0.0);
            }
            const double alpha = 2.0 / cutoff;
            const double radial =
                    2.0 * alpha * std::exp(-alpha * alpha * r * r) / (std::sqrt(Physics::pi) * r)
                    + (1.0 - std::erf(alpha * r)) / (r * r);
            return (charge / (4.0 * Physics::pi * Physics::epsilon_0) * radial / r) * separation;
        }

        // Switches per-particle charge storage on for one test and restores the previous mode.
        class QMAttributesScope {
        public:
            QMAttributesScope() : previous_m(Options::useQMAttributes) {
                Options::useQMAttributes = true;
            }
            ~QMAttributesScope() { Options::useQMAttributes = previous_m; }

        private:
            bool previous_m;
        };

        // Lab-frame E and B at a target separated by @p separation from a source moving with
        // @p momentum (beta gamma). The erfc kernel is evaluated at the separation in the source
        // rest frame and Lorentz-transformed to the lab, B = beta x E / c.
        void movingSourcePairField(
                const Vector_t<double, 3>& separation, const Vector_t<double, 3>& momentum,
                double charge, double cutoff, Vector_t<double, 3>& electric,
                Vector_t<double, 3>& magnetic) {
            using Vector       = Vector_t<double, 3>;
            const double gamma = std::sqrt(1.0 + momentum.dot(momentum));
            const Vector restSeparation =
                    separation + (momentum.dot(separation) / (gamma + 1.0)) * momentum;
            const Vector restField = shortRangePairField(restSeparation, charge, cutoff);
            electric = gamma * restField - (momentum.dot(restField) / (gamma + 1.0)) * momentum;
            magnetic = (1.0 / (gamma * Physics::c)) * cross(momentum, electric);
        }

        TEST_F(CartesianPIC3DFieldsTest, P3MMovingSourceFieldUsesEachSourceVelocity) {
            using Container     = ::ParticleContainer<double, 3>;
            using Vector        = Vector_t<double, 3>;
            auto setup          = storage();
            setup.layoutType    = ParticleLayoutType::SpatialOverlap;
            setup.overlapCutoff = 0.5;
            CartesianDomain<double, 3> domain(setup);
            const QMAttributesScope perParticleCharges;

            // Both partners inside each other's rest-frame cutoff, then a separation that only
            // the slower source reaches: s = 0.43 for pz = 0.4, but 0.62 for pz = 1.2.
            const std::array<Vector, 2> momenta{Vector(0.3, -0.1, 1.2), Vector(0.0, 0.0, 0.4)};
            const std::array<double, 2> charges{1.0e-15, 3.0e-15};
            for (const Vector& separation : {Vector(0.05, 0.02, 0.1), Vector(0.0, 0.0, 0.4)}) {
                Container particles(
                        domain.mesh(), domain.layout(), false,
                        Container::LayoutType::SpatialOverlap, setup.overlapCutoff, ippl::BC::NO);
                particles.setBunchStateHandler(std::make_shared<BunchStateHandler>());
                particles.createParticles(2);
                auto positions = particles.R.getHostMirror();
                auto momentum  = particles.P.getHostMirror();
                auto charge    = Kokkos::create_mirror_view(particles.getQView());
                positions(0)   = Vector(0.1, -0.05, -0.2);
                positions(1)   = positions(0) + separation;
                for (std::size_t i = 0; i < 2; ++i) {
                    momentum(i) = momenta[i];
                    charge(i)   = charges[i];
                }
                Kokkos::deep_copy(particles.R.getView(), positions);
                Kokkos::deep_copy(particles.P.getView(), momentum);
                Kokkos::deep_copy(particles.getQView(), charge);
                particles.update();
                ASSERT_EQ(particles.getLocalNum(), 2u);
                Kokkos::deep_copy(particles.E.getView(), Vector(0.0));
                Kokkos::deep_copy(particles.B.getView(), Vector(0.0));

                P3MShortRangeInteraction::addMovingSourceField(particles);

                // The layout sorts particles by cell, so read the order back.
                auto host = [](const auto& view) {
                    return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), view);
                };
                const auto sortedR  = host(particles.R.getView());
                const auto sortedP  = host(particles.P.getView());
                const auto sortedQ  = host(particles.getQView());
                const auto electric = host(particles.E.getView());
                const auto magnetic = host(particles.B.getView());
                std::size_t reached = 0;
                for (std::size_t i = 0; i < 2; ++i) {
                    const std::size_t j = 1 - i;
                    Vector expectedE, expectedB;
                    movingSourcePairField(
                            sortedR(i) - sortedR(j), sortedP(j), sortedQ(j), setup.overlapCutoff,
                            expectedE, expectedB);
                    const double scaleE = std::sqrt(expectedE.dot(expectedE));
                    const double scaleB = std::sqrt(expectedB.dot(expectedB));
                    reached += scaleE > 0.0 ? 1 : 0;
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_NEAR(electric(i)[d], expectedE[d], 1.0e-12 * scaleE)
                                << "particle " << i << ", component " << d;
                        EXPECT_NEAR(magnetic(i)[d], expectedB[d], 1.0e-12 * scaleB)
                                << "particle " << i << ", component " << d;
                    }
                }
                EXPECT_EQ(reached, separation[0] > 0.0 ? 2u : 1u);

                // Sources at rest give the electrostatic field of apply() and no B.
                Kokkos::deep_copy(particles.P.getView(), Vector(0.0));
                Kokkos::deep_copy(particles.E.getView(), Vector(0.0));
                Kokkos::deep_copy(particles.B.getView(), Vector(0.0));
                P3MShortRangeInteraction::addMovingSourceField(particles);
                const auto restingE = host(particles.E.getView());
                const auto restingB = host(particles.B.getView());
                Kokkos::deep_copy(particles.E.getView(), Vector(0.0));
                P3MShortRangeInteraction().apply(particles);
                const auto electrostatic = host(particles.E.getView());
                for (std::size_t i = 0; i < 2; ++i) {
                    const double scale = std::sqrt(electrostatic(i).dot(electrostatic(i)));
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_NEAR(restingE(i)[d], electrostatic(i)[d], 1.0e-12 * scale);
                        EXPECT_EQ(restingB(i)[d], 0.0);
                    }
                }
            }
        }

    }  // namespace
}  // namespace opalx::spacecharge
