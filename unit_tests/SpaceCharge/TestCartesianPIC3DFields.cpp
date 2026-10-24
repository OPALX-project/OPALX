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
 * - `P3MShortRangeInteraction`: the lab-frame pair field is finite and nonzero. The rest-frame
 *   field matches the analytic erfc pair field at the separation (dx, dy, gamma dz), is exactly
 *   zero for a pair inside the lab cutoff but outside the rest-frame cutoff, and equals the
 *   lab-frame field at gamma = 1.
 * - `RelativisticFieldComposer::addParticleContribution`: per-particle fields get the same
 *   Lorentz conversion as mesh fields, including B = sign * v x E / c^2.
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
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <cmath>
#include <limits>
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

        TEST_F(CartesianPIC3DFieldsTest, P3MRestFrameFieldStretchesLongitudinalSeparation) {
            using Container     = ::ParticleContainer<double, 3>;
            using Vector        = Vector_t<double, 3>;
            auto setup          = storage();
            setup.layoutType    = ParticleLayoutType::SpatialOverlap;
            setup.overlapCutoff = 0.5;
            CartesianDomain<double, 3> domain(setup);
            const double charge = 1.0e-15;
            const double gamma  = 2.0;

            // Inside the rest-frame cutoff, then inside the lab cutoff only.
            for (const Vector& separation : {Vector(0.05, 0.02, 0.1), Vector(0.0, 0.0, 0.3)}) {
                Container particles(
                        domain.mesh(), domain.layout(), false,
                        Container::LayoutType::SpatialOverlap, setup.overlapCutoff, ippl::BC::NO);
                particles.setBunchStateHandler(std::make_shared<BunchStateHandler>());
                particles.createParticles(2);
                particles.setQ(charge);
                auto positions = particles.R.getHostMirror();
                positions(0)   = Vector(0.1, -0.05, 0.2);
                positions(1)   = positions(0) + separation;
                Kokkos::deep_copy(particles.R.getView(), positions);
                particles.update();
                ASSERT_EQ(particles.getLocalNum(), 2u);

                P3MShortRangeInteraction interaction;
                const auto restField = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), interaction.computeRestFrameField(particles, gamma));
                const auto sorted = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), particles.R.getView());
                for (std::size_t i = 0; i < 2; ++i) {
                    Vector restSeparation = sorted(i) - sorted(1 - i);
                    restSeparation[2] *= gamma;
                    const Vector expected =
                            shortRangePairField(restSeparation, charge, setup.overlapCutoff);
                    const double scale = std::sqrt(expected.dot(expected));
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_NEAR(restField(i)[d], expected[d], 1.0e-12 * scale);
                    }
                }

                // At gamma = 1 the rest-frame field is the electrostatic one.
                Kokkos::deep_copy(particles.E.getView(), Vector(0.0));
                interaction.apply(particles);
                const auto electrostatic = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), particles.E.getView());
                const auto unboosted = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), interaction.computeRestFrameField(particles, 1.0));
                for (std::size_t i = 0; i < 2; ++i) {
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_DOUBLE_EQ(unboosted(i)[d], electrostatic(i)[d]);
                    }
                }
            }

            Container particles(
                    domain.mesh(), domain.layout(), false, Container::LayoutType::SpatialOverlap,
                    setup.overlapCutoff, ippl::BC::NO);
            particles.setBunchStateHandler(std::make_shared<BunchStateHandler>());
            particles.createParticles(2);
            particles.update();
            P3MShortRangeInteraction interaction;
            EXPECT_THROW(interaction.computeRestFrameField(particles, 0.5), OpalException);
            EXPECT_THROW(
                    interaction.computeRestFrameField(
                            particles, std::numeric_limits<double>::quiet_NaN()),
                    OpalException);
        }

        TEST_F(CartesianPIC3DFieldsTest, ParticleContributionMatchesMeshConversion) {
            using Vector = Vector_t<double, 3>;
            CartesianDomain<double, 3> domain(storage());
            CartesianPIC3DFieldStorage<double, 3> fields(domain);
            fields.initializeFields(PoissonSolverType::Open);
            const Vector electricRest(1.0e3, -2.0e3, 5.0e2);
            fields.electricField() = electricRest;

            RelativisticFieldComposer composer;
            RelativisticFieldComposer::Policy policy;
            policy.meanMomentum = {0.3, -0.2, 1.5};
            policy.gamma        = std::sqrt(1.0 + 0.3 * 0.3 + 0.2 * 0.2 + 1.5 * 1.5);
            policy.magneticSign = -1.0;
            composer.clearAccumulation(fields);
            composer.accumulate(fields, policy);

            using Container = ::ParticleContainer<double, 3>;
            Container particles(domain.mesh(), domain.layout());
            particles.setBunchStateHandler(std::make_shared<BunchStateHandler>());
            const std::size_t count = 3;
            particles.createParticles(count);
            Kokkos::deep_copy(particles.E.getView(), Vector(0.0));
            Kokkos::deep_copy(particles.B.getView(), Vector(0.0));
            RelativisticFieldComposer::VectorView restField("rest-frame field", count);
            Kokkos::deep_copy(restField, electricRest);
            composer.addParticleContribution(restField, particles.E, particles.B, count, policy);
            Kokkos::fence();

            const int ghost  = fields.electricField().getNghost();
            const auto meshE = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), fields.accumulatedElectricField().getView());
            const auto meshB = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), fields.accumulatedMagneticField().getView());
            const auto particleE =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.E.getView());
            const auto particleB =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.B.getView());
            const Vector velocity =
                    Physics::c / policy.gamma
                    * Vector(
                            policy.meanMomentum[0], policy.meanMomentum[1], policy.meanMomentum[2]);
            for (std::size_t i = 0; i < count; ++i) {
                // B = sign * v x E_lab / c^2 holds for the full boost.
                const Vector expectedB = (policy.magneticSign / (Physics::c * Physics::c))
                                         * cross(velocity, Vector(particleE(i)));
                const double scaleB = std::sqrt(expectedB.dot(expectedB));
                for (unsigned d = 0; d < 3; ++d) {
                    EXPECT_DOUBLE_EQ(particleE(i)[d], meshE(ghost, ghost, ghost)[d]);
                    EXPECT_DOUBLE_EQ(particleB(i)[d], meshB(ghost, ghost, ghost)[d]);
                    EXPECT_NEAR(particleB(i)[d], expectedB[d], 1.0e-12 * scaleB);
                }
            }

            policy.sourceRule = FieldSourceRule::ShiftedGreenImageZ;
            EXPECT_THROW(
                    composer.addParticleContribution(
                            restField, particles.E, particles.B, count, policy),
                    OpalException);
        }

        TEST_F(CartesianPIC3DFieldsTest, CentreOfMomentumPreservesVelocityAndFourMomentum) {
            using V = Vector_t<double, 3>;
            const std::array<std::array<V, 2>, 5> samples{
                    {{V(0, 0, 12), V(0, 0, 12)},
                     {V(0.3, 0, 12), V(-0.3, 0, 12)},
                     {V(0, 0, 8), V(0, 0, 16)},
                     {V(0, 0, -12), V(0, 0, 12)},
                     {V(1.3, -0.7, 12), V(1.7, -0.3, 8)}}};
            for (const auto& sample : samples) {
                const V mean        = (sample[0] + sample[1]) / 2.0;
                const double energy = (std::sqrt(1 + sample[0].dot(sample[0]))
                                       + std::sqrt(1 + sample[1].dot(sample[1])))
                                      / 2.0;
                const auto policy = FieldCompositionPolicy::centreOfMomentum(
                        {mean[0], mean[1], mean[2]}, energy);
                const auto boost = LorentzBoost::from(policy);
                const V beta     = mean / energy;
                EXPECT_NEAR(boost.gamma * boost.gamma * (1 - beta.dot(beta)), 1.0, 1e-12);
                V totalRestMomentum(0.0);
                for (const auto& p : sample) {
                    totalRestMomentum +=
                            p + (boost.gamma - 1) * p.dot(boost.velocityUnit) * boost.velocityUnit
                            - boost.gamma * std::sqrt(1 + p.dot(p)) * beta;
                }
                for (unsigned d = 0; d < 3; ++d) {
                    EXPECT_NEAR(boost.velocity[d] / Physics::c, beta[d], 1e-15);
                    EXPECT_NEAR(totalRestMomentum[d], 0.0, 1e-12);
                }
                const V restE(1.3, -0.7, 0.4);
                V e, b;
                boost.toLab(restE, e, b);
                const V cb = Physics::c * b;
                EXPECT_NEAR(e.dot(e) - cb.dot(cb), restE.dot(restE), 1e-10);
                EXPECT_NEAR(e.dot(cb), 0.0, 1e-10);
                const V force    = e + cross(beta, cb);
                const V parallel = restE.dot(boost.velocityUnit) * boost.velocityUnit;
                const V expected = parallel + (restE - parallel) / boost.gamma;
                for (unsigned d = 0; d < 3; ++d)
                    EXPECT_NEAR(force[d], expected[d], 1e-12);
            }
        }

        TEST_F(CartesianPIC3DFieldsTest, StationaryWarmCMHasExactUnitBoost) {
            const double meanGamma = std::sqrt(226.0);  // Equal populations with pz = +/-15.
            const auto policy      = FieldCompositionPolicy::centreOfMomentum({0, 0, 0}, meanGamma);
            const auto boost       = LorentzBoost::from(policy);
            EXPECT_EQ(policy.gamma, 1.0);
            EXPECT_DOUBLE_EQ(policy.meanParticleGamma, meanGamma);
            const Vector_t<double, 3> restE(1.3, -0.7, 0.4);
            Vector_t<double, 3> e, b;
            boost.toLab(restE, e, b);
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_DOUBLE_EQ(boost.velocity[d], 0.0);
                EXPECT_DOUBLE_EQ(e[d], restE[d]);
                EXPECT_DOUBLE_EQ(b[d], 0.0);
            }
            for (const double meanPz : {1e-12, 1e-9}) {
                SCOPED_TRACE(meanPz);
                const auto nearRest =
                        FieldCompositionPolicy::centreOfMomentum({0, 0, meanPz}, meanGamma);
                // These velocities are too small for gamma-1 to be representable in binary64.
                EXPECT_EQ(nearRest.gamma, 1.0);
                EXPECT_DOUBLE_EQ(nearRest.meanParticleGamma, meanGamma);
                EXPECT_NEAR(
                        LorentzBoost::from(nearRest).velocity[2] / Physics::c, meanPz / meanGamma,
                        1e-24);
            }
        }

        TEST_F(CartesianPIC3DFieldsTest, CentreOfMomentumRejectsInvalidMoments) {
            EXPECT_THROW(
                    static_cast<void>(FieldCompositionPolicy::centreOfMomentum({0, 0, 2}, 1)),
                    OpalException);
            EXPECT_THROW(
                    static_cast<void>(FieldCompositionPolicy::centreOfMomentum({0, 0, 0}, 0.5)),
                    OpalException);
            EXPECT_THROW(
                    static_cast<void>(FieldCompositionPolicy::centreOfMomentum(
                            {0, 0, std::numeric_limits<double>::infinity()}, 2)),
                    OpalException);
            EXPECT_THROW(
                    static_cast<void>(FieldCompositionPolicy::centreOfMomentum(
                            {0, 0, 0}, std::numeric_limits<double>::quiet_NaN())),
                    OpalException);
        }

        TEST_F(CartesianPIC3DFieldsTest, CMChargeNormalizationConservesChargeOnStretchedCells) {
            Options::useQMAttributes = false;
            CartesianDomain<double, 3> domain(storage());
            CartesianPIC3DFieldStorage<double, 3> fields(domain);
            fields.initializeFields(PoissonSolverType::Open);
            using Container = ::ParticleContainer<double, 3>;
            Container particles(domain.mesh(), domain.layout());
            particles.setBunchStateHandler(std::make_shared<BunchStateHandler>());
            particles.createParticles(2);
            particles.setQ(2e-15);
            Kokkos::deep_copy(particles.R.getView(), (domain.lower() + domain.upper()) / 2.0);
            Kokkos::deep_copy(particles.dt.getView(), 1e-12);
            const double g   = std::sqrt(1 + 0.3 * 0.3 + 12 * 12);
            const auto frame = FieldCompositionPolicy::centreOfMomentum({0, 0, 12}, g);
            ParticleMeshFieldTransfer transfer;
            ParticleMeshFieldTransfer::ChargeNormalization normalization;
            normalization.timeStep              = 1e-12;
            normalization.gamma                 = frame.gamma;
            normalization.selectedCharge        = 4e-15;
            normalization.couplingConstant      = 1.0;
            normalization.normalizeByCellVolume = true;
            transfer.depositCharge(
                    particles, fields, ParticleMeshFieldTransfer::DepositKind::Primary,
                    ParticleMeshFieldTransfer::Selection::direct(0, 2), normalization, {});
            const auto h                = fields.spacing();
            const double restCellVolume = h[0] * h[1] * h[2] * frame.gamma;
            EXPECT_NEAR(fields.chargeDensity().sum() * restCellVolume, 4e-15, 4e-27);
            const double densitySum = fields.chargeDensity().sum();
            normalization.gamma     = 1.0;
            transfer.depositCharge(
                    particles, fields, ParticleMeshFieldTransfer::DepositKind::Primary,
                    ParticleMeshFieldTransfer::Selection::direct(0, 2), normalization, {});
            EXPECT_NEAR(
                    densitySum * frame.gamma, fields.chargeDensity().sum(),
                    std::abs(fields.chargeDensity().sum()) * 1e-12);
            const auto dt = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), particles.dt.getView());
            EXPECT_DOUBLE_EQ(dt(0), 1e-12);
        }

    }  // namespace
}  // namespace opalx::spacecharge
