/**
 * @file TestP3MShortRangeParallel.cpp
 * @brief Two-rank test of the P3M short-range field with a ghost partner.
 *
 * A close particle pair is placed across the boundary between two ranks, so each rank owns one
 * particle and sees the other only as a ghost copy. The partners carry different charges and
 * momenta. On each rank, the electrostatic and the moving-source short-range fields are compared
 * with the analytic pair field of the partner. The checks fail if the ghost's own charge or
 * momentum is not used, for example if only the local range of an attribute is read.
 *
 * The test needs exactly two ranks and is skipped otherwise, so CTest only runs it on one rank.
 * Run it by hand from the build directory:
 * `mpirun -n 2 ./unit_tests/SpaceCharge/TestP3MShortRangeParallel`
 */

#include <gtest/gtest.h>

#include "PartBunch/BunchStateHandler.h"
#include "PartBunch/CartesianDomain.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"
#include "SpaceCharge/Poisson/P3MAdapters.h"
#include "Utilities/Options.h"

#include <array>
#include <cmath>
#include <memory>

namespace opalx::spacecharge {
    namespace {

        using Container = ::ParticleContainer<double, 3>;
        using Vector    = Vector_t<double, 3>;

        class P3MShortRangeParallelTest : public ::testing::Test {
        protected:
            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
            }

            static void TearDownTestSuite() { ippl::finalize(); }

            void SetUp() override {
                if (ippl::Comm->size() != 2) {
                    GTEST_SKIP() << "Needs 2 MPI ranks for a pair across the rank boundary.";
                }
                previousQMAttributes_m   = Options::useQMAttributes;
                Options::useQMAttributes = true;
            }

            void TearDown() override { Options::useQMAttributes = previousQMAttributes_m; }

        private:
            bool previousQMAttributes_m = Options::useQMAttributes;
        };

        // Field of a source charge at a target separated by @p separation in the frame of the
        // separation, from the P3M short-range kernel erfc(alpha r)/r with alpha = 2 / cutoff.
        Vector shortRangePairField(const Vector& separation, double charge, double cutoff) {
            const double r = std::sqrt(separation.dot(separation));
            if (r >= cutoff) {
                return Vector(0.0);
            }
            const double alpha = 2.0 / cutoff;
            const double radial =
                    2.0 * alpha * std::exp(-alpha * alpha * r * r) / (std::sqrt(Physics::pi) * r)
                    + (1.0 - std::erf(alpha * r)) / (r * r);
            return (charge / (4.0 * Physics::pi * Physics::epsilon_0) * radial / r) * separation;
        }

        // Each rank owns one particle of the pair and sees the other only as a ghost, so the
        // ghost's charge and momentum must enter the pair field.
        TEST_F(P3MShortRangeParallelTest, GhostPartnerUsesItsChargeAndMomentum) {
            CartesianDomainConfig3D setup;
            setup.meshSize      = {8, 8, 8};
            setup.decomposition = {false, false, true};
            setup.layoutType    = ParticleLayoutType::SpatialOverlap;
            setup.overlapCutoff = 0.5;
            // The domain spans [-2, 2) m per axis, and the two ranks meet at z = 0.
            CartesianDomain<double, 3> domain(setup);
            Container particles(
                    domain.mesh(), domain.layout(), false, Container::LayoutType::SpatialOverlap,
                    setup.overlapCutoff, ippl::BC::NO);
            particles.setBunchStateHandler(std::make_shared<BunchStateHandler>());

            const std::array<Vector, 2> placed{Vector(0.1, 0.05, -0.05), Vector(0.15, 0.0, 0.1)};
            const std::array<Vector, 2> momenta{Vector(0.0, 0.1, 0.5), Vector(0.2, 0.0, 1.5)};
            const std::array<double, 2> charges{1.0e-15, 4.0e-15};
            particles.createParticles(ippl::Comm->rank() == 0 ? placed.size() : 0);
            if (ippl::Comm->rank() == 0) {
                auto positions = particles.R.getHostMirror();
                auto momentum  = particles.P.getHostMirror();
                auto charge    = Kokkos::create_mirror_view(particles.getQView());
                for (std::size_t i = 0; i < placed.size(); ++i) {
                    positions(i) = placed[i];
                    momentum(i)  = momenta[i];
                    charge(i)    = charges[i];
                }
                Kokkos::deep_copy(particles.R.getView(), positions);
                Kokkos::deep_copy(particles.P.getView(), momentum);
                Kokkos::deep_copy(particles.getQView(), charge);
            }
            particles.update();
            ASSERT_EQ(particles.getLocalNum(), 1u);

            const auto positions =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.R.getView());
            const Vector own           = positions(0);
            const std::size_t partner  = own[2] < 0.0 ? 1 : 0;
            const Vector separation    = own - placed[partner];
            const double partnerCharge = charges[partner];
            const Vector& u            = momenta[partner];

            // The electrostatic field uses the ghost's charge.
            Kokkos::deep_copy(particles.E.getView(), Vector(0.0));
            P3MShortRangeInteraction().apply(particles);
            const auto labField =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.E.getView());
            const Vector labExpected =
                    shortRangePairField(separation, partnerCharge, setup.overlapCutoff);
            const double labScale = std::sqrt(labExpected.dot(labExpected));
            EXPECT_GT(labScale, 0.0);
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_NEAR(labField(0)[d], labExpected[d], 1.0e-12 * labScale);
            }

            // The moving-source field uses the ghost's momentum: the erfc field at the separation
            // in the partner's rest frame, Lorentz-transformed to the lab.
            Kokkos::deep_copy(particles.E.getView(), Vector(0.0));
            Kokkos::deep_copy(particles.B.getView(), Vector(0.0));
            P3MShortRangeInteraction::addMovingSourceField(particles);
            const auto electric =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.E.getView());
            const auto magnetic =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.B.getView());
            const double gamma     = std::sqrt(1.0 + u.dot(u));
            const Vector restField = shortRangePairField(
                    separation + (u.dot(separation) / (gamma + 1.0)) * u, partnerCharge,
                    setup.overlapCutoff);
            const Vector expectedE = gamma * restField - (u.dot(restField) / (gamma + 1.0)) * u;
            const Vector expectedB = (1.0 / (gamma * Physics::c)) * cross(u, expectedE);
            const double scaleE    = std::sqrt(expectedE.dot(expectedE));
            const double scaleB    = std::sqrt(expectedB.dot(expectedB));
            EXPECT_GT(scaleE, 0.0);
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_NEAR(electric(0)[d], expectedE[d], 1.0e-12 * scaleE);
                EXPECT_NEAR(magnetic(0)[d], expectedB[d], 1.0e-12 * scaleB);
            }
        }

    }  // namespace
}  // namespace opalx::spacecharge
