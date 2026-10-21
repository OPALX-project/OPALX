/**
 * @file TestP3MRestFrameParallel.cpp
 * @brief Two-rank test of the P3M short-range field in the bunch rest frame.
 *
 * A close particle pair is placed across the boundary between two ranks, so each rank owns one
 * particle and sees the other only as a ghost copy. On each rank, the lab-frame and the
 * rest-frame short-range fields are compared with the analytic pair field. The rest-frame check
 * fails if only the local positions, and not the ghost positions, are stretched by gamma.
 *
 * The test needs exactly two ranks and is skipped otherwise, so CTest only runs it on one rank.
 * Run it by hand from the build directory:
 * `mpirun -n 2 ./unit_tests/SpaceCharge/TestP3MRestFrameParallel`
 */

#include <gtest/gtest.h>

#include "PartBunch/BunchStateHandler.h"
#include "PartBunch/CartesianDomain.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"
#include "SpaceCharge/Poisson/P3MAdapters.h"

#include <array>
#include <cmath>
#include <memory>

namespace opalx::spacecharge {
    namespace {

        using Container = ::ParticleContainer<double, 3>;
        using Vector    = Vector_t<double, 3>;

        class P3MRestFrameParallelTest : public ::testing::Test {
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
            }
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
        // ghost positions must receive the same longitudinal stretch as the local ones.
        TEST_F(P3MRestFrameParallelTest, GhostPartnerUsesRestFrameSeparation) {
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
            const double charge = 1.0e-15;
            const double gamma  = 2.0;
            particles.createParticles(ippl::Comm->rank() == 0 ? placed.size() : 0);
            particles.setQ(charge);
            if (ippl::Comm->rank() == 0) {
                auto positions = particles.R.getHostMirror();
                for (std::size_t i = 0; i < placed.size(); ++i)
                    positions(i) = placed[i];
                Kokkos::deep_copy(particles.R.getView(), positions);
            }
            particles.update();
            ASSERT_EQ(particles.getLocalNum(), 1u);

            P3MShortRangeInteraction interaction;
            const auto positions =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.R.getView());
            const Vector own      = positions(0);
            const Vector& partner = own[2] < 0.0 ? placed[1] : placed[0];

            // The electrostatic field confirms that the ghost partner is in place.
            Kokkos::deep_copy(particles.E.getView(), Vector(0.0));
            interaction.apply(particles);
            const auto labField =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.E.getView());
            const Vector labExpected =
                    shortRangePairField(own - partner, charge, setup.overlapCutoff);
            const double labScale = std::sqrt(labExpected.dot(labExpected));
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_NEAR(labField(0)[d], labExpected[d], 1.0e-12 * labScale);
            }

            const auto restField = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), interaction.computeRestFrameField(particles, gamma));
            Vector restSeparation = own - partner;
            restSeparation[2] *= gamma;
            const Vector expected =
                    shortRangePairField(restSeparation, charge, setup.overlapCutoff);
            const double scale = std::sqrt(expected.dot(expected));
            EXPECT_GT(scale, 0.0);
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_NEAR(restField(0)[d], expected[d], 1.0e-12 * scale);
            }
        }

    }  // namespace
}  // namespace opalx::spacecharge
