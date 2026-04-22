/**
 * @file TestParticleContainer.cpp
 * @brief Unit tests for `ParticleContainer<double, 3>`.
 *
 * Covers:
 * - Charge / mass storage in both `SingleValue` and `Attributes` QM modes.
 * - Round-trip `scaleDtByCharge` / `unscaleDtByCharge` in both QM modes.
 * - Moment computation stability (empty and known-data cases).
 * - Min / max position computation.
 * - `deleteParticlesOutside` edge cases and known-outlier removal.
 *
 * The fixture constructs a lightweight `ParticleContainer` on a minimal 8^3
 * periodic mesh without a full `PartBunch`, so no field-solver or DataSink
 * infrastructure is required.
 */

#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <cmath>
#include <memory>
#include <vector>

#include "Ippl.h"
#include "PartBunch/FieldContainer.hpp"
#include "PartBunch/ParticleContainer.hpp"
#include "Utilities/Options.h"

namespace {

    using PC_t = ParticleContainer<double, 3>;

    class ParticleContainerTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }

        static void TearDownTestSuite() { ippl::finalize(); }

        /// Build a minimal ParticleContainer on an 8^3 periodic mesh over [-4, 4]^3.
        /// Reads `Options::useQMAttributes` at construction time.
        std::shared_ptr<PC_t> makeContainer() {
            ippl::Vector<int, 3> nr        = 8;
            ippl::Vector<double, 3> rmin   = -4.0;
            ippl::Vector<double, 3> rmax   = 4.0;
            ippl::Vector<double, 3> origin = rmin;
            ippl::Vector<double, 3> hr     = (rmax - rmin) / ippl::Vector<double, 3>(nr);
            std::array<bool, 3> decomp     = {true, true, true};

            ippl::NDIndex<3> domain;
            for (unsigned i = 0; i < 3; i++) {
                domain[i] = ippl::Index(nr[i]);
            }

            Mesh_t<3> mesh(domain, hr, origin);
            FieldLayout_t<3> fl(MPI_COMM_WORLD, domain, decomp, true);

            std::shared_ptr<PC_t> pc = std::make_shared<PC_t>(mesh, fl);
            pc->setBunchStateHandler(std::make_shared<BunchStateHandler>());
            return pc;
        }

        /// Create local particles at the given positions with uniform dt and P_z.
        void createParticlesAt(
                std::shared_ptr<PC_t>& pc, const std::vector<std::array<double, 3>>& positions,
                double dtVal = 1e-12, double pz = 0.1) {
            const size_t n = positions.size();
            if (n == 0) return;

            pc->create(n);

            auto R_host  = pc->R.getHostMirror();
            auto P_host  = pc->P.getHostMirror();
            auto dt_host = pc->dt.getHostMirror();

            for (size_t i = 0; i < n; ++i) {
                R_host(i)[0] = positions[i][0];
                R_host(i)[1] = positions[i][1];
                R_host(i)[2] = positions[i][2];
                P_host(i)[0] = 0.0;
                P_host(i)[1] = 0.0;
                P_host(i)[2] = pz;
                dt_host(i)   = dtVal;
            }

            Kokkos::deep_copy(pc->R.getView(), R_host);
            Kokkos::deep_copy(pc->P.getView(), P_host);
            Kokkos::deep_copy(pc->dt.getView(), dt_host);
            Kokkos::fence();
        }
    };

    // ================================================================
    // Charge / Mass – SingleValue mode
    // ================================================================

    TEST_F(ParticleContainerTest, ChargeMass_SingleValueMode_SetAndRead) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();

        ASSERT_EQ(pc->getQMStorageMode(), PC_t::QMStorageMode::SingleValue);

        const double qExpected = 1.6e-19;
        const double mExpected = 0.938272;

        pc->setQ(qExpected);
        pc->setM(mExpected);

        auto qView = pc->getQView();
        auto mView = pc->getMView();

        ASSERT_EQ(qView.extent(0), 1u);
        ASSERT_EQ(mView.extent(0), 1u);

        auto q_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qView);
        auto m_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), mView);

        EXPECT_DOUBLE_EQ(q_host(0), qExpected);
        EXPECT_DOUBLE_EQ(m_host(0), mExpected);
    }

    // ================================================================
    // Charge / Mass – Attributes mode
    // ================================================================

    TEST_F(ParticleContainerTest, ChargeMass_AttributeMode_SetAndRead) {
        Options::useQMAttributes = true;
        auto pc                  = makeContainer();

        ASSERT_EQ(pc->getQMStorageMode(), PC_t::QMStorageMode::Attributes);

        constexpr size_t nPart = 8;
        pc->create(nPart);
        Kokkos::fence();

        const double qExpected = -1.6e-19;
        const double mExpected = 0.000511;

        pc->setQ(qExpected);
        pc->setM(mExpected);

        auto q_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->getQView());
        auto m_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->getMView());

        for (size_t i = 0; i < nPart; ++i) {
            EXPECT_DOUBLE_EQ(q_host(i), qExpected) << "particle " << i;
            EXPECT_DOUBLE_EQ(m_host(i), mExpected) << "particle " << i;
        }
    }

    // ================================================================
    // dt scaling round-trip – SingleValue mode
    // ================================================================

    TEST_F(ParticleContainerTest, DtScaling_RoundTrip_SingleValueMode) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();

        const double q = 1.6e-19;
        pc->setQ(q);

        constexpr size_t nPart = 16;
        const double dtOrig    = 1e-12;
        std::vector<std::array<double, 3>> positions(nPart, {0.0, 0.0, 0.0});
        createParticlesAt(pc, positions, dtOrig);

        pc->scaleDtByCharge();

        auto dt_scaled = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->dt.getView());
        for (size_t i = 0; i < nPart; ++i) {
            EXPECT_NEAR(dt_scaled(i), dtOrig * q, 1e-44) << "particle " << i;
        }

        pc->unscaleDtByCharge();

        auto dt_back = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->dt.getView());
        for (size_t i = 0; i < nPart; ++i) {
            EXPECT_NEAR(dt_back(i), dtOrig, 1e-20) << "particle " << i;
        }
    }

    // ================================================================
    // dt scaling round-trip – Attributes mode
    // ================================================================

    TEST_F(ParticleContainerTest, DtScaling_RoundTrip_AttributeMode) {
        Options::useQMAttributes = true;
        auto pc                  = makeContainer();

        constexpr size_t nPart = 16;
        const double dtOrig    = 1e-12;
        std::vector<std::array<double, 3>> positions(nPart, {0.0, 0.0, 0.0});
        createParticlesAt(pc, positions, dtOrig);

        const double q = 1.6e-19;
        pc->setQ(q);

        pc->scaleDtByCharge();

        auto dt_scaled = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->dt.getView());
        for (size_t i = 0; i < nPart; ++i) {
            EXPECT_NEAR(dt_scaled(i), dtOrig * q, 1e-44) << "particle " << i;
        }

        pc->unscaleDtByCharge();

        auto dt_back = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->dt.getView());
        for (size_t i = 0; i < nPart; ++i) {
            EXPECT_NEAR(dt_back(i), dtOrig, 1e-20) << "particle " << i;
        }
    }

    // ================================================================
    // Moments & MinMax – empty container (no NaN / crash)
    // ================================================================

    TEST_F(ParticleContainerTest, MomentsAndMinMax_EmptyContainer_NoNaN) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();

        ASSERT_EQ(pc->getLocalNum(), 0u);

        pc->updateMoments();
        pc->computeMinMaxR();

        auto meanR = pc->getMeanR();
        auto rmsR  = pc->getRmsR();
        auto minR  = pc->getMinR();
        auto maxR  = pc->getMaxR();

        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_TRUE(std::isfinite(meanR[d])) << "meanR[" << d << "]";
            EXPECT_TRUE(std::isfinite(rmsR[d])) << "rmsR[" << d << "]";
            EXPECT_TRUE(std::isfinite(minR[d])) << "minR[" << d << "]";
            EXPECT_TRUE(std::isfinite(maxR[d])) << "maxR[" << d << "]";
            EXPECT_LE(minR[d], maxR[d]) << "min > max at dim " << d;
        }
    }

    // ================================================================
    // Moments – known symmetric data → mean ≈ 0, rms ≈ 1
    // ================================================================

    TEST_F(ParticleContainerTest, Moments_KnownSymmetricData) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();

        std::vector<std::array<double, 3>> positions;
        if (ippl::Comm->rank() == 0) {
            positions = {{1.0, 1.0, 1.0}, {-1.0, -1.0, -1.0}, {1.0, -1.0, 1.0}, {-1.0, 1.0, -1.0}};
        }
        createParticlesAt(pc, positions);
        pc->setM(1.0);

        pc->updateMoments();

        auto meanR = pc->getMeanR();
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_NEAR(meanR[d], 0.0, 1e-14) << "meanR[" << d << "]";
        }

        auto rmsR = pc->getRmsR();
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_NEAR(rmsR[d], 1.0, 1e-14) << "rmsR[" << d << "]";
        }
    }

    // ================================================================
    // MinMax R – known particle positions
    // ================================================================

    TEST_F(ParticleContainerTest, MinMaxR_KnownPositions) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();

        std::vector<std::array<double, 3>> positions;
        if (ippl::Comm->rank() == 0) {
            positions = {{1.0, 2.0, 3.0}, {-0.5, -1.5, 0.0}, {0.0, 0.0, -2.0}};
        }
        createParticlesAt(pc, positions);

        pc->computeMinMaxR();

        auto minR = pc->getMinR();
        auto maxR = pc->getMaxR();

        EXPECT_NEAR(minR[0], -0.5, 1e-14);
        EXPECT_NEAR(minR[1], -1.5, 1e-14);
        EXPECT_NEAR(minR[2], -2.0, 1e-14);

        EXPECT_NEAR(maxR[0], 1.0, 1e-14);
        EXPECT_NEAR(maxR[1], 2.0, 1e-14);
        EXPECT_NEAR(maxR[2], 3.0, 1e-14);
    }

    // ================================================================
    // deleteParticlesOutside – edge cases
    // ================================================================

    TEST_F(ParticleContainerTest, DeleteParticlesOutside_NegativeSigma_ReturnsZero) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();

        // Every rank gets a particle to avoid mismatched collectives inside
        // deleteParticlesOutside (getTotalNum is an allreduce).
        std::vector<std::array<double, 3>> positions = {{{0.0, 0.0, 0.0}}};
        createParticlesAt(pc, positions);
        pc->setM(1.0);

        size_t localBefore = pc->getLocalNum();
        auto destroyed     = pc->deleteParticlesOutside(-1.0);

        EXPECT_EQ(destroyed, 0u);
        EXPECT_EQ(pc->getLocalNum(), localBefore);
    }

    TEST_F(ParticleContainerTest, DeleteParticlesOutside_ZeroSigma_ReturnsZero) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();

        std::vector<std::array<double, 3>> positions = {{{0.0, 0.0, 0.0}}};
        createParticlesAt(pc, positions);
        pc->setM(1.0);

        size_t localBefore = pc->getLocalNum();
        auto destroyed     = pc->deleteParticlesOutside(0.0);

        EXPECT_EQ(destroyed, 0u);
        EXPECT_EQ(pc->getLocalNum(), localBefore);
    }

    TEST_F(ParticleContainerTest, DeleteParticlesOutside_EmptyContainer_ReturnsZero) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();

        ASSERT_EQ(pc->getLocalNum(), 0u);
        auto destroyed = pc->deleteParticlesOutside(3.0);
        EXPECT_EQ(destroyed, 0u);
    }

    // ================================================================
    // deleteParticlesOutside – known outlier removal
    // ================================================================

    TEST_F(ParticleContainerTest, DeleteParticlesOutside_KnownOutlierRemoval) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();

        // 10 cluster particles near origin + 1 outlier on rank 0.
        // Non-rank-0 processes get a single particle at the origin so that all
        // ranks satisfy nLocal > 0 (avoids mismatched collective in getTotalNum).
        std::vector<std::array<double, 3>> positions;
        if (ippl::Comm->rank() == 0) {
            positions = {
                    {0.10, 0.10, 0.10},    {-0.10, -0.10, -0.10}, {0.05, 0.05, 0.05},
                    {-0.05, -0.05, -0.05}, {0.10, -0.10, 0.05},   {-0.10, 0.10, -0.05},
                    {0.08, 0.03, 0.07},    {-0.08, -0.03, -0.07}, {0.02, 0.09, 0.04},
                    {-0.02, -0.09, -0.04}, {3.50, 3.50, 3.50}  // outlier – well beyond 2σ
            };
        } else {
            positions = {{{0.0, 0.0, 0.0}}};
        }
        createParticlesAt(pc, positions);
        pc->setM(1.0);

        size_t totalBefore = pc->getTotalNum();
        auto destroyed     = pc->deleteParticlesOutside(2.0);

        EXPECT_EQ(destroyed, 1u);

        size_t totalAfter = pc->getTotalNum();
        EXPECT_EQ(totalAfter + destroyed, totalBefore);
    }

}  // namespace
