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
#include "PartBunch/ImageChargeScatterController.h"
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

            pc->createParticles(n);

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
        pc->createParticles(nPart);
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

    TEST_F(ParticleContainerTest, ImageChargeMirrorTransform_AllParticles_RoundTrip) {
        Options::useQMAttributes                     = false;
        auto pc                                      = makeContainer();
        std::vector<std::array<double, 3>> positions = {
                {0.1, -0.2, 0.0}, {0.0, 0.3, 0.4}, {-0.2, 0.1, -0.5}};
        const double dtOrig = 2.5e-12;
        const double qOrig  = 1.6e-19;
        createParticlesAt(pc, positions, dtOrig);
        pc->setQ(qOrig);

        // Build dedicated rho fields for baseline and image-charge deposition.
        ippl::Vector<int, 3> nr        = 8;
        ippl::Vector<double, 3> rmin   = -4.0;
        ippl::Vector<double, 3> rmax   = 4.0;
        ippl::Vector<double, 3> origin = rmin;
        ippl::Vector<double, 3> hr     = (rmax - rmin) / ippl::Vector<double, 3>(nr);
        std::array<bool, 3> decomp     = {true, true, true};
        ippl::NDIndex<3> domain;
        for (unsigned i = 0; i < 3; ++i) {
            domain[i] = ippl::Index(nr[i]);
        }
        Mesh_t<3> mesh(domain, hr, origin);
        FieldLayout_t<3> fl(MPI_COMM_WORLD, domain, decomp, true);
        Field_t<3> rhoBaseline;
        Field_t<3> rhoImage;
        rhoBaseline.initialize(mesh, fl);
        rhoImage.initialize(mesh, fl);
        rhoBaseline = 0.0;
        rhoImage    = 0.0;

        constexpr double zPlane = 0.125;
        ImageChargeScatterController<double, 3> baselineController(false, zPlane);
        ImageChargeScatterController<double, 3> imageController(true, zPlane);

        auto R_before_view =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
        auto dt_before_view =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->dt.getView());
        auto q_before_view =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->getQView());
        std::vector<std::array<double, 3>> R_before(pc->getLocalNum());
        std::vector<double> dt_before(pc->getLocalNum(), 0.0);
        for (size_t i = 0; i < pc->getLocalNum(); ++i) {
            R_before[i][0] = R_before_view(i)[0];
            R_before[i][1] = R_before_view(i)[1];
            R_before[i][2] = R_before_view(i)[2];
            dt_before[i]   = dt_before_view(i);
        }
        const double q_before = q_before_view(0);

        // Run the exact production path: primary-only and primary+image scatter.
        baselineController.scatterPrimaryAndImage(pc, pc->R, rhoBaseline);
        imageController.scatterPrimaryAndImage(pc, pc->R, rhoImage);
        Kokkos::fence();

        // Particle state must be restored after image scatter.
        auto R_after  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
        auto dt_after = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->dt.getView());
        auto q_after  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->getQView());
        for (size_t i = 0; i < pc->getLocalNum(); ++i) {
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_NEAR(R_after(i)[d], R_before[i][d], 1e-14);
            }
            EXPECT_NEAR(dt_after(i), dt_before[i], 1e-20);
        }
        EXPECT_NEAR(q_after(0), q_before, 1e-30);

        // Verify image scatter actually changes deposited rho versus baseline.
        auto rhoBaseHost =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rhoBaseline.getView());
        auto rhoImgHost =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rhoImage.getView());
        bool rhoDiffers = false;
        for (size_t i = 0; i < rhoBaseHost.extent(0) && !rhoDiffers; ++i) {
            for (size_t j = 0; j < rhoBaseHost.extent(1) && !rhoDiffers; ++j) {
                for (size_t k = 0; k < rhoBaseHost.extent(2); ++k) {
                    if (std::abs(rhoBaseHost(i, j, k) - rhoImgHost(i, j, k)) > 0.0) {
                        rhoDiffers = true;
                        break;
                    }
                }
            }
        }
        EXPECT_TRUE(rhoDiffers);
    }

    // ================================================================
    // allocateParticles — pre-allocates capacity, throws on non-empty container
    // ================================================================

    TEST_F(ParticleContainerTest, AllocateParticles_ReservesCapacityAndKeepsLocalNumZero) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();
        ASSERT_EQ(pc->R.size(), 0u);

        constexpr size_t nReserve = 1024;
        pc->allocateParticles(nReserve);

        // Logical particle count stays at 0 — alloc only reserves underlying capacity.
        EXPECT_EQ(pc->getLocalNum(), 0u);
        // Underlying view capacity is at least the requested size (overalloc may give more).
        EXPECT_GE(pc->R.size(), nReserve);
        EXPECT_GE(pc->P.size(), nReserve);
        EXPECT_GE(pc->dt.size(), nReserve);

        // Calling allocateParticles on a non-empty container must throw — alloc is destructive.
        EXPECT_THROW(pc->allocateParticles(nReserve), OpalException);
    }

    // ================================================================
    // createParticles — non-destructive grow preserves existing data
    // ================================================================

    TEST_F(ParticleContainerTest, CreateParticles_GrowPreservesExistingData) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();

        // Step 1: pre-allocate a small capacity, fill with a known pattern.
        constexpr size_t nFirst = 8;
        pc->allocateParticles(nFirst);
        pc->createParticles(nFirst);
        ASSERT_EQ(pc->getLocalNum(), nFirst);

        auto P_host = pc->P.getHostMirror();
        for (size_t i = 0; i < nFirst; ++i) {
            P_host(i) = ippl::Vector<double, 3>(double(i), double(2 * i), double(3 * i));
        }
        Kokkos::deep_copy(pc->P.getView(), P_host);
        Kokkos::fence();

        // Step 2: force a capacity grow by adding more particles than the buffer holds.
        const size_t capBefore = pc->R.size();
        const size_t nSecond   = capBefore + 16;
        pc->createParticles(nSecond);

        EXPECT_EQ(pc->getLocalNum(), nFirst + nSecond);
        EXPECT_GE(pc->R.size(), nFirst + nSecond);

        // Step 3: verify the original pattern in the first nFirst slots survived the grow.
        auto P_after = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
        for (size_t i = 0; i < nFirst; ++i) {
            EXPECT_DOUBLE_EQ(P_after(i)[0], double(i)) << "particle " << i;
            EXPECT_DOUBLE_EQ(P_after(i)[1], double(2 * i)) << "particle " << i;
            EXPECT_DOUBLE_EQ(P_after(i)[2], double(3 * i)) << "particle " << i;
        }
    }

    // ================================================================
    // destroyParticles — removes marked entries and validates inputs
    // ================================================================

    TEST_F(ParticleContainerTest, DestroyParticles_RemovesMarkedAndThrowsOnInvalidInput) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();

        constexpr size_t nPart = 10;
        std::vector<std::array<double, 3>> positions(nPart);
        for (size_t i = 0; i < nPart; ++i) {
            positions[i] = {double(i) * 0.1, 0.0, 0.0};
        }
        createParticlesAt(pc, positions);
        ASSERT_EQ(pc->getLocalNum(), nPart);

        // Mark every other particle invalid (5 of 10) — kill those at even indices.
        Kokkos::View<bool*> invalid("DestroyParticlesTest::invalid", nPart);
        auto invalid_host         = Kokkos::create_mirror_view(invalid);
        size_type expectedDestroy = 0;
        for (size_t i = 0; i < nPart; ++i) {
            const bool kill = (i % 2 == 0);
            invalid_host(i) = kill;
            if (kill) ++expectedDestroy;
        }
        Kokkos::deep_copy(invalid, invalid_host);

        pc->destroyParticles(invalid, expectedDestroy);
        EXPECT_EQ(pc->getLocalNum(), nPart - expectedDestroy);

        // localDestroyNum exceeding the local particle count must throw.
        Kokkos::View<bool*> oversize_invalid(
                "DestroyParticlesTest::oversize_invalid", pc->getLocalNum());
        EXPECT_THROW(
                pc->destroyParticles(oversize_invalid, pc->getLocalNum() + 1u), OpalException);

        // invalid mask smaller than the local particle count must throw.
        ASSERT_GT(pc->getLocalNum(), 0u);
        Kokkos::View<bool*> short_invalid(
                "DestroyParticlesTest::short_invalid", pc->getLocalNum() - 1);
        EXPECT_THROW(pc->destroyParticles(short_invalid, 0u), OpalException);
    }

}  // namespace
