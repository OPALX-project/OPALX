/**
 * @file TestParticleContainer.cpp
 * @brief Unit tests for `ParticleContainer<double, 3>`.
 *
 * Covers:
 * - Charge / mass storage in both `SingleValue` and `Attributes` QM modes.
 * - Round-trip `scaleDtByCharge` / `unscaleDtByCharge` in both QM modes.
 * - Moment computation stability (empty and known-data cases).
 * - Min / max position computation.
 * - `markParticlesOutside` and `deleteInvalidParticles` mask-driven deletion.
 * - Propagation of periodic and open boundary conditions to the P3M overlap layout.
 * - Rotation of `Pol` together with `P` in every frame change (`transformBunch`,
 *   `updateRefToLabCSTrafo`, space-charge solve frame round trip).
 *
 * The fixture constructs a lightweight `ParticleContainer` on a minimal 8^3
 * periodic mesh without a full `PartBunch`, so no field-solver or DataSink
 * infrastructure is required.
 */

#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

#include "Ippl.h"
#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/SpaceChargeFrames.h"
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

        /// Build a minimal ParticleContainer on an 8^3 mesh over [-4, 4]^3 with selectable layout,
        /// particle boundary conditions and spin storage.
        std::shared_ptr<PC_t> makeContainer(
                PC_t::LayoutType layoutType = PC_t::LayoutType::Spatial,
                ippl::BC particleBC = ippl::BC::PERIODIC, bool spinEnabled = false) {
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

            constexpr double overlapCutoff = 1.0;
            std::shared_ptr<PC_t> pc       = std::make_shared<PC_t>(
                    mesh, fl, spinEnabled, layoutType, overlapCutoff, particleBC);
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

    TEST_F(ParticleContainerTest, P3MLayoutUsesSelectedParticleBoundaryCondition) {
        auto periodic = makeContainer(PC_t::LayoutType::SpatialOverlap, ippl::BC::PERIODIC);
        auto open     = makeContainer(PC_t::LayoutType::SpatialOverlap, ippl::BC::NO);

        ASSERT_TRUE(periodic->hasP3MLayout());
        ASSERT_TRUE(open->hasP3MLayout());

        for (const auto bc : periodic->getP3MLayout().getParticleBC()) {
            EXPECT_EQ(bc, ippl::BC::PERIODIC);
        }
        for (const auto bc : open->getP3MLayout().getParticleBC()) {
            EXPECT_EQ(bc, ippl::BC::NO);
        }
    }

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
    // markParticlesOutside – edge cases
    // ================================================================

    TEST_F(ParticleContainerTest, MarkParticlesOutside_NegativeSigma_ReturnsZero) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();

        std::vector<std::array<double, 3>> positions = {{{0.0, 0.0, 0.0}}};
        createParticlesAt(pc, positions);
        pc->setM(1.0);

        size_t localBefore = pc->getLocalNum();
        auto marked        = pc->markParticlesOutside(-1.0);

        EXPECT_EQ(marked, 0u);
        EXPECT_EQ(pc->getLocalNum(), localBefore);
    }

    TEST_F(ParticleContainerTest, MarkParticlesOutside_ZeroSigma_ReturnsZero) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();

        std::vector<std::array<double, 3>> positions = {{{0.0, 0.0, 0.0}}};
        createParticlesAt(pc, positions);
        pc->setM(1.0);

        size_t localBefore = pc->getLocalNum();
        auto marked        = pc->markParticlesOutside(0.0);

        EXPECT_EQ(marked, 0u);
        EXPECT_EQ(pc->getLocalNum(), localBefore);
    }

    TEST_F(ParticleContainerTest, MarkParticlesOutside_EmptyContainer_ReturnsZero) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();

        ASSERT_EQ(pc->getLocalNum(), 0u);
        auto marked = pc->markParticlesOutside(3.0);
        EXPECT_EQ(marked, 0u);
    }

    // ================================================================
    // markParticlesOutside + deleteInvalidParticles – known outlier removal
    // ================================================================

    TEST_F(ParticleContainerTest, MarkAndDeleteInvalidParticles_KnownOutlierRemoval) {
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
        auto marked        = pc->markParticlesOutside(2.0);

        EXPECT_EQ(marked, 1u);

        auto destroyed = pc->deleteInvalidParticles();
        EXPECT_EQ(destroyed, marked);

        size_t totalAfter = pc->getTotalNum();
        EXPECT_EQ(totalAfter + destroyed, totalBefore);
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
    // deleteInvalidParticles — removes marked InvalidMask entries and resets the mask
    // ================================================================

    TEST_F(ParticleContainerTest, DeleteInvalidParticles_RemovesMarkedAndResetsMask) {
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
        auto invalid_host         = pc->InvalidMask.getHostMirror();
        size_type expectedDestroy = 0;
        for (size_t i = 0; i < nPart; ++i) {
            const bool kill = (i % 2 == 0);
            invalid_host(i) = kill;
            if (kill) ++expectedDestroy;
        }
        Kokkos::deep_copy(pc->InvalidMask.getView(), invalid_host);

        const size_type expectedGlobalDestroy = pc->getTotalNum() / 2;
        const size_type destroyed             = pc->deleteInvalidParticles();
        EXPECT_EQ(destroyed, expectedGlobalDestroy);
        EXPECT_EQ(pc->getLocalNum(), nPart - expectedDestroy);

        auto resetHost =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->InvalidMask.getView());
        for (size_t i = 0; i < pc->getLocalNum(); ++i) {
            EXPECT_FALSE(resetHost(i));
        }
    }

    TEST_F(ParticleContainerTest, CreateParticles_ResetsInvalidMaskForNewParticles) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();

        constexpr size_t nPart = 6;
        pc->createParticles(nPart);
        const size_t totalBeforeDelete = pc->getTotalNum();
        pc->InvalidMask                = true;
        Kokkos::fence();
        ASSERT_EQ(pc->deleteInvalidParticles(), totalBeforeDelete);
        ASSERT_EQ(pc->getLocalNum(), 0u);

        pc->createParticles(nPart);
        auto invalidHost =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->InvalidMask.getView());
        for (size_t i = 0; i < pc->getLocalNum(); ++i) {
            EXPECT_FALSE(invalidHost(i));
        }
    }

    // ================================================================
    // ID assignment — IDs follow the IPPL strided rank-interleave scheme
    // ================================================================

    TEST_F(ParticleContainerTest, CreateParticles_AssignsStridedIDs) {
        Options::useQMAttributes = false;
        auto pc                  = makeContainer();

        constexpr size_t nPart = 32;
        pc->createParticles(nPart);
        ASSERT_EQ(pc->getLocalNum(), nPart);

        auto idHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->ID.getView());

        const auto rank   = static_cast<std::int64_t>(ippl::Comm->rank());
        const auto nRanks = static_cast<std::int64_t>(ippl::Comm->size());

        for (size_t i = 0; i < nPart; ++i) {
            EXPECT_EQ(idHost(i), rank + nRanks * static_cast<std::int64_t>(i)) << "particle " << i;
        }
    }

    // ================================================================
    // Pol follows every frame change that P follows
    // ================================================================

    /// Spin-enabled container with one particle at R, with momentum P and spin Pol.
    std::shared_ptr<PC_t> makeSpinParticle(
            std::shared_ptr<PC_t> pc, const Vector_t<double, 3>& R, const Vector_t<double, 3>& P,
            const Vector_t<double, 3>& Pol) {
        pc->createParticles(1);
        auto rHost   = pc->R.getHostMirror();
        auto pHost   = pc->P.getHostMirror();
        auto polHost = pc->Pol.getHostMirror();
        rHost(0)     = R;
        pHost(0)     = P;
        polHost(0)   = Pol;
        Kokkos::deep_copy(pc->R.getView(), rHost);
        Kokkos::deep_copy(pc->P.getView(), pHost);
        Kokkos::deep_copy(pc->Pol.getView(), polHost);
        return pc;
    }

    Vector_t<double, 3> firstOf(const ippl::ParticleAttrib<Vector_t<double, 3>>& attrib) {
        auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), attrib.getView());
        return host(0);
    }

    TEST_F(ParticleContainerTest, TransformBunch_RotatesPolLikeP) {
        Options::useQMAttributes = false;
        const Vector_t<double, 3> P0(0.3, 0.1, 0.5), Pol0(0.6, 0.0, 0.8);
        auto pc = makeSpinParticle(
                makeContainer(PC_t::LayoutType::Spatial, ippl::BC::PERIODIC, true),
                Vector_t<double, 3>(0.1, 0.2, 0.3), P0, Pol0);

        // 45 degree turn about y plus a shift.
        const CoordinateSystemTrafo trafo(
                Vector_t<double, 3>(0.05, 0.0, 0.1),
                getQuaternion(Vector_t<double, 3>(1, 0, 1), Vector_t<double, 3>(0, 0, 1)));
        pc->transformBunch(trafo);

        const Vector_t<double, 3> P = firstOf(pc->P), Pol = firstOf(pc->Pol);
        const Vector_t<double, 3> expected = trafo.rotateTo(Pol0);
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_NEAR(Pol[d], expected[d], 1e-15) << "component " << d;
        }
        EXPECT_NEAR(Pol.dot(P), Pol0.dot(P0), 1e-15);
        EXPECT_GT(std::abs(Pol[0] - Pol0[0]), 0.1);  // the rotation is not the identity
    }

    TEST_F(ParticleContainerTest, UpdateRefToLabCSTrafo_KeepsPolAlongP) {
        Options::useQMAttributes = false;
        // Reference particle in the lab frame; the bunch frame starts as the lab frame.
        const Vector_t<double, 3> refR(0.01, 0.0, 0.02), refP(0.1, 0.05, 0.4);
        const double pAbs = std::sqrt(refP.dot(refP));
        auto pc           = makeSpinParticle(
                makeContainer(PC_t::LayoutType::Spatial, ippl::BC::PERIODIC, true), refR, refP,
                refP / pAbs);
        pc->setRefPartR(refR);
        pc->setRefPartP(refP);

        pc->updateRefToLabCSTrafo(1e-12);

        // The frame now sits on the reference particle with z along its momentum, so a particle
        // on the reference with its spin along P has both P and Pol along z.
        const Vector_t<double, 3> P = firstOf(pc->P), Pol = firstOf(pc->Pol);
        EXPECT_NEAR(P[0], 0.0, 1e-14);
        EXPECT_NEAR(P[1], 0.0, 1e-14);
        EXPECT_NEAR(P[2], pAbs, 1e-14);
        EXPECT_NEAR(Pol[0], 0.0, 1e-14);
        EXPECT_NEAR(Pol[1], 0.0, 1e-14);
        EXPECT_NEAR(Pol[2], 1.0, 1e-14);
    }

    TEST_F(ParticleContainerTest, FollowFrameOf_KeepsLabPositionsAndCopiesRoot) {
        Options::useQMAttributes = false;
        // Root: one particle on its reference, which is off the origin and turned.
        const Vector_t<double, 3> refR(0.01, 0.0, 0.02), refP(0.1, 0.05, 0.4);
        auto root = makeSpinParticle(
                makeContainer(PC_t::LayoutType::Spatial, ippl::BC::PERIODIC, true), refR, refP,
                Vector_t<double, 3>(0, 0, 1));
        root->setRefPartR(refR);
        root->setRefPartP(refP);
        root->set_sPos(1.0);

        // Follower: a decay daughter somewhere else, in the same (lab) frame as the root.
        const Vector_t<double, 3> R0(0.03, -0.02, 0.05), P0(-0.02, 0.1, 0.3);
        auto follower = makeSpinParticle(
                makeContainer(PC_t::LayoutType::Spatial, ippl::BC::PERIODIC, true), R0, P0,
                Vector_t<double, 3>(0, 0, 1));

        follower->followFrameOf(*root, root->updateRefToLabCSTrafo(1e-12));

        const CoordinateSystemTrafo& toLab = follower->getToLabTrafo();
        const Vector_t<double, 3> labR     = toLab.transformTo(firstOf(follower->R));
        const Vector_t<double, 3> labP     = toLab.rotateTo(firstOf(follower->P));
        const Vector_t<double, 3> probe(0.3, -0.7, 1.1);
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_NEAR(labR[d], R0[d], 1e-15) << "component " << d;
            EXPECT_NEAR(labP[d], P0[d], 1e-15) << "component " << d;
            EXPECT_EQ(toLab.transformTo(probe)[d], root->getToLabTrafo().transformTo(probe)[d]);
            EXPECT_EQ(follower->getRefPartR()[d], root->getRefPartR()[d]);
            EXPECT_EQ(follower->getRefPartP()[d], root->getRefPartP()[d]);
        }
        EXPECT_EQ(follower->get_sPos(), root->get_sPos());
        EXPECT_GT(std::abs(firstOf(follower->R)[0] - R0[0]), 1e-3);  // the frame did change
    }

    TEST_F(ParticleContainerTest, EnterLeaveSolveFrame_KeepsPol) {
        Options::useQMAttributes = false;
        const Vector_t<double, 3> Pol0(0.0, 0.6, 0.8);
        auto pc = makeSpinParticle(
                makeContainer(PC_t::LayoutType::Spatial, ippl::BC::PERIODIC, true),
                Vector_t<double, 3>(0.1, 0.2, 0.3), Vector_t<double, 3>(0.3, 0.1, 0.5), Pol0);

        const CoordinateSystemTrafo solveToTracker(
                Vector_t<double, 3>(0.0, 0.0, 1.5),
                getQuaternion(Vector_t<double, 3>(1, 1, 1), Vector_t<double, 3>(0, 0, 1)));
        const opalx::spacecharge::CoordinateFrameTransforms frames{
                solveToTracker.inverted(), solveToTracker};

        opalx::spacecharge::enterSolveFrame(frames, *pc);
        EXPECT_GT(std::abs(firstOf(pc->Pol)[1] - Pol0[1]), 0.01);  // Pol is in the solve axes
        opalx::spacecharge::leaveSolveFrame(frames, *pc);

        const Vector_t<double, 3> Pol = firstOf(pc->Pol);
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_NEAR(Pol[d], Pol0[d], 1e-15) << "component " << d;
        }
    }

}  // namespace
