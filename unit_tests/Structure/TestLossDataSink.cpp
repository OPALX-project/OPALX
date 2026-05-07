 /**
 * \file TestLossDataSink.cpp
 * \brief Unit tests for LossDataSink statistics and bookkeeping.
 *
 * This test suite validates the core behavior of LossDataSink using:
 *  - small deterministic OpalParticle samples
 *  - direct access to internal bookkeeping where needed
 *  - a single-rank initialized IPPL/MPI environment
 *
 * ---------------------------------------------------------------------------
 * Coverage
 * ---------------------------------------------------------------------------
 *
 * 1. Particle bookkeeping
 *    - plain particles and turn/bunch particles cannot be mixed
 *    - inconsistent turn/bunch information is rejected
 *
 * 2. Basic statistics
 *    - total particle count
 *    - centroids in position and momentum
 *    - RMS position and momentum
 *    - total charge and total mass
 *
 * 3. Time statistics
 *    - mean particle time
 *    - RMS particle time
 *
 * 4. Spatial extent statistics
 *    - rmin / rmax are computed correctly
 *    - maxR uses the maximum absolute position extent
 *
 * 5. Energy statistics
 *    - mean kinetic energy is computed from particle momentum and mass
 *    - kinetic-energy RMS is stable for identical particle energies
 *
 * 6. Set splitting
 *    - splitSets() clears stale state
 *    - split boundaries are monotonic
 *    - split boundaries cover all local particles
 *
 * ---------------------------------------------------------------------------
 * Notes
 * ---------------------------------------------------------------------------
 *
 * - Options::computePercentiles is disabled in the fixture so these tests focus
 *   on the core LossDataSink statistics and do not depend on percentile
 *   histogram behavior.
 *
 * - The tests initialize IPPL in SetUpTestSuite(), because LossDataSink uses
 *   ippl::Comm reductions and barriers even in single-rank unit tests.
 *
 * - The private/public macro is used only to test internal set-splitting state
 *   and other bookkeeping that is not exposed through the public interface.
 */

#include <gtest/gtest.h>

#include "Ippl.h"

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#define private public
#include "Structure/LossDataSink.h"
#undef private

#include "Utilities/GeneralOpalException.h"
#include "Utilities/Options.h"

namespace {

    Vector_t<double, 3> makeVector(double x, double y, double z) {
        Vector_t<double, 3> v(0.0);
        v(0) = x;
        v(1) = y;
        v(2) = z;
        return v;
    }

    OpalParticle makeParticle(
            std::size_t id,
            const Vector_t<double, 3>& r,
            const Vector_t<double, 3>& p,
            double time = 0.0,
            double charge = -1.0e-17,
            double mass = 0.51099895) {
        return OpalParticle(id, r, p, time, charge, mass);
    }

    class LossDataSinkTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }

        static void TearDownTestSuite() {
            ippl::finalize();
        }

        void SetUp() override {
            oldComputePercentiles_m = Options::computePercentiles;
            Options::computePercentiles = false;
        }

        void TearDown() override {
            Options::computePercentiles = oldComputePercentiles_m;
        }

    private:
        bool oldComputePercentiles_m = false;
    };

    TEST_F(LossDataSinkTest, AddParticleRejectsTurnInformationAfterPlainParticle) {
        LossDataSink sink;

        sink.addParticle(makeParticle(
            1,
            makeVector(0.0, 0.0, 0.0),
            makeVector(0.0, 0.0, 1.0)));

        EXPECT_THROW(
            sink.addParticle(
                makeParticle(
                    2,
                    makeVector(0.0, 0.0, 0.0),
                    makeVector(0.0, 0.0, 1.0)),
                std::make_pair(3, static_cast<short>(4))),
            GeneralOpalException);
    }

    TEST_F(LossDataSinkTest, AddParticleRejectsPlainParticleAfterTurnInformation) {
        LossDataSink sink;

        sink.addParticle(
            makeParticle(
                1,
                makeVector(0.0, 0.0, 0.0),
                makeVector(0.0, 0.0, 1.0)),
            std::make_pair(3, static_cast<short>(4)));

        EXPECT_THROW(
            sink.addParticle(makeParticle(
                2,
                makeVector(0.0, 0.0, 0.0),
                makeVector(0.0, 0.0, 1.0))),
            GeneralOpalException);
    }

    TEST_F(LossDataSinkTest, ComputeStatisticsSingleParticle) {
        LossDataSink sink;

        const Vector_t<double, 3> r = makeVector(1.0, -2.0, 3.0);
        const Vector_t<double, 3> p = makeVector(0.1, -0.2, 4.0);

        sink.addParticle(makeParticle(1, r, p, 5.0, -2.0e-17, 0.51099895));

        const auto stats = sink.computeStatistics(1);

        ASSERT_EQ(stats.size(), 1u);

        const SetStatistics& stat = *stats.begin();

        EXPECT_EQ(stat.nTotal_m, 1u);

        EXPECT_DOUBLE_EQ(stat.rmean_m(0), 1.0);
        EXPECT_DOUBLE_EQ(stat.rmean_m(1), -2.0);
        EXPECT_DOUBLE_EQ(stat.rmean_m(2), 3.0);

        EXPECT_DOUBLE_EQ(stat.pmean_m(0), 0.1);
        EXPECT_DOUBLE_EQ(stat.pmean_m(1), -0.2);
        EXPECT_DOUBLE_EQ(stat.pmean_m(2), 4.0);

        EXPECT_DOUBLE_EQ(stat.rrms_m(0), 0.0);
        EXPECT_DOUBLE_EQ(stat.rrms_m(1), 0.0);
        EXPECT_DOUBLE_EQ(stat.rrms_m(2), 0.0);

        EXPECT_DOUBLE_EQ(stat.prms_m(0), 0.0);
        EXPECT_DOUBLE_EQ(stat.prms_m(1), 0.0);
        EXPECT_DOUBLE_EQ(stat.prms_m(2), 0.0);

        EXPECT_DOUBLE_EQ(stat.tmean_m, 5.0);
        EXPECT_DOUBLE_EQ(stat.trms_m, 0.0);

        EXPECT_DOUBLE_EQ(stat.totalCharge_m, -2.0e-17);
        EXPECT_DOUBLE_EQ(stat.totalMass_m, 0.51099895);
    }

    TEST_F(LossDataSinkTest, MaxRUsesMaximumAbsoluteExtent) {
        LossDataSink sink;

        sink.addParticle(makeParticle(
            1,
            makeVector(-2.0, 5.0, -7.0),
            makeVector(0.0, 0.0, 1.0)));

        sink.addParticle(makeParticle(
            2,
            makeVector(1.0, -6.0, 4.0),
            makeVector(0.0, 0.0, 1.0)));

        const auto stats = sink.computeStatistics(1);

        ASSERT_EQ(stats.size(), 1u);

        const SetStatistics& stat = *stats.begin();

        EXPECT_DOUBLE_EQ(stat.rmin_m(0), -2.0);
        EXPECT_DOUBLE_EQ(stat.rmax_m(0), 1.0);
        EXPECT_DOUBLE_EQ(stat.maxR_m(0), 2.0);

        EXPECT_DOUBLE_EQ(stat.rmin_m(1), -6.0);
        EXPECT_DOUBLE_EQ(stat.rmax_m(1), 5.0);
        EXPECT_DOUBLE_EQ(stat.maxR_m(1), 6.0);

        EXPECT_DOUBLE_EQ(stat.rmin_m(2), -7.0);
        EXPECT_DOUBLE_EQ(stat.rmax_m(2), 4.0);
        EXPECT_DOUBLE_EQ(stat.maxR_m(2), 7.0);
    }

    TEST_F(LossDataSinkTest, MeanAndRmsTimeAreComputedCorrectly) {
        LossDataSink sink;

        sink.addParticle(makeParticle(
            1,
            makeVector(0.0, 0.0, 0.0),
            makeVector(0.0, 0.0, 1.0),
            1.0));

        sink.addParticle(makeParticle(
            2,
            makeVector(0.0, 0.0, 0.0),
            makeVector(0.0, 0.0, 1.0),
            3.0));

        const auto stats = sink.computeStatistics(1);

        ASSERT_EQ(stats.size(), 1u);

        const SetStatistics& stat = *stats.begin();

        EXPECT_DOUBLE_EQ(stat.tmean_m, 2.0);
        EXPECT_DOUBLE_EQ(stat.trms_m, 1.0);
    }

    TEST_F(LossDataSinkTest, StdKineticEnergyIsStableForIdenticalParticles) {
        LossDataSink sink;

        const Vector_t<double, 3> r = makeVector(0.0, 0.0, 0.0);
        const Vector_t<double, 3> p = makeVector(0.0, 0.0, 13.87990113);

        constexpr std::size_t numParticles = 10000;

        for (std::size_t i = 0; i < numParticles; ++i) {
            sink.addParticle(makeParticle(i, r, p, 1.0, -1.5e-17, 0.51099895));
        }

        const auto stats = sink.computeStatistics(1);

        ASSERT_EQ(stats.size(), 1u);

        const SetStatistics& stat = *stats.begin();

        EXPECT_EQ(stat.nTotal_m, numParticles);
        EXPECT_NEAR(stat.meanKineticEnergy_m, 6.6, 1.0e-6);
        EXPECT_NEAR(stat.stdKineticEnergy_m, 0.0, 1.0e-12);
    }

    TEST_F(LossDataSinkTest, SplitSetsClearsPreviousStateForSingleSet) {
        LossDataSink sink;

        sink.startSet_m = {0, 1, 2};

        sink.addParticle(makeParticle(
            1,
            makeVector(0.0, 0.0, 0.0),
            makeVector(0.0, 0.0, 1.0),
            1.0));

        sink.splitSets(1);

        EXPECT_TRUE(sink.startSet_m.empty());
    }

    TEST_F(LossDataSinkTest, SplitSetsCreatesMonotonicBoundaries) {
        LossDataSink sink;

        sink.addParticle(makeParticle(
            1,
            makeVector(0.0, 0.0, 0.0),
            makeVector(0.0, 0.0, 1.0),
            0.0));

        sink.addParticle(makeParticle(
            2,
            makeVector(0.0, 0.0, 0.0),
            makeVector(0.0, 0.0, 1.0),
            1.0));

        sink.addParticle(makeParticle(
            3,
            makeVector(0.0, 0.0, 0.0),
            makeVector(0.0, 0.0, 1.0),
            10.0));

        sink.addParticle(makeParticle(
            4,
            makeVector(0.0, 0.0, 0.0),
            makeVector(0.0, 0.0, 1.0),
            11.0));

        sink.splitSets(2);

        ASSERT_EQ(sink.startSet_m.size(), 3u);

        EXPECT_EQ(sink.startSet_m.front(), 0u);
        EXPECT_EQ(sink.startSet_m.back(), sink.particles_m.size());

        EXPECT_LE(sink.startSet_m[0], sink.startSet_m[1]);
        EXPECT_LE(sink.startSet_m[1], sink.startSet_m[2]);

        EXPECT_EQ(sink.startSet_m[1], 2u);
    }

    TEST_F(LossDataSinkTest, KineticEnergyUsesParticleMass) {
        LossDataSink sink;

        const auto r = makeVector(0.0, 0.0, 0.0);
        const auto p = makeVector(0.0, 0.0, 13.87990113);

        sink.addParticle(makeParticle(1, r, p, 1.0, -1.0e-17, 0.51099895));

        const auto stats = sink.computeStatistics(1);
        ASSERT_EQ(stats.size(), 1u);

        const SetStatistics& stat = *stats.begin();

        EXPECT_NEAR(stat.meanKineticEnergy_m, 6.6, 1.0e-6);
        EXPECT_NEAR(stat.totalMass_m, 0.51099895, 1.0e-12);
    }

}  // namespace