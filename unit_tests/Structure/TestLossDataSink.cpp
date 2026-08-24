/**
 * \file TestLossDataSink.cpp
 * \brief Unit tests for LossDataSink statistics and bookkeeping.
 *
 * This test suite validates the core behavior of LossDataSink using:
 *  - small deterministic OpalParticle samples
 *  - direct access to internal bookkeeping where needed
 *  - an initialized IPPL/MPI environment
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
 * 7. Monitor HDF5 checkpoint rewind
 *    - records at and after the checkpoint global step are removed
 *    - earlier records are preserved
 *    - repeated rewind and append does not accumulate duplicate records
 *    - malformed records without GlobalTrackStep are rejected
 *
 * 8. HDF5 file-space management
 *    - rewind physically compacts discarded steps
 *    - repeated rewind and append cycles do not grow the file
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
#include <filesystem>
#include <stdexcept>
#include <string>
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
            std::size_t id, const Vector_t<double, 3>& r, const Vector_t<double, 3>& p,
            double time = 0.0, double charge = -1.0e-17, double mass = 0.51099895) {
        return OpalParticle(id, r, p, time, charge, mass);
    }

    h5_file_t openCollectiveH5(const std::string& fileName, h5_int32_t mode) {
        h5_prop_t props = H5CreateFileProp();
        if (props == H5_ERR) {
            throw std::runtime_error("Could not create HDF5 file properties");
        }

        MPI_Comm comm = ippl::Comm->getCommunicator();
        if (H5SetPropFileMPIOCollective(props, &comm) == H5_ERR) {
            H5CloseProp(props);
            throw std::runtime_error("Could not configure collective HDF5 access");
        }

        h5_file_t file = H5OpenFile(fileName.c_str(), mode, props);
        H5CloseProp(props);
        if (file == static_cast<h5_file_t>(H5_ERR)) {
            throw std::runtime_error("Could not open HDF5 file " + fileName);
        }
        return file;
    }

    void writeGlobalTrackSteps(
            const std::string& fileName, const std::vector<h5_int64_t>& globalTrackSteps,
            std::size_t payloadSize = 0) {
        h5_file_t file                = openCollectiveH5(fileName, H5_O_WRONLY);
        const h5_int64_t rootMetadata = 17;
        if (H5WriteFileAttribInt64(file, "TestRootMetadata", &rootMetadata, 1) != H5_SUCCESS) {
            H5CloseFile(file);
            throw std::runtime_error("Could not write HDF5 root metadata test data");
        }
        for (std::size_t step = 0; step < globalTrackSteps.size(); ++step) {
            if (H5SetStep(file, static_cast<h5_int64_t>(step)) != H5_SUCCESS
                || H5WriteStepAttribInt64(file, "GlobalTrackStep", &globalTrackSteps[step], 1)
                           != H5_SUCCESS) {
                H5CloseFile(file);
                throw std::runtime_error("Could not write GlobalTrackStep test data");
            }
            if (payloadSize > 0) {
                std::vector<h5_float64_t> payload(
                        payloadSize, static_cast<h5_float64_t>(globalTrackSteps[step]));
                if (H5PartSetNumParticles(file, static_cast<h5_size_t>(payload.size()))
                            != H5_SUCCESS
                    || H5PartWriteDataFloat64(file, "x", payload.data()) != H5_SUCCESS) {
                    H5CloseFile(file);
                    throw std::runtime_error("Could not write HDF5 payload test data");
                }
            }
        }
        if (H5CloseFile(file) != H5_SUCCESS) {
            throw std::runtime_error("Could not close HDF5 test file after writing");
        }
    }

    void appendGlobalTrackStep(
            const std::string& fileName, h5_int64_t globalTrackStep, std::size_t payloadSize = 0) {
        h5_file_t file            = openCollectiveH5(fileName, H5_O_APPENDONLY);
        const h5_ssize_t numSteps = H5GetNumSteps(file);
        if (numSteps < 0 || H5SetStep(file, numSteps) != H5_SUCCESS
            || H5WriteStepAttribInt64(file, "GlobalTrackStep", &globalTrackStep, 1) != H5_SUCCESS) {
            H5CloseFile(file);
            throw std::runtime_error("Could not append GlobalTrackStep test data");
        }
        if (payloadSize > 0) {
            std::vector<h5_float64_t> payload(
                    payloadSize, static_cast<h5_float64_t>(globalTrackStep));
            if (H5PartSetNumParticles(file, static_cast<h5_size_t>(payload.size())) != H5_SUCCESS
                || H5PartWriteDataFloat64(file, "x", payload.data()) != H5_SUCCESS) {
                H5CloseFile(file);
                throw std::runtime_error("Could not append HDF5 payload test data");
            }
        }
        if (H5CloseFile(file) != H5_SUCCESS) {
            throw std::runtime_error("Could not close HDF5 test file after appending");
        }
    }

    std::vector<h5_int64_t> readGlobalTrackSteps(const std::string& fileName) {
        h5_file_t file            = openCollectiveH5(fileName, H5_O_RDONLY);
        const h5_ssize_t numSteps = H5GetNumSteps(file);
        if (numSteps < 0) {
            H5CloseFile(file);
            throw std::runtime_error("Could not read HDF5 test step count");
        }

        std::vector<h5_int64_t> globalTrackSteps(static_cast<std::size_t>(numSteps));
        for (h5_ssize_t step = 0; step < numSteps; ++step) {
            if (H5SetStep(file, step) != H5_SUCCESS
                || H5ReadStepAttribInt64(
                           file, "GlobalTrackStep",
                           &globalTrackSteps[static_cast<std::size_t>(step)])
                           != H5_SUCCESS) {
                H5CloseFile(file);
                throw std::runtime_error("Could not read GlobalTrackStep test data");
            }
        }
        if (H5CloseFile(file) != H5_SUCCESS) {
            throw std::runtime_error("Could not close HDF5 test file after reading");
        }
        return globalTrackSteps;
    }

    h5_int64_t readRootMetadata(const std::string& fileName) {
        h5_file_t file          = openCollectiveH5(fileName, H5_O_RDONLY);
        h5_int64_t rootMetadata = 0;
        if (H5ReadFileAttribInt64(file, "TestRootMetadata", &rootMetadata) != H5_SUCCESS) {
            H5CloseFile(file);
            throw std::runtime_error("Could not read HDF5 root metadata test data");
        }
        if (H5CloseFile(file) != H5_SUCCESS) {
            throw std::runtime_error("Could not close HDF5 test file after reading metadata");
        }
        return rootMetadata;
    }

    class LossDataSinkTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }

        static void TearDownTestSuite() { ippl::finalize(); }

        void SetUp() override {
            oldComputePercentiles_m     = Options::computePercentiles;
            Options::computePercentiles = false;
        }

        void TearDown() override {
            Options::computePercentiles = oldComputePercentiles_m;

            ippl::Comm->barrier();
            if (ippl::Comm->rank() == 0) {
                std::filesystem::remove("test_monitor_checkpoint_rewind.h5");
                std::filesystem::remove("test_monitor_checkpoint_keep_all.h5");
                std::filesystem::remove("test_monitor_checkpoint_missing_step.h5");
                std::filesystem::remove("test_monitor_checkpoint_compaction.h5");
                std::filesystem::remove("test_monitor_checkpoint_compaction.h5.opalx-rewind.tmp");
            }
            ippl::Comm->barrier();
        }

    private:
        bool oldComputePercentiles_m = false;
    };

    TEST_F(LossDataSinkTest, AddParticleRejectsTurnInformationAfterPlainParticle) {
        LossDataSink sink;

        sink.addParticle(makeParticle(1, makeVector(0.0, 0.0, 0.0), makeVector(0.0, 0.0, 1.0)));

        EXPECT_THROW(
                sink.addParticle(
                        makeParticle(2, makeVector(0.0, 0.0, 0.0), makeVector(0.0, 0.0, 1.0)),
                        std::make_pair(3, static_cast<short>(4))),
                GeneralOpalException);
    }

    TEST_F(LossDataSinkTest, AddParticleRejectsPlainParticleAfterTurnInformation) {
        LossDataSink sink;

        sink.addParticle(
                makeParticle(1, makeVector(0.0, 0.0, 0.0), makeVector(0.0, 0.0, 1.0)),
                std::make_pair(3, static_cast<short>(4)));

        EXPECT_THROW(
                sink.addParticle(
                        makeParticle(2, makeVector(0.0, 0.0, 0.0), makeVector(0.0, 0.0, 1.0))),
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

        sink.addParticle(makeParticle(1, makeVector(-2.0, 5.0, -7.0), makeVector(0.0, 0.0, 1.0)));

        sink.addParticle(makeParticle(2, makeVector(1.0, -6.0, 4.0), makeVector(0.0, 0.0, 1.0)));

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

        sink.addParticle(
                makeParticle(1, makeVector(0.0, 0.0, 0.0), makeVector(0.0, 0.0, 1.0), 1.0));

        sink.addParticle(
                makeParticle(2, makeVector(0.0, 0.0, 0.0), makeVector(0.0, 0.0, 1.0), 3.0));

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

        sink.addParticle(
                makeParticle(1, makeVector(0.0, 0.0, 0.0), makeVector(0.0, 0.0, 1.0), 1.0));

        sink.splitSets(1);

        EXPECT_TRUE(sink.startSet_m.empty());
    }

    TEST_F(LossDataSinkTest, SplitSetsCreatesMonotonicBoundaries) {
        LossDataSink sink;

        sink.addParticle(
                makeParticle(1, makeVector(0.0, 0.0, 0.0), makeVector(0.0, 0.0, 1.0), 0.0));

        sink.addParticle(
                makeParticle(2, makeVector(0.0, 0.0, 0.0), makeVector(0.0, 0.0, 1.0), 1.0));

        sink.addParticle(
                makeParticle(3, makeVector(0.0, 0.0, 0.0), makeVector(0.0, 0.0, 1.0), 10.0));

        sink.addParticle(
                makeParticle(4, makeVector(0.0, 0.0, 0.0), makeVector(0.0, 0.0, 1.0), 11.0));

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

    TEST_F(LossDataSinkTest, RewindH5RemovesCheckpointAndLaterStepsBeforeAppend) {
        const std::string fileName = "test_monitor_checkpoint_rewind.h5";
        writeGlobalTrackSteps(fileName, {4, 8, 12});

        LossDataSink::rewindH5ToGlobalTrackStep(fileName, 8);
        EXPECT_EQ(readGlobalTrackSteps(fileName), std::vector<h5_int64_t>({4}));

        appendGlobalTrackStep(fileName, 8);
        EXPECT_EQ(readGlobalTrackSteps(fileName), std::vector<h5_int64_t>({4, 8}));

        // Restarting from the same checkpoint again removes the abandoned replacement
        // before appending it once more, so repeated restart cannot accumulate duplicates.
        LossDataSink::rewindH5ToGlobalTrackStep(fileName, 8);
        appendGlobalTrackStep(fileName, 8);
        EXPECT_EQ(readGlobalTrackSteps(fileName), std::vector<h5_int64_t>({4, 8}));
    }

    TEST_F(LossDataSinkTest, RewindH5KeepsEveryStepBeforeCheckpoint) {
        const std::string fileName = "test_monitor_checkpoint_keep_all.h5";
        writeGlobalTrackSteps(fileName, {2, 5});

        LossDataSink::rewindH5ToGlobalTrackStep(fileName, 8);

        EXPECT_EQ(readGlobalTrackSteps(fileName), std::vector<h5_int64_t>({2, 5}));
    }

    TEST_F(LossDataSinkTest, RewindH5RejectsStepWithoutGlobalTrackStep) {
        const std::string fileName = "test_monitor_checkpoint_missing_step.h5";
        h5_file_t file             = openCollectiveH5(fileName, H5_O_WRONLY);
        h5_float64_t time          = 1.0e-12;
        ASSERT_EQ(H5SetStep(file, 0), H5_SUCCESS);
        ASSERT_EQ(H5WriteStepAttribFloat64(file, "TIME", &time, 1), H5_SUCCESS);
        ASSERT_EQ(H5CloseFile(file), H5_SUCCESS);

        EXPECT_THROW(LossDataSink::rewindH5ToGlobalTrackStep(fileName, 1), GeneralOpalException);
    }

    TEST_F(LossDataSinkTest, RewindH5CompactsPhysicalStorageAcrossRestartCycles) {
        const std::string fileName        = "test_monitor_checkpoint_compaction.h5";
        constexpr std::size_t payloadSize = 32768;
        writeGlobalTrackSteps(fileName, {4, 8, 12}, payloadSize);
        const auto originalSize = std::filesystem::file_size(fileName);

        LossDataSink::rewindH5ToGlobalTrackStep(fileName, 8);
        const auto firstCompactedSize = std::filesystem::file_size(fileName);
        EXPECT_LT(firstCompactedSize, originalSize);
        EXPECT_EQ(readRootMetadata(fileName), 17);
        EXPECT_EQ(readGlobalTrackSteps(fileName), std::vector<h5_int64_t>({4}));

        appendGlobalTrackStep(fileName, 8, payloadSize);
        appendGlobalTrackStep(fileName, 12, payloadSize);
        const auto firstRestartSize = std::filesystem::file_size(fileName);

        LossDataSink::rewindH5ToGlobalTrackStep(fileName, 8);
        const auto secondCompactedSize = std::filesystem::file_size(fileName);
        EXPECT_EQ(secondCompactedSize, firstCompactedSize);

        appendGlobalTrackStep(fileName, 8, payloadSize);
        appendGlobalTrackStep(fileName, 12, payloadSize);
        EXPECT_EQ(std::filesystem::file_size(fileName), firstRestartSize);
        EXPECT_EQ(readGlobalTrackSteps(fileName), std::vector<h5_int64_t>({4, 8, 12}));
    }

}  // namespace
