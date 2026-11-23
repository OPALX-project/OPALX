#include <gtest/gtest.h>

#include "Ippl.h"

#include "Algorithms/StepSizeConfig.h"
#include "Utilities/OpalException.h"

#include <string>

TEST(StepSizeConfigTest, AdvanceToIndexSelectsSegmentAndEnd) {
    StepSizeConfig config;
    config.push_back(1.0e-12, 1.0, 10);
    config.push_back(2.0e-12, 2.0, 20);
    config.push_back(3.0e-12, 3.0, 30);

    config.advanceToIndex(1);
    EXPECT_EQ(config.getCurrentIndex(), 1u);
    EXPECT_DOUBLE_EQ(config.getdT(), 2.0e-12);
    EXPECT_DOUBLE_EQ(config.getSStop(), 2.0);
    EXPECT_EQ(config.getNumSteps(), 20u);

    config.advanceToIndex(3);
    EXPECT_TRUE(config.reachedEnd());
    EXPECT_EQ(config.getCurrentIndex(), 3u);
}

TEST(StepSizeConfigTest, AdvanceToIndexRejectsOutOfRangeIndex) {
    StepSizeConfig config;
    config.push_back(1.0e-12, 1.0, 10);

    EXPECT_THROW(config.advanceToIndex(2), OpalException);
}

TEST(StepSizeConfigTest, AdvanceToResumePositionUsesStoredSegmentAfterEarlyStop) {
    StepSizeConfig config;
    config.push_back(1.0e-12, 1.0, 10);
    config.push_back(2.0e-12, 2.0, 20);

    // Segment 0 may have reached ZSTOP after only five of its ten allowed steps. A checkpoint
    // after two more global steps must resume segment 1 at offset two, even though global step
    // seven would still fall inside segment 0 if inferred only from MAXSTEPS.
    config.advanceToResumePosition({1, 2});

    EXPECT_EQ(config.getCurrentIndex(), 1u);
    EXPECT_DOUBLE_EQ(config.getdT(), 2.0e-12);
    EXPECT_EQ(config.getNumSteps(), 20u);
}

TEST(StepSizeConfigTest, AdvanceToResumePositionAcceptsEndOfSchedule) {
    StepSizeConfig config;
    config.push_back(1.0e-12, 1.0, 10);
    config.push_back(2.0e-12, 2.0, 20);

    config.advanceToResumePosition({2, 0});

    EXPECT_TRUE(config.reachedEnd());
}

TEST(StepSizeConfigTest, AdvanceToResumePositionRejectsInvalidOffsets) {
    StepSizeConfig config;
    config.push_back(1.0e-12, 1.0, 10);
    config.push_back(2.0e-12, 2.0, 20);

    EXPECT_THROW(config.advanceToResumePosition({0, 10}), OpalException);
    EXPECT_THROW(config.advanceToResumePosition({2, 1}), OpalException);
    EXPECT_THROW(config.advanceToResumePosition({3, 0}), OpalException);
}

// resumeFromCheckpoint() logs through Inform, which needs IPPL.
class StepSizeConfigResumeTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }
    static void TearDownTestSuite() { ippl::finalize(); }

    static void fill(StepSizeConfig& config) {
        config.push_back(1.0e-12, 1.0, 10);
        config.push_back(2.0e-12, 2.0, 20);
        config.resetIterator();
    }
};

TEST_F(StepSizeConfigResumeTest, ResumesAtTheCheckpointSegment) {
    StepSizeConfig config;
    fill(config);

    config.resumeFromCheckpoint({1, 4}, 2.0e-12, 14);
    EXPECT_EQ(config.getCurrentIndex(), 1u);
    EXPECT_EQ(config.getdT(), 2.0e-12);
    EXPECT_EQ(config.getSStop(), 2.0);

    // A checkpoint time step that differs from the schedule's only by rounding is the same.
    StepSizeConfig rounded;
    fill(rounded);
    EXPECT_NO_THROW(rounded.resumeFromCheckpoint({1, 4}, 2.0e-12 * (1.0 + 1.0e-15), 14));
    EXPECT_EQ(rounded.getCurrentIndex(), 1u);
}

TEST_F(StepSizeConfigResumeTest, RejectsACheckpointWhoseTimeStepTheInputChanged) {
    StepSizeConfig config;
    fill(config);

    try {
        config.resumeFromCheckpoint({1, 4}, 1.0e-12, 14);
        FAIL() << "a changed time step at the checkpoint position was accepted";
    } catch (const OpalException& error) {
        EXPECT_NE(error.what().find("changes the time-step schedule"), std::string::npos)
                << error.what();
    }
}

TEST_F(StepSizeConfigResumeTest, PositionAtTheEndOfTheScheduleIsFinished) {
    StepSizeConfig config;
    fill(config);

    // No segment is left, so there is no time step to compare with.
    EXPECT_NO_THROW(config.resumeFromCheckpoint({2, 0}, 5.0e-12, 30));
    EXPECT_TRUE(config.reachedEnd());
}
