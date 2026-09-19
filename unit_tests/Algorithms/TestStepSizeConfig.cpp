#include <gtest/gtest.h>

#include "Algorithms/StepSizeConfig.h"
#include "Utilities/OpalException.h"

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
