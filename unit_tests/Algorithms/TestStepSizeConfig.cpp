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

TEST(StepSizeConfigTest, AdvanceToGlobalStepMapsIntoExtendedSingleSegment) {
    StepSizeConfig config;
    config.push_back(1.0e-12, 1.0, 20);

    const auto position = config.advanceToGlobalStep(15);

    EXPECT_EQ(position.segment, 0u);
    EXPECT_EQ(position.stepsCompletedInSegment, 15u);
    EXPECT_EQ(config.getCurrentIndex(), 0u);
    EXPECT_DOUBLE_EQ(config.getdT(), 1.0e-12);
    EXPECT_EQ(config.getNumSteps(), 20u);
}

TEST(StepSizeConfigTest, AdvanceToGlobalStepMapsAcrossSegmentBoundaries) {
    StepSizeConfig config;
    config.push_back(1.0e-12, 1.0, 10);
    config.push_back(2.0e-12, 2.0, 20);
    config.push_back(3.0e-12, 3.0, 30);

    auto position = config.advanceToGlobalStep(10);
    EXPECT_EQ(position.segment, 1u);
    EXPECT_EQ(position.stepsCompletedInSegment, 0u);
    EXPECT_EQ(config.getCurrentIndex(), 1u);
    EXPECT_DOUBLE_EQ(config.getdT(), 2.0e-12);

    position = config.advanceToGlobalStep(12);
    EXPECT_EQ(position.segment, 1u);
    EXPECT_EQ(position.stepsCompletedInSegment, 2u);
    EXPECT_EQ(config.getCurrentIndex(), 1u);
    EXPECT_DOUBLE_EQ(config.getdT(), 2.0e-12);

    position = config.advanceToGlobalStep(60);
    EXPECT_EQ(position.segment, 3u);
    EXPECT_EQ(position.stepsCompletedInSegment, 0u);
    EXPECT_TRUE(config.reachedEnd());
}

TEST(StepSizeConfigTest, AdvanceToGlobalStepRejectsStepsPastSchedule) {
    StepSizeConfig config;
    config.push_back(1.0e-12, 1.0, 10);

    EXPECT_THROW(config.advanceToGlobalStep(11), OpalException);
}
