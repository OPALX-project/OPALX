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
