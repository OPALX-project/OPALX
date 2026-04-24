#include "gtest/gtest.h"

#include "Algorithms/ReferencePathModel.h"

TEST(ReferencePathModelTest, PreservesLegacyElementEdgeMetadata) {
    ReferencePathModel model;
    model.addSegment(ReferencePathSegment(1.0, 2.5, 1.25));
    model.addSegment(ReferencePathSegment(2.5, 4.0));

    ASSERT_EQ(model.size(), 2u);
    EXPECT_DOUBLE_EQ(model[0].getBegin(), 1.0);
    EXPECT_DOUBLE_EQ(model[0].getEnd(), 2.5);
    ASSERT_TRUE(model[0].hasLegacyElementEdge());
    EXPECT_DOUBLE_EQ(*model[0].getLegacyElementEdge(), 1.25);
    EXPECT_FALSE(model[1].hasLegacyElementEdge());
}
