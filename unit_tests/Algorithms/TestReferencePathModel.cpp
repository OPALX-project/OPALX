#include "gtest/gtest.h"

#include "Algorithms/ReferencePathModel.h"
#include "BeamlineCore/DriftRep.h"

TEST(ReferencePathModelTest, PreservesLegacyElementEdgeMetadata) {
    ReferencePathModel model;
    auto drift = std::make_shared<DriftRep>("D1");
    model.addSegment(ReferencePathSegment(1.0, 2.5, {drift}, 1.25));
    model.addSegment(ReferencePathSegment(2.5, 4.0));

    ASSERT_EQ(model.size(), 2u);
    EXPECT_DOUBLE_EQ(model[0].getBegin(), 1.0);
    EXPECT_DOUBLE_EQ(model[0].getEnd(), 2.5);
    ASSERT_EQ(model[0].getActiveElements().size(), 1u);
    EXPECT_EQ((*model[0].getActiveElements().begin())->getName(), "D1");
    ASSERT_TRUE(model[0].hasLegacyElementEdge());
    EXPECT_DOUBLE_EQ(*model[0].getLegacyElementEdge(), 1.25);
    EXPECT_FALSE(model[1].hasLegacyElementEdge());
}
