#include "gtest/gtest.h"

#include "Algorithms/IndexMap.h"
#include "BeamlineCore/DriftRep.h"

#include <limits>

TEST(IndexMapTest, RebuildsReferencePathModelFromOrderedRanges) {
    auto drift1 = std::make_shared<DriftRep>("D1");
    auto drift2 = std::make_shared<DriftRep>("D2");

    IndexMap map;
    map.add(0.0, 1.0, IndexMap::value_t{drift1});
    map.add(1.0, 2.0, IndexMap::value_t{drift1, drift2});

    const ReferencePathModel& model = map.getReferencePathModel();

    ASSERT_EQ(model.size(), 2u);
    EXPECT_DOUBLE_EQ(model[0].getBegin(), 0.0);
    EXPECT_DOUBLE_EQ(model[0].getEnd(), 1.0 * (1.0 - std::numeric_limits<double>::epsilon()));
    ASSERT_EQ(model[0].getActiveElements().size(), 1u);
    EXPECT_EQ((*model[0].getActiveElements().begin())->getName(), "D1");

    EXPECT_DOUBLE_EQ(model[1].getBegin(), 1.0);
    EXPECT_DOUBLE_EQ(model[1].getEnd(), 2.0 * (1.0 - std::numeric_limits<double>::epsilon()));
    ASSERT_EQ(model[1].getActiveElements().size(), 2u);
}
