#include "gtest/gtest.h"

#include "Algorithms/IndexMap.h"
#include "BeamlineCore/DriftRep.h"

#include <limits>
#include <set>
#include <string>

namespace {
    std::set<std::string> names(const IndexMap::value_t& elements) {
        std::set<std::string> result;
        for (const auto& element : elements) {
            result.insert(element->getName());
        }
        return result;
    }
}  // namespace

// query(s, ds) returns the union of all elements whose range overlaps [s - ds, s + ds].
TEST(IndexMapTest, QueryReturnsElementsOverlappingTheInterval) {
    auto drift1 = std::make_shared<DriftRep>("D1");
    auto drift2 = std::make_shared<DriftRep>("D2");

    IndexMap map;
    map.add(0.0, 1.0, IndexMap::value_t{drift1});
    map.add(1.0, 2.0, IndexMap::value_t{drift1, drift2});

    EXPECT_EQ(names(map.query(0.5, 0.1)), (std::set<std::string>{"D1"}));
    EXPECT_EQ(names(map.query(1.5, 0.1)), (std::set<std::string>{"D1", "D2"}));
}

// getRange(element) returns the element's single (merged) range. Consecutive contiguous adds
// for the same element are merged into one interval.
TEST(IndexMapTest, GetRangeReturnsElementSingleRange) {
    auto drift1 = std::make_shared<DriftRep>("D1");
    auto drift2 = std::make_shared<DriftRep>("D2");

    IndexMap map;
    map.add(0.0, 1.0, IndexMap::value_t{drift1});
    map.add(1.0, 2.0, IndexMap::value_t{drift1, drift2});

    const double end = 2.0 * (1.0 - std::numeric_limits<double>::epsilon());

    IndexMap::key_t r1 = map.getRange(drift1);
    EXPECT_DOUBLE_EQ(r1.begin, 0.0);
    EXPECT_DOUBLE_EQ(r1.end, end);

    IndexMap::key_t r2 = map.getRange(drift2);
    EXPECT_DOUBLE_EQ(r2.begin, 1.0);
    EXPECT_DOUBLE_EQ(r2.end, end);
}
