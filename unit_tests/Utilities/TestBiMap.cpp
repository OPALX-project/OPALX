/**
 * \file TestBiMap.cpp
 * \brief Unit tests for BiMap class (replacement for boost::bimap)
 *
 * This file contains comprehensive unit tests for the BiMap bidirectional map
 * implementation. The BiMap class provides a two-way mapping between left and
 * right value types, allowing efficient lookup in both directions.
 *
 * \test BiMapTest::BasicInsertion
 * Tests basic insertion and bidirectional lookup functionality. Verifies that
 * values can be inserted and retrieved from both left-to-right and right-to-left
 * directions.
 *
 * \test BiMapTest::LeftViewAccess
 * Tests the left view accessor, which provides STL-compatible iterator interface
 * for accessing the left-to-right mapping. Verifies find() and end() operations.
 *
 * \test BiMapTest::RightViewAccess
 * Tests the right view accessor, which provides STL-compatible iterator interface
 * for accessing the right-to-left mapping. Verifies find() and end() operations.
 *
 * \test BiMapTest::OutOfRangeException
 * Tests exception handling when accessing non-existent keys. Verifies that
 * \c std::out_of_range exceptions are thrown for both left_at() and right_at()
 * when keys are not found.
 *
 * \test BiMapTest::MakeBiMapHelper
 * Tests the make_bimap() helper function for creating a BiMap from an
 * initializer list. Verifies that the helper function correctly constructs
 * the bidirectional mapping.
 *
 * \test BiMapTest::FindMethods
 * Tests the find methods (left_find() and right_find()) and their corresponding
 * end iterators. Verifies that find operations return correct iterators for
 * existing keys and end iterators for non-existent keys.
 *
 * \test BiMapTest::StringToString
 * Tests BiMap with string types for both left and right values. Verifies that
 * the bidirectional mapping works correctly with string keys and values.
 */

#include <gtest/gtest.h>
#include <string>
#include "Utilities/BiMap.h"

class BiMapTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Test setup if needed
    }
};

TEST_F(BiMapTest, BasicInsertion) {
    BiMap<int, std::string> bimap;
    bimap.insert(1, "one");
    bimap.insert(2, "two");
    bimap.insert(3, "three");

    EXPECT_EQ(bimap.left_at(1), "one");
    EXPECT_EQ(bimap.left_at(2), "two");
    EXPECT_EQ(bimap.left_at(3), "three");

    EXPECT_EQ(bimap.right_at("one"), 1);
    EXPECT_EQ(bimap.right_at("two"), 2);
    EXPECT_EQ(bimap.right_at("three"), 3);
}

TEST_F(BiMapTest, LeftViewAccess) {
    BiMap<int, std::string> bimap;
    bimap.insert(1, "one");
    bimap.insert(2, "two");

    auto it = bimap.left.find(1);
    EXPECT_NE(it, bimap.left.end());
    EXPECT_EQ(it->second, "one");

    it = bimap.left.find(3);
    EXPECT_EQ(it, bimap.left.end());
}

TEST_F(BiMapTest, RightViewAccess) {
    BiMap<int, std::string> bimap;
    bimap.insert(1, "one");
    bimap.insert(2, "two");

    auto it = bimap.right.find("one");
    EXPECT_NE(it, bimap.right.end());
    EXPECT_EQ(it->second, 1);

    it = bimap.right.find("three");
    EXPECT_EQ(it, bimap.right.end());
}

TEST_F(BiMapTest, OutOfRangeException) {
    BiMap<int, std::string> bimap;
    bimap.insert(1, "one");

    EXPECT_THROW(bimap.left_at(2), std::out_of_range);
    EXPECT_THROW(bimap.right_at("two"), std::out_of_range);
}

TEST_F(BiMapTest, MakeBiMapHelper) {
    auto bimap = make_bimap<int, std::string>({{1, "one"}, {2, "two"}, {3, "three"}});

    EXPECT_EQ(bimap.left_at(1), "one");
    EXPECT_EQ(bimap.left_at(2), "two");
    EXPECT_EQ(bimap.left_at(3), "three");

    EXPECT_EQ(bimap.right_at("one"), 1);
    EXPECT_EQ(bimap.right_at("two"), 2);
    EXPECT_EQ(bimap.right_at("three"), 3);
}

TEST_F(BiMapTest, FindMethods) {
    BiMap<int, std::string> bimap;
    bimap.insert(1, "one");
    bimap.insert(2, "two");

    auto left_it = bimap.left_find(1);
    EXPECT_NE(left_it, bimap.left_end());
    EXPECT_EQ(left_it->second, "one");

    auto right_it = bimap.right_find("two");
    EXPECT_NE(right_it, bimap.right_end());
    EXPECT_EQ(right_it->second, 2);

    left_it = bimap.left_find(99);
    EXPECT_EQ(left_it, bimap.left_end());

    right_it = bimap.right_find("ninety-nine");
    EXPECT_EQ(right_it, bimap.right_end());
}

TEST_F(BiMapTest, StringToString) {
    BiMap<std::string, std::string> bimap;
    bimap.insert("key1", "value1");
    bimap.insert("key2", "value2");

    EXPECT_EQ(bimap.left_at("key1"), "value1");
    EXPECT_EQ(bimap.right_at("value2"), "key2");
}
