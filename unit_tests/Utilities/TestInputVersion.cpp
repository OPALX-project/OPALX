#include "Utilities/InputVersion.h"
#include "Utilities/OpalException.h"

#include <gtest/gtest.h>

#include <array>
#include <optional>

using InputVersion::Change;
using InputVersion::Status;
using InputVersion::Version;

namespace {
    constexpr Version build{1, 3};
    constexpr std::array<Change, 3> table{
            Change{{1, 1}, "first"}, Change{{1, 2}, "second"}, Change{{1, 3}, "third"}};
}  // namespace

TEST(InputVersionTest, ParsesMajorMinorAndIgnoresPatch) {
    EXPECT_EQ(InputVersion::parse("1.1"), (Version{1, 1}));
    EXPECT_EQ(InputVersion::parse("2.10"), (Version{2, 10}));
    EXPECT_EQ(InputVersion::parse("1.1.7"), (Version{1, 1}));
    EXPECT_EQ(InputVersion::toString(Version{2, 10}), "2.10");
}

TEST(InputVersionTest, RejectsOtherText) {
    for (const char* text :
         {"", "1", "10900", "1.", ".1", "1.1.", "1.1.1.1", "v1.1", "1.1 ", " 1.1", "1,1", "-1.1",
          "1.x", "12345.1"}) {
        EXPECT_THROW(InputVersion::parse(text), OpalException) << "text: '" << text << "'";
    }
}

TEST(InputVersionTest, OrdersMinorNumerically) {
    EXPECT_LT((Version{1, 9}), (Version{1, 10}));
    EXPECT_LT((Version{1, 10}), (Version{2, 0}));
}

TEST(InputVersionTest, MissingVersionIsReported) {
    EXPECT_EQ(InputVersion::evaluate(std::nullopt, build, table).status, Status::NOTGIVEN);
}

TEST(InputVersionTest, CurrentVersionPasses) {
    const auto result = InputVersion::evaluate(build, build, table);
    EXPECT_EQ(result.status, Status::CURRENT);
    EXPECT_TRUE(result.changes.empty());
}

TEST(InputVersionTest, NewerVersionIsReported) {
    EXPECT_EQ(InputVersion::evaluate(Version{1, 4}, build, table).status, Status::NEWER);
    EXPECT_EQ(InputVersion::evaluate(Version{2, 0}, build, table).status, Status::NEWER);
}

TEST(InputVersionTest, OlderVersionListsOnlyNewerChanges) {
    const auto result = InputVersion::evaluate(Version{1, 1}, build, table);
    ASSERT_EQ(result.status, Status::OLDERAFFECTED);
    ASSERT_EQ(result.changes.size(), 2u);
    EXPECT_EQ(result.changes[0].text, "second");
    EXPECT_EQ(result.changes[1].text, "third");
}

TEST(InputVersionTest, OlderVersionWithoutNewerChangesPasses) {
    constexpr std::array<Change, 1> early{Change{{1, 1}, "first"}};
    const auto result = InputVersion::evaluate(Version{1, 2}, build, early);
    EXPECT_EQ(result.status, Status::OLDERUNAFFECTED);
    EXPECT_TRUE(result.changes.empty());
}

// A release bump alone must never stop an input: with an empty table every version runs.
TEST(InputVersionTest, EmptyTableNeverAffectsAnInput) {
    const std::span<const Change> none;
    for (const Version input : {Version{0, 1}, Version{1, 0}, Version{1, 2}}) {
        EXPECT_EQ(InputVersion::evaluate(input, build, none).status, Status::OLDERUNAFFECTED);
    }
}

TEST(InputVersionTest, TableGuardRejectsEntriesAheadOfTheBuildOrOutOfOrder) {
    EXPECT_TRUE(InputVersion::isValidTable(table, build));
    EXPECT_TRUE(InputVersion::isValidTable({}, build));
    EXPECT_FALSE(InputVersion::isValidTable(table, Version{1, 2}));

    constexpr std::array<Change, 2> unordered{Change{{1, 2}, "second"}, Change{{1, 1}, "first"}};
    EXPECT_FALSE(InputVersion::isValidTable(unordered, build));
}

TEST(InputVersionTest, RegisteredTableMatchesTheBuild) {
    EXPECT_TRUE(
            InputVersion::isValidTable(
                    InputVersion::registeredChanges(), InputVersion::currentVersion()));
}
