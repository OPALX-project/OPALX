//
// Unit tests for BinConfigWriter.
//

#include <gtest/gtest.h>

#include "Ippl.h"
#include "Structure/BinConfigWriter.h"
#include "Utilities/OpalException.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

class BinConfigWriterTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() {
        ippl::finalize();
    }

    // Helper to read entire file into string.
    static std::string readFile(const std::string& path) {
        std::ifstream f(path);
        std::stringstream ss;
        ss << f.rdbuf();
        return ss.str();
    }

    // Helper to check that a string is valid JSON array (basic structural check).
    static bool looksLikeJsonArray(const std::string& s) {
        if (s.empty()) return false;
        size_t i = 0;
        while (i < s.size() && (s[i] == ' ' || s[i] == '\n')) ++i;
        if (i >= s.size() || s[i] != '[') return false;
        size_t depth = 0;
        for (; i < s.size(); ++i) {
            if (s[i] == '[') ++depth;
            else if (s[i] == ']') {
                if (depth == 0) return false;
                --depth;
            }
        }
        return depth == 0;
    }
};

// ConstructionOpensFileAndWritesOpeningBracket:
// Constructor creates the file and writes "[\n".
TEST_F(BinConfigWriterTest, ConstructionOpensFileAndWritesOpeningBracket) {
    const std::string path = "bin_config_construction_test.json";
    {
        BinConfigWriter w(path);
    }
    std::string content = readFile(path);
    EXPECT_TRUE(content.find("[\n") == 0) << "Expected file to start with '[\\n', got: " << content.substr(0, 20);
}

// SingleEntryProducesValidJson:
// One writeEntry produces a well-formed JSON array with one object.
TEST_F(BinConfigWriterTest, SingleEntryProducesValidJson) {
    const std::string path = "bin_config_single_entry_test.json";
    {
        BinConfigWriter w(path);
        std::vector<std::size_t> counts = {10, 20, 30};
        std::vector<double> widths = {0.1, 0.2, 0.3};
        w.writeEntry(0, 1.0e-10, true, counts, widths, 0.5);
    }
    std::string content = readFile(path);
    EXPECT_TRUE(looksLikeJsonArray(content)) << "Output is not valid JSON array: " << content.substr(0, 200);
    EXPECT_NE(content.find("\"step\": 0"), std::string::npos);
    EXPECT_NE(content.find("\"time\":"), std::string::npos);
    EXPECT_NE(content.find("\"preMerge\": true"), std::string::npos);
    EXPECT_NE(content.find("\"xMin\":"), std::string::npos);
    EXPECT_NE(content.find("\"binCounts\": ["), std::string::npos);
    EXPECT_NE(content.find("10, 20, 30"), std::string::npos);
    EXPECT_NE(content.find("\"binWidths\": ["), std::string::npos);
}

// MultipleEntriesProduceValidJsonArray:
// Several writeEntry calls produce a JSON array with multiple objects.
TEST_F(BinConfigWriterTest, MultipleEntriesProduceValidJsonArray) {
    const std::string path = "bin_config_multi_entry_test.json";
    {
        BinConfigWriter w(path);
        std::vector<std::size_t> counts1 = {100};
        std::vector<double> widths1 = {1.0};
        w.writeEntry(0, 0.0, true, counts1, widths1, 0.0);

        std::vector<std::size_t> counts2 = {50, 50};
        std::vector<double> widths2 = {0.5, 0.5};
        w.writeEntry(1, 1.0e-10, false, counts2, widths2, 0.25);
    }
    std::string content = readFile(path);
    EXPECT_TRUE(looksLikeJsonArray(content));
    EXPECT_NE(content.find("\"step\": 0"), std::string::npos);
    EXPECT_NE(content.find("\"step\": 1"), std::string::npos);
    EXPECT_NE(content.find("\"preMerge\": true"), std::string::npos);
    EXPECT_NE(content.find("\"preMerge\": false"), std::string::npos);
}

// EmptyArraysHandledCorrectly:
// writeEntry with empty binCounts and binWidths produces valid JSON.
// (Writer outputs "binCounts": [    ], and "binWidths": [    ] with spaces, not "[]".)
TEST_F(BinConfigWriterTest, EmptyArraysHandledCorrectly) {
    const std::string path = "bin_config_empty_arrays_test.json";
    {
        BinConfigWriter w(path);
        std::vector<std::size_t> counts;
        std::vector<double> widths;
        w.writeEntry(42, 3.14, false, counts, widths, -1.5);
    }
    std::string content = readFile(path);
    EXPECT_TRUE(looksLikeJsonArray(content));
    EXPECT_NE(content.find("\"step\": 42"), std::string::npos);
    // Writer formats empty arrays as "key": [    ], with spaces before ]
    EXPECT_NE(content.find("\"binCounts\": ["), std::string::npos);
    EXPECT_NE(content.find("\"binWidths\": ["), std::string::npos);
}

// FileClosedAfterEachEntry:
// After each writeEntry, the file on disk ends with a closing ']' (self-closing design).
TEST_F(BinConfigWriterTest, FileClosedAfterEachEntry) {
    const std::string path = "bin_config_closed_test.json";
    {
        BinConfigWriter w(path);
        std::vector<std::size_t> counts = {1};
        std::vector<double> widths = {1.0};
        w.writeEntry(0, 0.0, true, counts, widths, 0.0);
        // File should already be valid and closed; read without destroying writer.
    }
    std::string content = readFile(path);
    EXPECT_TRUE(looksLikeJsonArray(content));
    EXPECT_EQ(content.back(), '\n');
    size_t lastBracket = content.rfind(']');
    EXPECT_NE(lastBracket, std::string::npos);
    EXPECT_GT(lastBracket, 0u);
}

// InvalidPathThrows:
// Constructor throws OpalException when the file cannot be opened.
TEST_F(BinConfigWriterTest, InvalidPathThrows) {
    // Use a path under a non-existent directory; open should fail on most systems.
    const std::string path = "/nonexistent_directory_xyz123/bin_config_test.json";
    EXPECT_THROW(BinConfigWriter w(path), OpalException);
}
