#include "gtest/gtest.h"

#include "AbstractObjects/OpalData.h"
#include "Algorithms/IndexMap.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/RFCavityRep.h"
#include "Ippl.h"

#include <filesystem>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

extern Inform* gmsg;

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

// Crossing an element a second time (ring) keeps the first crossing in getRange, while query
// still finds the element at every crossing.
TEST(IndexMapTest, ReentryKeepsFirstRangeAndQueryFindsBothCrossings) {
    auto drift1 = std::make_shared<DriftRep>("D1");

    IndexMap map;
    map.add(0.0, 1.0, IndexMap::value_t{drift1});
    map.add(1.0, 2.0, IndexMap::value_t{drift1});  // continuation, merges
    map.add(5.0, 6.0, IndexMap::value_t{drift1});  // re-entry
    map.add(6.0, 7.0, IndexMap::value_t{drift1});

    IndexMap::key_t r1 = map.getRange(drift1);
    EXPECT_DOUBLE_EQ(r1.begin, 0.0);
    EXPECT_DOUBLE_EQ(r1.end, 2.0 * (1.0 - std::numeric_limits<double>::epsilon()));

    EXPECT_EQ(names(map.query(0.5, 0.1)), (std::set<std::string>{"D1"}));
    EXPECT_EQ(names(map.query(6.5, 0.1)), (std::set<std::string>{"D1"}));
}

TEST(IndexMapTest, PeriodicQueryReusesOneTurnMap) {
    auto start = std::make_shared<DriftRep>("START");
    auto end   = std::make_shared<DriftRep>("END");

    IndexMap map;
    map.setPeriod(10.0, 4.0);
    map.add(10.0, 11.0, IndexMap::value_t{start});
    map.add(11.0, 13.0, IndexMap::value_t{});
    map.add(13.0, 14.0, IndexMap::value_t{end});

    EXPECT_EQ(names(map.query(10.5, 0.1)), (std::set<std::string>{"START"}));
    EXPECT_EQ(names(map.query(22.5, 0.1)), (std::set<std::string>{"START"}));
    EXPECT_EQ(names(map.query(9.5, 0.1)), (std::set<std::string>{"END"}));
}

TEST(IndexMapTest, PeriodicQueryUnionsBothSidesOfTurnSeam) {
    auto start = std::make_shared<DriftRep>("START");
    auto end   = std::make_shared<DriftRep>("END");

    IndexMap map;
    map.setPeriod(0.0, 4.0);
    map.add(0.0, 1.0, IndexMap::value_t{start});
    map.add(1.0, 3.0, IndexMap::value_t{});
    map.add(3.0, 4.0, IndexMap::value_t{end});

    EXPECT_EQ(names(map.query(4.0, 0.2)), (std::set<std::string>{"END", "START"}));
    EXPECT_EQ(names(map.query(12.0, 0.2)), (std::set<std::string>{"END", "START"}));
    EXPECT_EQ(names(map.query(2.0, 2.0)), (std::set<std::string>{"END", "START"}));
}

TEST(IndexMapTest, PeriodicCoordinateReportsTurnAndTopologicalPhase) {
    IndexMap map;
    map.setPeriod(10.0, 4.0);

    EXPECT_EQ(map.getTurnNumber(9.9), -1);
    EXPECT_EQ(map.getTurnNumber(10.0), 0);
    EXPECT_EQ(map.getTurnNumber(13.99), 0);
    EXPECT_EQ(map.getTurnNumber(14.0), 1);
    EXPECT_NEAR(map.getPhase(11.0), 0.5 * std::acos(-1.0), 1.0e-14);
    EXPECT_NEAR(map.getPhase(15.0), 0.5 * std::acos(-1.0), 1.0e-14);
}

// saveSDDS needs OpalData (input file name, output directory) and gmsg.
class IndexMapSDDSTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        gmsg = new Inform(nullptr, -1);
        std::filesystem::create_directories(OpalData::getInstance()->getAuxiliaryOutputDirectory());
    }

    static void TearDownTestSuite() {
        delete gmsg;
        gmsg = nullptr;
        ippl::finalize();
    }

    void SetUp() override {
        OpalData::getInstance()->storeInputFn("indexmap_sdds.opal");
        OpalData::getInstance()->setOpenMode(OpalData::OpenMode::WRITE);
    }

    void TearDown() override { std::filesystem::remove(sddsPath()); }

    static std::filesystem::path sddsPath() {
        return std::filesystem::path(OpalData::getInstance()->getAuxiliaryOutputDirectory())
               / "indexmap_sdds_ElementPositions.sdds";
    }
};

// With two crossings per element the file must contain the first crossing only: the map
// stores one interval per element, so rows past the first re-entry would be drawn with
// wrong intervals.
TEST_F(IndexMapSDDSTest, SaveSDDSWritesOnlyFirstCrossingOfEachElement) {
    auto drift  = std::make_shared<DriftRep>("D1");
    auto cavity = std::make_shared<RFCavityRep>("RF1");

    IndexMap map;
    // Turn 1: drift over [0,1], cavity over [1,2].
    map.add(0.0, 1.0, IndexMap::value_t{drift});
    map.add(1.0, 2.0, IndexMap::value_t{cavity});
    // Turn 2: same elements again.
    map.add(2.0, 3.0, IndexMap::value_t{drift});
    map.add(3.0, 4.0, IndexMap::value_t{cavity});
    map.saveSDDS(0.0);

    std::ifstream sdds(sddsPath());
    ASSERT_TRUE(sdds.is_open());

    bool foundCavityRow = false;
    std::vector<double> sValues;
    std::string line;
    while (std::getline(sdds, line)) {
        // Data rows start with the s value and end with the quoted element_names column;
        // header lines start with '&' and parameter lines carry no quotes.
        if (line.empty() || line[0] == '&' || line.find('"') == std::string::npos) {
            continue;
        }
        std::istringstream row(line);
        double s = 0.0;
        if (!(row >> s)) {
            continue;
        }
        sValues.push_back(s);
        if (line.find("\"RF1\"") != std::string::npos) {
            foundCavityRow = true;
        }
    }
    EXPECT_TRUE(foundCavityRow);
    // Exactly the two first-turn sectors, four rows each: (s_i, s_i, s_f, s_f). The re-entry
    // at s = 2 and everything after it must not be written.
    const std::vector<double> expected{0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0};
    ASSERT_EQ(sValues.size(), expected.size());
    for (size_t i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(sValues[i], expected[i], 1e-9) << "row " << i;
    }
}
