#include <filesystem>
#include <fstream>
#include <limits>
#include <vector>
#include "AbsBeamline/EndFieldModel/EndFieldModelManager.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Ippl.h"
#include "OpalConfigure/Configure.h"
#include "Track/CofCmd.h"
#include "Track/MatchCmd.h"
#include "Utilities/ResultJson.h"
#include "gtest/gtest.h"

namespace {
    std::string read(const std::filesystem::path& path) {
        std::ifstream input(path);
        return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
    }
    class MatchCommandTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }
        static void TearDownTestSuite() {
            // Match Main.cpp: Configure registers Kokkos-backed field models.
            endfieldmodel::EndFieldModelManager::clearEFMManager();
            ippl::finalize();
        }
        void SetUp() override { Configure::configure(); }
        void TearDown() override { OpalData::deleteInstance(); }
    };
}  // namespace

TEST_F(MatchCommandTest, RegistrationAttributesAndUnexecutedClones) {
    auto* data  = OpalData::getInstance();
    auto* model = dynamic_cast<MatchCmd*>(data->find("MATCH"));
    ASSERT_NE(model, nullptr);
    for (const char* name : {"ORBIT", "MODEEMITTANCES_N", "OUTPUT"})
        EXPECT_NE(model->findAttribute(name), nullptr) << name;
    EXPECT_TRUE(model->findAttribute("OUTPUT")->defaultUsed());
    EXPECT_TRUE(Attributes::getString(*model->findAttribute("OUTPUT")).empty());
    EXPECT_TRUE(Attributes::getRealArray(*model->findAttribute("MODEEMITTANCES_N")).empty());
    auto* first = model->clone("M1");
    data->define(first);
    Attributes::setString(*first->findAttribute("ORBIT"), "UNSOLVED");
    Attributes::setRealArray(*first->findAttribute("MODEEMITTANCES_N"), {2e-6, 7e-6});
    auto* copy = first->clone("M2");
    data->define(copy);
    Attributes::setRealArray(*copy->findAttribute("MODEEMITTANCES_N"), {3e-6, 8e-6});
    EXPECT_EQ(
            Attributes::getRealArray(*first->findAttribute("MODEEMITTANCES_N")),
            (std::vector<double>{2e-6, 7e-6}));
    EXPECT_THROW(MatchCmd::findResult("MATCH"), OpalException);
    EXPECT_THROW(MatchCmd::findResult("M1"), OpalException);
    EXPECT_THROW(MatchCmd::findResult("M2"), OpalException);
    auto* cof = dynamic_cast<CofCmd*>(data->find("COF"));
    ASSERT_NE(cof, nullptr);
    data->define(cof->clone("UNSOLVED"));
    EXPECT_THROW(CofCmd::findOptics("UNSOLVED"), OpalException);
    EXPECT_THROW(first->execute(), OpalException);
    EXPECT_THROW(MatchCmd::findResult("M1"), OpalException);
    EXPECT_THROW(MatchCmd::findResult("M2"), OpalException);
}

TEST_F(MatchCommandTest, RejectsMissingWrongTypeAndUnexecutedNamedResults) {
    auto* data  = OpalData::getInstance();
    auto* model = dynamic_cast<MatchCmd*>(data->find("MATCH"));
    ASSERT_NE(model, nullptr);
    auto* command = model->clone("CANDIDATE");
    data->define(command);
    Attributes::setRealArray(*command->findAttribute("MODEEMITTANCES_N"), {2e-6, 7e-6});
    for (const char* name : {"MISSING", "COF", "MATCH", "CANDIDATE"}) {
        Attributes::setString(*command->findAttribute("ORBIT"), name);
        EXPECT_THROW(command->execute(), OpalException) << name;
        EXPECT_THROW(CofCmd::findOptics(name), OpalException) << name;
        EXPECT_THROW(MatchCmd::findResult(name), OpalException) << name;
        EXPECT_THROW(MatchCmd::findResult("CANDIDATE"), OpalException);
    }
}

TEST_F(MatchCommandTest, JsonEscapingAndFiniteNumbers) {
    opalx::report::ResultJsonWriter json;
    json.object(nullptr, [&] {
        json.text("quote", "a\"b\\c\n");
        json.number("value", .12345678901234567);
    });
    EXPECT_NE(json.str().find("a\\\"b\\\\c\\u000a"), std::string::npos);
    EXPECT_NE(json.str().find("0.12345678901234566"), std::string::npos);
    EXPECT_THROW(json.number("bad", std::numeric_limits<double>::infinity()), OpalException);
}

TEST_F(MatchCommandTest, AtomicReplacementAndFailureKeepPreviousFile) {
    // Independent C++ filesystem fixture on each rank; no OPAL input or ray solve.
    std::string pattern = (std::filesystem::temp_directory_path() / "opalx-match-XXXXXX").string();
    const auto created  = ::mkdtemp(pattern.data());
    ASSERT_NE(created, nullptr);
    const std::filesystem::path directory(created), target = directory / "old.json";
    struct Cleanup {
        std::filesystem::path path;
        ~Cleanup() {
            std::error_code ignored;
            std::filesystem::remove_all(path, ignored);
        }
    } cleanup{directory};
    opalx::report::writeResult("TEST", target.string(), "old");
    opalx::report::writeResult("TEST", target.string(), "new");
    EXPECT_EQ(read(target), "new");
    EXPECT_THROW(
            opalx::report::writeResult("TEST", (target / "impossible").string(), "lost"),
            OpalException);
    EXPECT_EQ(read(target), "new");
    EXPECT_THROW(opalx::report::writeResult("TEST", directory.string(), "lost"), OpalException);
    for (const auto& entry : std::filesystem::directory_iterator(directory))
        EXPECT_EQ(entry.path().filename(), "old.json");
}
