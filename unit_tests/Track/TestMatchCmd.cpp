#include <filesystem>
#include <fstream>
#include <limits>
#include <vector>
#include "AbsBeamline/EndFieldModel/EndFieldModelManager.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Distribution/Distribution.h"
#include "Ippl.h"
#include "OpalConfigure/Configure.h"
#include "OpalParser/SimpleStatement.h"
#include "OpalParser/StringStream.h"
#include "Track/CofCmd.h"
#include "Track/MatchCmd.h"
#include "Track/MatchedSource.h"
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

TEST_F(MatchCommandTest, MomentsAllowsForwardNamesAndOverwritesAllCorrelationInputs) {
    auto* model = dynamic_cast<Distribution*>(OpalData::getInstance()->find("DISTRIBUTION"));
    ASSERT_NE(model, nullptr);
    std::unique_ptr<Distribution> distribution(model->clone("SHAPE"));
    auto* attribute = distribution->findAttribute("MOMENTS");
    ASSERT_NE(attribute, nullptr);
    StringStream stream("future_match");
    SimpleStatement statement("TestMatchCmd", 1);
    for (Token token = stream.readToken(); !token.isEOF(); token = stream.readToken())
        statement.append(token);
    statement.start();
    EXPECT_NO_THROW(attribute->parse(statement, true));
    EXPECT_EQ(distribution->getMomentsName(), "FUTURE_MATCH");
    EXPECT_EQ(OpalData::getInstance()->find("FUTURE_MATCH"), nullptr);
    Attributes::setPredefinedString(*distribution->findAttribute("TYPE"), "MULTIVARIATEGAUSS");
    Attributes::setReal(*distribution->findAttribute("NPARTDIST"), 256);
    Attributes::setReal(*distribution->findAttribute("SIGMAZ"), .001);
    Attributes::setReal(*distribution->findAttribute("SIGMAPZ"), .002);
    Attributes::setRealArray(
            *distribution->findAttribute("CORR"),
            {99.});  // Ignored, not a malformed active covariance.
    Attributes::setReal(*distribution->findAttribute("CORRX"), 99.);
    EXPECT_NO_THROW(distribution->execute());
    EXPECT_EQ(distribution->getNumParticles(), 256u);
    EXPECT_DOUBLE_EQ(distribution->getSigmaR()[2], .001);
    EXPECT_DOUBLE_EQ(distribution->getSigmaP()[2], .002);
    Attributes::setReal(*distribution->findAttribute("NPARTDIST"), 256.5);
    EXPECT_THROW(distribution->execute(), OpalException);
    Attributes::setReal(*distribution->findAttribute("NPARTDIST"), 256);
    Attributes::setReal(*distribution->findAttribute("SIGMAZ"), 0);
    EXPECT_THROW(distribution->execute(), OpalException);
    Attributes::setReal(*distribution->findAttribute("SIGMAZ"), .001);
    Attributes::setReal(*distribution->findAttribute("SIGMAPZ"), 0);
    EXPECT_THROW(distribution->execute(), OpalException);
    Attributes::setReal(*distribution->findAttribute("SIGMAPZ"), .002);
    Attributes::setBool(*distribution->findAttribute("EMITTED"), true);
    EXPECT_THROW(distribution->execute(), OpalException);
    Attributes::setBool(*distribution->findAttribute("EMITTED"), false);
    Attributes::setReal(*distribution->findAttribute("NPARTDIST"), 0);
    EXPECT_THROW(distribution->execute(), OpalException);
    Attributes::setReal(*distribution->findAttribute("NPARTDIST"), 256);
    Attributes::setPredefinedString(*distribution->findAttribute("TYPE"), "GAUSS");
    EXPECT_THROW(distribution->execute(), OpalException);
}

TEST_F(MatchCommandTest, MomentsSourceRejectsDelayedEmissionOffsetsAndImagePlanes) {
    auto* model = dynamic_cast<EmissionSource*>(OpalData::getInstance()->find("EMISSIONSOURCE"));
    ASSERT_NE(model, nullptr);
    std::unique_ptr<EmissionSource> source(model->clone("S"));
    EXPECT_NO_THROW(validateMatchedSource(*source));
    for (const char* name : {"T0", "EKIN", "R0X", "R0Y", "R0Z", "P0X", "P0Y", "P0Z"}) {
        Attributes::setReal(*source->findAttribute(name), 1.);
        EXPECT_THROW(validateMatchedSource(*source), OpalException) << name;
        Attributes::setReal(*source->findAttribute(name), 0.);
    }
    Attributes::setPredefinedString(*source->findAttribute("EMISSIONMODEL"), "ASTRA");
    EXPECT_THROW(validateMatchedSource(*source), OpalException);
    Attributes::setPredefinedString(*source->findAttribute("EMISSIONMODEL"), "NONE");
    for (const char* name : {"ZEROFACE_R0Z", "SHIFTED_GREENS_FUNCTION"}) {
        Attributes::setBool(*source->findAttribute(name), true);
        EXPECT_THROW(validateMatchedSource(*source), OpalException);
        Attributes::setBool(*source->findAttribute(name), false);
    }
}
