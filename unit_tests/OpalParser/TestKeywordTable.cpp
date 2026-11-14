//
// Tests for the keyword table: keywords are kept apart from the input objects, and the input
// forms that refer to them keep working.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#include <gtest/gtest.h>

#include <memory>
#include <sstream>
#include <string>

#include "AbsBeamline/EndFieldModel/EndFieldModelManager.h"
#include "AbstractObjects/Object.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Ippl.h"
#include "OpalConfigure/Configure.h"
#include "OpalParser/OpalParser.h"
#include "OpalParser/Statement.h"
#include "OpalParser/StringStream.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/ParseError.h"
#include "Utility/Inform.h"

extern Inform* gmsg;

namespace {

    class KeywordTableTest : public ::testing::Test {
    protected:
        void SetUp() override {
            gmsg = new Inform("", output_m, INFORM_ALL_NODES);
            Configure::configure();
        }

        void TearDown() override {
            OpalData::deleteInstance();
            delete gmsg;
            gmsg = nullptr;
        }

        void run(const std::string& text) {
            std::string line = text;
            for (char& c : line) {
                if (c == '\n') c = ' ';
            }
            OpalParser parser;
            parser.run(new StringStream(line));
        }

        /// Parse and execute one statement; a ParseError is thrown, not reported (run() aborts).
        static void parseOne(const std::string& text) {
            OpalParser parser;
            StringStream stream(text);
            std::unique_ptr<Statement> statement(parser.readStatement(&stream));
            statement->execute(parser);
        }

        static double real(const std::string& object, const std::string& attribute) {
            Object* found = OpalData::getInstance()->find(object);
            if (found == nullptr) throw std::runtime_error("no object " + object);
            return Attributes::getReal(*found->findAttribute(attribute));
        }

        std::ostringstream output_m;
    };

}  // namespace

TEST_F(KeywordTableTest, KeywordsAreNotObjects) {
    OpalData* opal = OpalData::getInstance();
    EXPECT_TRUE(opal->isKeyword("QUADRUPOLE"));
    EXPECT_EQ(opal->find("QUADRUPOLE"), nullptr);
    EXPECT_FALSE(opal->isKeyword("Q1"));
}

TEST_F(KeywordTableTest, MakeObjectGivesANewObjectWithTheKeywordDefaults) {
    OpalData* opal = OpalData::getInstance();
    std::unique_ptr<Object> quad(opal->makeObject("QUADRUPOLE", "Q9"));
    ASSERT_NE(quad, nullptr);
    EXPECT_EQ(quad->getOpalName(), "Q9");
    EXPECT_EQ(quad->getKeyword(), "QUADRUPOLE");
    EXPECT_FALSE(quad->isBuiltin());
    EXPECT_EQ(opal->find("Q9"), nullptr);  // not stored
    EXPECT_EQ(opal->makeObject("NOT_A_KEYWORD", "X"), nullptr);
}

TEST_F(KeywordTableTest, CopyOfAnObjectKeepsItsValuesUntilEitherIsSetAgain) {
    run(R"(
        Q1: QUADRUPOLE, L = 0.1, K1 = 1.5, ELEMEDGE = 0;
        Q2: Q1, K1 = 2.5;
        Q1->L = 0.4;
    )");
    EXPECT_DOUBLE_EQ(real("Q1", "K1"), 1.5);
    EXPECT_DOUBLE_EQ(real("Q2", "K1"), 2.5);
    EXPECT_DOUBLE_EQ(real("Q1", "L"), 0.4);
    EXPECT_DOUBLE_EQ(real("Q2", "L"), 0.1);  // taken from Q1 when Q2 was defined
    EXPECT_EQ(OpalData::getInstance()->find("Q2")->getKeyword(), "QUADRUPOLE");
}

TEST_F(KeywordTableTest, InputObjectsCanBeRedefinedButKeywordsCannot) {
    run("D1: DRIFT, L = 0.5, ELEMEDGE = 0; D1: DRIFT, L = 0.7, ELEMEDGE = 0;");
    EXPECT_DOUBLE_EQ(real("D1", "L"), 0.7);
    EXPECT_THROW(run("QUADRUPOLE: DRIFT, L = 1, ELEMEDGE = 0;"), OpalException);
    EXPECT_EQ(OpalData::getInstance()->find("QUADRUPOLE"), nullptr);
}

TEST_F(KeywordTableTest, SharedIsAccepted) {
    run("SHARED D2: DRIFT, L = 0.3, ELEMEDGE = 0; SHARED D2, L = 0.9;");
    EXPECT_DOUBLE_EQ(real("D2", "L"), 0.3);  // "SHARED NAME, ...;" has no effect
}

TEST_F(KeywordTableTest, LabelledActionRunsAgainWithItsOwnValues) {
    run("O1: OPTION, PSDUMPFREQ = 5;");
    EXPECT_EQ(Options::psDumpFreq, 5);
    run("OPTION, PSDUMPFREQ = 9;");
    EXPECT_EQ(Options::psDumpFreq, 9);
    run("O1;");
    EXPECT_EQ(Options::psDumpFreq, 5);
}

TEST_F(KeywordTableTest, SecondOptionKeepsTheValuesOfTheFirst) {
    run("OPTION, PSDUMPFREQ = 7; OPTION, STATDUMPFREQ = 3;");
    EXPECT_EQ(Options::psDumpFreq, 7);
    EXPECT_EQ(Options::statDumpFreq, 3);
}

TEST_F(KeywordTableTest, KeywordNameIsAcceptedAsString) {
    run("STRING S1 = RING;");
    Object* s1 = OpalData::getInstance()->find("S1");
    ASSERT_NE(s1, nullptr);
    EXPECT_EQ(Attributes::getString(*s1->findAttribute("VALUE")), "RING");
}

TEST_F(KeywordTableTest, KeywordsCannotBeUsedAsObjects) {
    EXPECT_THROW(parseOne("QUADRUPOLE->L = 1;"), ParseError);
    EXPECT_THROW(parseOne("L1: LINE = (DRIFT);"), ParseError);
}

TEST_F(KeywordTableTest, MacroIsCalledOnlyWithItsArguments) {
    run("M(): MACRO { };");
    EXPECT_THROW(parseOne("M;"), ParseError);
    EXPECT_THROW(parseOne("X: M;"), ParseError);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ippl::initialize(argc, argv);
    const int result = RUN_ALL_TESTS();
    endfieldmodel::EndFieldModelManager::clearEFMManager();
    ippl::finalize();
    return result;
}
