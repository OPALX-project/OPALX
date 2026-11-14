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
#include "AbstractObjects/ObjectFunction.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Ippl.h"
#include "OpalConfigure/Configure.h"
#include "OpalParser/OpalParser.h"
#include "OpalParser/StringStream.h"
#include "Utilities/OpalException.h"
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

        static double real(const std::string& object, const std::string& attribute) {
            Object* found = OpalData::getInstance()->find(object);
            if (found == nullptr) throw std::runtime_error("no object " + object);
            return Attributes::getReal(*found->findAttribute(attribute));
        }

        std::ostringstream output_m;
    };

}  // namespace

TEST_F(KeywordTableTest, KeywordsAreFoundButNotStoredWithTheInputObjects) {
    OpalData* opal  = OpalData::getInstance();
    Object* keyword = opal->findKeyword("QUADRUPOLE");
    ASSERT_NE(keyword, nullptr);
    EXPECT_EQ(opal->find("QUADRUPOLE"), keyword);

    // OpalData::apply() visits the stored objects only.
    bool visited = false;
    struct Visit : ObjectFunction {
        explicit Visit(bool& seen, const Object* keyword) : seen_m(seen), keyword_m(keyword) {}
        void operator()(Object* object) const override { seen_m |= object == keyword_m; }
        bool& seen_m;
        const Object* keyword_m;
    } visit(visited, keyword);
    opal->apply(visit);
    EXPECT_FALSE(visited);
}

TEST_F(KeywordTableTest, MakeObjectGivesANewObjectWithTheKeywordDefaults) {
    OpalData* opal = OpalData::getInstance();
    std::unique_ptr<Object> quad(opal->makeObject("QUADRUPOLE", "Q9"));
    ASSERT_NE(quad, nullptr);
    EXPECT_EQ(quad->getOpalName(), "Q9");
    EXPECT_FALSE(quad->isBuiltin());
    EXPECT_EQ(opal->find("Q9"), nullptr);  // not stored
    EXPECT_EQ(opal->makeObject("NOT_A_KEYWORD", "X"), nullptr);
}

TEST_F(KeywordTableTest, InputObjectsUseKeywordsAndOtherObjects) {
    run(R"(
        Q1: QUADRUPOLE, L = 0.1, K1 = 1.5, ELEMEDGE = 0;
        Q2: Q1, K1 = 2.5;
    )");
    EXPECT_DOUBLE_EQ(real("Q1", "K1"), 1.5);
    EXPECT_DOUBLE_EQ(real("Q2", "K1"), 2.5);
    EXPECT_DOUBLE_EQ(real("Q2", "L"), 0.1);  // taken from Q1
    EXPECT_EQ(OpalData::getInstance()->find("Q2")->getKeyword(), "QUADRUPOLE");
}

TEST_F(KeywordTableTest, InputObjectsCanBeRedefinedButKeywordsCannot) {
    run("D1: DRIFT, L = 0.5, ELEMEDGE = 0; D1: DRIFT, L = 0.7, ELEMEDGE = 0;");
    EXPECT_DOUBLE_EQ(real("D1", "L"), 0.7);
    EXPECT_THROW(run("QUADRUPOLE: DRIFT, L = 1, ELEMEDGE = 0;"), OpalException);
    EXPECT_EQ(
            OpalData::getInstance()->find("QUADRUPOLE"),
            OpalData::getInstance()->findKeyword("QUADRUPOLE"));
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ippl::initialize(argc, argv);
    const int result = RUN_ALL_TESTS();
    endfieldmodel::EndFieldModelManager::clearEFMManager();
    ippl::finalize();
    return result;
}
