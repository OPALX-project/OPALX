#include "Util/OpalInputFileParser.h"
#include "Util/Types.h"
#include "gtest/gtest.h"

#include "Expression/Parser/function.hpp"

#include "boost/smart_ptr.hpp"
#include "boost/tuple/tuple.hpp"

namespace {

    // The fixture for testing class Foo.
    class InputFileParserTest : public ::testing::Test {
    protected:

        InputFileParserTest() {
            // You can do set-up work for each test here.
            functionDictionary_t funcs;
            opalp.reset(new OpalInputFileParser("resources/test.in", funcs));
        }

        virtual ~InputFileParserTest() {
            // You can do clean-up work that doesn't throw exceptions here.
            //delete opalp;
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:

        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
            opalp->doParse();
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        // Objects declared here can be used by all tests in the test case
        boost::scoped_ptr<OpalInputFileParser> opalp;
    };

    TEST_F(InputFileParserTest, ProblemSize) {

        Expressions::Named_t obj, constr;
        DVarContainer_t dvars;
        opalp->getProblem(obj, constr, dvars);

        EXPECT_EQ(static_cast<size_t>(3), obj.size());
        EXPECT_EQ(static_cast<size_t>(1), constr.size());
        EXPECT_EQ(static_cast<size_t>(1), dvars.size());
    }

    TEST_F(InputFileParserTest, ParsedCorrectDesignVars) {

        Expressions::Named_t obj, constr;
        DVarContainer_t dvars;
        opalp->getProblem(obj, constr, dvars);

        DVarContainer_t::iterator itr;
        for(itr = dvars.begin(); itr != dvars.end(); itr++) {
            EXPECT_EQ("d1", itr->first);
            EXPECT_EQ("KS", boost::get<VAR_NAME>(itr->second));
            EXPECT_EQ(0.0,  boost::get<LOWER_BOUND>(itr->second));
            EXPECT_EQ(0.5,  boost::get<UPPER_BOUND>(itr->second));
        }
    }

    TEST_F(InputFileParserTest, ParsedCorrectExpressionObjectives) {

        Expressions::Named_t obj, constr;
        DVarContainer_t dvars;
        opalp->getProblem(obj, constr, dvars);

        std::vector<std::string> names;
        names.push_back("obj1");
        names.push_back("ob3");
        names.push_back("o2");
        for(Expressions::SingleNamed_t expr : obj) {
            std::string expected_name = names.back();
            EXPECT_EQ(expected_name, expr.first);
            names.pop_back();
        }
    }

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
