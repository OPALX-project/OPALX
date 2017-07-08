#include <set>
#include <string>

#include "Util/OptPilotException.h"
#include "Expression/Expression.h"
#include "Expression/Parser/function.hpp"
#include "Expression/FromFile.h"

#include "gtest/gtest.h"

#include "boost/smart_ptr.hpp"
#include "boost/tuple/tuple.hpp"
#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"


namespace {

    // The fixture for testing class Foo.
    class FromFileExpressionTest : public ::testing::Test {
    protected:

        FromFileExpressionTest() {
            // You can do set-up work for each test here.
        }

        virtual ~FromFileExpressionTest() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:

        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }
    };


    TEST_F(FromFileExpressionTest, EvaluateFromFileExpression) {

        variableDictionary_t vars;
        double expected = 0.3126 + 0.3561 + 0.4242;

        functionDictionary_t funcs;
        client::function::type ff;
        ff = FromFile();
        funcs.insert(std::pair<std::string, client::function::type>
                ("fromFile", ff));

        std::string testexpr = "fromFile(\"resources/fromfile_test.dat\", \"var\")";
        boost::scoped_ptr<Expression> e(new Expression(testexpr, funcs));
        Expressions::Result_t result;
        EXPECT_NO_THROW({
            result = e->evaluate(vars);
        });

        ASSERT_NEAR(expected, boost::get<0>(result), 1e-6);
        ASSERT_TRUE(boost::get<1>(result));
    }

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
