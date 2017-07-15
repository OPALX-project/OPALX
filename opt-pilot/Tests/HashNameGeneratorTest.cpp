#include "Util/HashNameGenerator.h"
#include "gtest/gtest.h"

namespace {

    // The fixture for testing class Foo.
    class HashNameGeneratorTest : public ::testing::Test {
    protected:

        HashNameGeneratorTest() {
            // You can do set-up work for each test here.
        }

        virtual ~HashNameGeneratorTest() {
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

    TEST_F(HashNameGeneratorTest, HashName) {

        std::vector<std::string> params;
        params.push_back("transferline");
        params.push_back("sigmax=5.05");
        params.push_back("sigmay=6.05");

        std::string hash = HashNameGenerator::generate(params);

        std::string expected = "2d6ea8be68514560d38cf107709163668906c8a9";
        ASSERT_STREQ(expected.c_str(), hash.c_str());
    }

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
