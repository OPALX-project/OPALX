
#include "Util/SDDSReader.h"
#include "gtest/gtest.h"


namespace {

    // The fixture for testing class Foo.
    class SDDSParserTest : public ::testing::Test {
    protected:

        SDDSParserTest() {
            // You can do set-up work for each test here.
            sddsr = new SDDSReader("resources/test.stat");
        }

        virtual ~SDDSParserTest() {
            // You can do clean-up work that doesn't throw exceptions here.
            delete sddsr;
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:

        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
            sddsr->parseFile();
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        // Objects declared here can be used by all tests in the test case
        SDDSReader *sddsr;
    };

    TEST_F(SDDSParserTest, ReadEnergy) {

        double energy = 0.0;
        std::string s = "energy";
        sddsr->getValue(1, s, energy);

        double expected = 2.220252349855340e-01;
        ASSERT_DOUBLE_EQ(expected, energy);
    }

    TEST_F(SDDSParserTest, ReadLastPosition) {

        double position = 0.0;
        std::string s = "s";
        sddsr->getValue(-1, s, position);

        double expected = 4.247573354842603e-03;
        ASSERT_DOUBLE_EQ(expected, position);
    }

    TEST_F(SDDSParserTest, InterpolateRms_x) {

        double spos = 4.0e-03;
        double rmsx_interp = 0.0;

        EXPECT_NO_THROW({
            sddsr->getInterpolatedValue(spos, "rms_x", rmsx_interp);
        });

        double spos_before = 3.786226707177705e-03;
        double spos_after  = 4.015012129060398e-03;
        double rmsx_before = 3.147090549966750e-04;
        double rmsx_after  = 3.166390884805550e-04;

        double expected = rmsx_before + (spos - spos_before) * (rmsx_after - rmsx_before) / (spos_after - spos_before);

        ASSERT_DOUBLE_EQ(expected, rmsx_interp);
    }

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

