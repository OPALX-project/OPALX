#include "gtest/gtest.h"

#include "Sample/SampleGaussianSequence.h"
#include "opal_test_utilities/SilenceTest.h"

TEST(GaussSampleTest, ChainTest) {
    OpalTestUtilities::SilenceTest silencer;
    unsigned int nSample = 101;
    SampleGaussianSequence seq(-5, 5, nSample);

    for (unsigned int i = 0; i * 2 < nSample - 1; ++ i) {
        seq.getNext();
    }
    double x = seq.getNext();

    EXPECT_NEAR(x, 0.0, 1e-8);
}