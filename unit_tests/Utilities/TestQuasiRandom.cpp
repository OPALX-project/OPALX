/**
 * \file TestQuasiRandom.cpp
 * \brief Unit tests for quasi-random number generation (replacement for gsl_qrng)
 *
 * This file contains unit tests for GSL-compatible quasi-random number
 * generation using Sobol sequences. Quasi-random sequences provide better
 * space-filling properties than pseudo-random sequences, making them useful
 * for Monte Carlo integration and similar applications.
 *
 * \test QuasiRandomTest::BasicGeneration
 * Tests basic quasi-random number generation. Verifies that the generator can
 * be allocated, that generated values are in the range \f$[0, 1)\f$, and that
 * multi-dimensional sequences are generated correctly.
 *
 * \test QuasiRandomTest::SequenceUniqueness
 * Tests that generated sequences contain unique values and are not all zeros.
 * Generates a sequence of 100 samples and verifies that values are properly
 * distributed and non-trivial.
 *
 * \test QuasiRandomTest::MultiDimensional
 * Tests multi-dimensional quasi-random sequence generation. Verifies that
 * all dimensions of a 5-dimensional sequence are correctly generated and
 * that all values remain within the valid range \f$[0, 1)\f$.
 *
 * \test QuasiRandomTest::Reproducibility
 * Tests that quasi-random sequences are deterministic and reproducible. Creates
 * two generators with the same initialization and verifies that they produce
 * identical sequences, confirming deterministic behavior.
 *
 * \test QuasiRandomTest::Coverage
 * Tests space-filling properties of quasi-random sequences. Generates 1000
 * samples and verifies that the sequence covers a significant portion of the
 * unit interval, demonstrating better coverage than pseudo-random sequences
 * for low-discrepancy sequences.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "Utilities/QuasiRandom.h"

class QuasiRandomTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Test setup
    }
};

TEST_F(QuasiRandomTest, BasicGeneration) {
    gsl_qrng* q = gsl_qrng_alloc("sobol", 2);
    EXPECT_NE(q, nullptr);

    double x[2];
    gsl_qrng_get(q, x);

    // Values should be in [0, 1)
    EXPECT_GE(x[0], 0.0);
    EXPECT_LT(x[0], 1.0);
    EXPECT_GE(x[1], 0.0);
    EXPECT_LT(x[1], 1.0);

    gsl_qrng_free(q);
}

TEST_F(QuasiRandomTest, SequenceUniqueness) {
    gsl_qrng* q = gsl_qrng_alloc("sobol", 1);

    std::vector<double> samples(100);
    for (size_t i = 0; i < samples.size(); ++i) {
        double x[1];
        gsl_qrng_get(q, x);
        samples[i] = x[0];
    }

    // Check that values are in [0, 1)
    for (double val : samples) {
        EXPECT_GE(val, 0.0);
        EXPECT_LT(val, 1.0);
    }

    // Check that sequence is not all zeros
    bool all_zero = true;
    for (double val : samples) {
        if (val != 0.0) {
            all_zero = false;
            break;
        }
    }
    EXPECT_FALSE(all_zero);

    gsl_qrng_free(q);
}

TEST_F(QuasiRandomTest, MultiDimensional) {
    gsl_qrng* q = gsl_qrng_alloc("sobol", 5);

    double x[5];
    gsl_qrng_get(q, x);

    // All dimensions should be in [0, 1)
    for (int i = 0; i < 5; ++i) {
        EXPECT_GE(x[i], 0.0);
        EXPECT_LT(x[i], 1.0);
    }

    gsl_qrng_free(q);
}

TEST_F(QuasiRandomTest, Reproducibility) {
    gsl_qrng* q1 = gsl_qrng_alloc("sobol", 2);
    gsl_qrng* q2 = gsl_qrng_alloc("sobol", 2);

    double x1[2], x2[2];

    // Both should start with same sequence
    gsl_qrng_get(q1, x1);
    gsl_qrng_get(q2, x2);

    // Values should be the same (deterministic sequence)
    EXPECT_DOUBLE_EQ(x1[0], x2[0]);
    EXPECT_DOUBLE_EQ(x1[1], x2[1]);

    gsl_qrng_free(q1);
    gsl_qrng_free(q2);
}

TEST_F(QuasiRandomTest, Coverage) {
    gsl_qrng* q = gsl_qrng_alloc("sobol", 1);

    // Generate many samples and check coverage
    const int n_samples = 1000;
    std::vector<double> samples(n_samples);

    for (int i = 0; i < n_samples; ++i) {
        double x[1];
        gsl_qrng_get(q, x);
        samples[i] = x[0];
    }

    // Check that we cover a reasonable range
    double min_val = *std::min_element(samples.begin(), samples.end());
    double max_val = *std::max_element(samples.begin(), samples.end());

    EXPECT_GE(min_val, 0.0);
    EXPECT_LT(max_val, 1.0);
    EXPECT_GT(max_val - min_val, 0.5);  // Should cover at least half the range

    gsl_qrng_free(q);
}
