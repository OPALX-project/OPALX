/**
 * \file TestRandom.cpp
 * \brief Unit tests for random number generation (replacement for gsl_rng)
 *
 * This file contains comprehensive unit tests for the GSL-compatible random
 * number generator implementation. The implementation uses \c std::mt19937_64
 * (Mersenne Twister) as the underlying engine and provides uniform and
 * Gaussian distributions.
 *
 * \test RandomTest::UniformDistribution
 * Tests uniform random number generation in the range [0, 1). Generates 1000
 * samples and verifies: (1) all values are within bounds, and (2) the mean
 * is approximately 0.5, confirming uniform distribution properties.
 *
 * \test RandomTest::GaussianDistribution
 * Tests Gaussian (normal) random number generation. Generates 10,000 samples
 * with \f$\sigma = 2.0\f$ and verifies: (1) the mean is approximately zero,
 * and (2) the variance is approximately \f$\sigma^2 = 4.0\f$.
 *
 * \test RandomTest::SeedReproducibility
 * Tests that seeding the RNG with the same value produces identical sequences.
 * Creates two RNG instances with the same seed and verifies that they generate
 * the same sequence of random numbers, ensuring deterministic behavior.
 *
 * \test RandomTest::DifferentSeedsDifferentSequences
 * Tests that different seeds produce different random sequences. Creates two
 * RNG instances with different seeds and verifies that they generate different
 * sequences, confirming proper seed handling.
 *
 * \test RandomTest::GaussianStatistics
 * Tests statistical properties of Gaussian distribution with large sample size
 * (100,000 samples). Verifies that the computed mean and standard deviation
 * match the expected values (\f$\mu = 0\f$, \f$\sigma = 1.0\f$) within
 * statistical tolerance.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include "Utilities/Random.h"

class RandomTest : public ::testing::Test {
protected:
    void SetUp() override {
        rng = gsl_rng_alloc(gsl_rng_default);
        rng->seed(12345);  // Fixed seed for reproducibility
    }

    void TearDown() override { gsl_rng_free(rng); }

    gsl_rng* rng;
};

TEST_F(RandomTest, UniformDistribution) {
    std::vector<double> samples(1000);
    for (size_t i = 0; i < samples.size(); ++i) {
        samples[i] = gsl_rng_uniform(rng);
    }

    // Check bounds
    for (double val : samples) {
        EXPECT_GE(val, 0.0);
        EXPECT_LT(val, 1.0);
    }

    // Check mean (should be around 0.5)
    double mean = std::accumulate(samples.begin(), samples.end(), 0.0) / samples.size();
    EXPECT_NEAR(mean, 0.5, 0.1);
}

TEST_F(RandomTest, GaussianDistribution) {
    double sigma = 2.0;
    std::vector<double> samples(10000);
    for (size_t i = 0; i < samples.size(); ++i) {
        samples[i] = gsl_ran_gaussian(rng, sigma);
    }

    // Check mean (should be around 0)
    double mean = std::accumulate(samples.begin(), samples.end(), 0.0) / samples.size();
    EXPECT_NEAR(mean, 0.0, 0.1);

    // Check variance (should be around sigma^2)
    double variance = 0.0;
    for (double val : samples) {
        variance += val * val;
    }
    variance /= samples.size();
    EXPECT_NEAR(variance, sigma * sigma, 0.5);
}

TEST_F(RandomTest, SeedReproducibility) {
    gsl_rng* rng1 = gsl_rng_alloc(gsl_rng_default);
    gsl_rng* rng2 = gsl_rng_alloc(gsl_rng_default);

    rng1->seed(42);
    rng2->seed(42);

    // Generate same sequence
    for (int i = 0; i < 100; ++i) {
        double val1 = gsl_rng_uniform(rng1);
        double val2 = gsl_rng_uniform(rng2);
        EXPECT_DOUBLE_EQ(val1, val2);
    }

    gsl_rng_free(rng1);
    gsl_rng_free(rng2);
}

TEST_F(RandomTest, DifferentSeedsDifferentSequences) {
    gsl_rng* rng1 = gsl_rng_alloc(gsl_rng_default);
    gsl_rng* rng2 = gsl_rng_alloc(gsl_rng_default);

    rng1->seed(100);
    rng2->seed(200);

    // Generate sequences - should be different
    bool different = false;
    for (int i = 0; i < 100; ++i) {
        double val1 = gsl_rng_uniform(rng1);
        double val2 = gsl_rng_uniform(rng2);
        if (std::abs(val1 - val2) > 1e-10) {
            different = true;
            break;
        }
    }
    EXPECT_TRUE(different);

    gsl_rng_free(rng1);
    gsl_rng_free(rng2);
}

TEST_F(RandomTest, GaussianStatistics) {
    double sigma  = 1.0;
    int n_samples = 100000;
    std::vector<double> samples(n_samples);

    for (int i = 0; i < n_samples; ++i) {
        samples[i] = gsl_ran_gaussian(rng, sigma);
    }

    // Calculate mean
    double mean = std::accumulate(samples.begin(), samples.end(), 0.0) / n_samples;

    // Calculate standard deviation
    double variance = 0.0;
    for (double val : samples) {
        double diff = val - mean;
        variance += diff * diff;
    }
    variance /= n_samples;
    double std_dev = std::sqrt(variance);

    EXPECT_NEAR(mean, 0.0, 0.01);
    EXPECT_NEAR(std_dev, sigma, 0.01);
}
