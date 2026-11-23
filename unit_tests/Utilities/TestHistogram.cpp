/**
 * \file TestHistogram.cpp
 * \brief Unit tests for histogram implementation (replacement for gsl_histogram)
 *
 * This file contains comprehensive unit tests for the 1D histogram
 * implementation. The histogram supports uniform bin ranges, increment
 * operations, and value retrieval.
 *
 * \test HistogramTest::IncrementAndGet
 * Tests incrementing histogram bins and retrieving bin values. Verifies that
 * values are correctly assigned to appropriate bins based on their position
 * within the range, and that bin counts can be retrieved accurately.
 *
 * \test HistogramTest::OutOfRange
 * Tests handling of values outside the histogram range. Verifies that values
 * below the minimum or above the maximum range do not increment any bins and
 * are silently ignored (no exceptions thrown).
 */

#include <gtest/gtest.h>
#include <cmath>
#include "Utilities/GSLHistogram.h"

class HistogramTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Test setup
    }
};

TEST_F(HistogramTest, IncrementAndGet) {
    gsl_histogram* h = gsl_histogram_alloc(10);
    gsl_histogram_set_ranges_uniform(h, 0.0, 10.0);

    // Increment some values
    gsl_histogram_increment(h, 1.5);
    gsl_histogram_increment(h, 2.5);
    gsl_histogram_increment(h, 1.7);

    EXPECT_DOUBLE_EQ(gsl_histogram_get(h, 1), 2.0);  // 1.5 and 1.7 in bin 1
    EXPECT_DOUBLE_EQ(gsl_histogram_get(h, 2), 1.0);  // 2.5 in bin 2

    gsl_histogram_free(h);
}

TEST_F(HistogramTest, OutOfRange) {
    gsl_histogram* h = gsl_histogram_alloc(10);
    gsl_histogram_set_ranges_uniform(h, 0.0, 10.0);

    // Values outside range should not increment any bin
    gsl_histogram_increment(h, -1.0);
    gsl_histogram_increment(h, 11.0);

    // All bins should be zero
    for (size_t i = 0; i < 10; ++i) {
        EXPECT_DOUBLE_EQ(gsl_histogram_get(h, i), 0.0);
    }

    gsl_histogram_free(h);
}
