/**
 * \file TestHistogram.cpp
 * \brief Unit tests for histogram implementation (replacement for gsl_histogram)
 *
 * This file contains comprehensive unit tests for 1D and 2D histogram
 * implementations. The histograms support uniform and custom bin ranges,
 * increment operations, and value retrieval.
 *
 * \test HistogramTest::UniformRanges
 * Tests creation of histogram with uniform bin ranges. Verifies that bins are
 * correctly allocated and that range boundaries are set correctly for uniform
 * spacing between specified minimum and maximum values.
 *
 * \test HistogramTest::IncrementAndGet
 * Tests incrementing histogram bins and retrieving bin values. Verifies that
 * values are correctly assigned to appropriate bins based on their position
 * within the range, and that bin counts can be retrieved accurately.
 *
 * \test HistogramTest::CustomRanges
 * Tests histogram with custom (non-uniform) bin ranges. Verifies that bins
 * with arbitrary range boundaries work correctly and that values are assigned
 * to the correct bins based on the custom ranges.
 *
 * \test HistogramTest::OutOfRange
 * Tests handling of values outside the histogram range. Verifies that values
 * below the minimum or above the maximum range do not increment any bins and
 * are silently ignored (no exceptions thrown).
 *
 * \test HistogramTest::Histogram2D
 * Tests 2D histogram functionality with uniform ranges. Verifies that 2D
 * histograms can be created, values can be incremented using (x, y) coordinates,
 * and bin values can be retrieved correctly.
 *
 * \test HistogramTest::Histogram2DCustomRanges
 * Tests 2D histogram with custom ranges for both x and y dimensions. Verifies
 * that non-uniform bin boundaries work correctly in both dimensions and that
 * values are assigned to the correct 2D bins.
 *
 * \test HistogramTest::Histogram2DOutOfRange
 * Tests handling of out-of-range values in 2D histograms. Verifies that values
 * outside the x or y ranges (or both) do not increment any bins.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "Utilities/GSLHistogram.h"

class HistogramTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Test setup
    }
};

TEST_F(HistogramTest, UniformRanges) {
    gsl_histogram* h = gsl_histogram_alloc(10);
    gsl_histogram_set_ranges_uniform(h, 0.0, 10.0);

    EXPECT_EQ(gsl_histogram_bins(h), 10);

    // Check range boundaries
    double lower, upper;
    gsl_histogram_get_range(h, 0, &lower, &upper);
    EXPECT_DOUBLE_EQ(lower, 0.0);
    EXPECT_DOUBLE_EQ(upper, 1.0);

    gsl_histogram_get_range(h, 9, &lower, &upper);
    EXPECT_DOUBLE_EQ(lower, 9.0);
    EXPECT_DOUBLE_EQ(upper, 10.0);

    gsl_histogram_free(h);
}

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

TEST_F(HistogramTest, CustomRanges) {
    gsl_histogram* h           = gsl_histogram_alloc(5);
    std::vector<double> ranges = {0.0, 1.0, 2.0, 5.0, 10.0, 20.0};
    gsl_histogram_set_ranges(h, ranges.data(), ranges.size());

    gsl_histogram_increment(h, 1.5);
    gsl_histogram_increment(h, 3.0);
    gsl_histogram_increment(h, 7.5);

    EXPECT_DOUBLE_EQ(gsl_histogram_get(h, 1), 1.0);  // 1.5 in bin 1
    EXPECT_DOUBLE_EQ(gsl_histogram_get(h, 2), 1.0);  // 3.0 in bin 2
    EXPECT_DOUBLE_EQ(gsl_histogram_get(h, 3), 1.0);  // 7.5 in bin 3

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

TEST_F(HistogramTest, Histogram2D) {
    gsl_histogram2d* h = gsl_histogram2d_alloc(5, 5);
    gsl_histogram2d_set_ranges_uniform(h, 0.0, 5.0, 0.0, 5.0);

    EXPECT_EQ(gsl_histogram2d_nx(h), 5);
    EXPECT_EQ(gsl_histogram2d_ny(h), 5);

    // Increment some values
    gsl_histogram2d_increment(h, 1.5, 2.5);
    gsl_histogram2d_increment(h, 1.7, 2.3);

    EXPECT_DOUBLE_EQ(gsl_histogram2d_get(h, 1, 2), 2.0);

    gsl_histogram2d_free(h);
}

TEST_F(HistogramTest, Histogram2DCustomRanges) {
    gsl_histogram2d* h = gsl_histogram2d_alloc(3, 3);

    std::vector<double> xranges = {0.0, 1.0, 3.0, 6.0};
    std::vector<double> yranges = {0.0, 2.0, 4.0, 8.0};

    gsl_histogram2d_set_ranges(h, xranges.data(), xranges.size(), yranges.data(), yranges.size());

    gsl_histogram2d_increment(h, 2.0, 3.0);
    gsl_histogram2d_increment(h, 4.0, 5.0);

    EXPECT_DOUBLE_EQ(gsl_histogram2d_get(h, 1, 1), 1.0);
    EXPECT_DOUBLE_EQ(gsl_histogram2d_get(h, 2, 2), 1.0);

    gsl_histogram2d_free(h);
}

TEST_F(HistogramTest, Histogram2DOutOfRange) {
    gsl_histogram2d* h = gsl_histogram2d_alloc(5, 5);
    gsl_histogram2d_set_ranges_uniform(h, 0.0, 5.0, 0.0, 5.0);

    // Values outside range
    gsl_histogram2d_increment(h, -1.0, 2.0);
    gsl_histogram2d_increment(h, 2.0, -1.0);
    gsl_histogram2d_increment(h, 6.0, 2.0);
    gsl_histogram2d_increment(h, 2.0, 6.0);

    // All bins should be zero
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            EXPECT_DOUBLE_EQ(gsl_histogram2d_get(h, i, j), 0.0);
        }
    }

    gsl_histogram2d_free(h);
}
