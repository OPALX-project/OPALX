/**
 * \file TestCubicSpline.cpp
 * \brief Unit tests for CubicSpline interpolation (replacement for gsl_spline)
 * 
 * This file contains comprehensive unit tests for the CubicSpline class, which
 * implements cubic spline interpolation using natural boundary conditions.
 * The implementation provides both a direct C++ interface and GSL-compatible
 * wrapper functions.
 * 
 * \test CubicSplineTest::BasicInterpolation
 * Tests basic cubic spline interpolation at data points and between points.
 * Verifies that the spline exactly reproduces input data points and provides
 * smooth interpolation between them. Uses test data \f$y = x^2\f$.
 * 
 * \test CubicSplineTest::Extrapolation
 * Tests extrapolation behavior beyond the data range. Verifies that the spline
 * can extrapolate values to the left (before first data point) and right
 * (after last data point) of the input data range.
 * 
 * \test CubicSplineTest::Accelerator
 * Tests the accelerator optimization for repeated evaluations. The accelerator
 * caches the last interval index to speed up subsequent evaluations in the
 * same interval. Verifies correct behavior when evaluating in the same interval
 * and when switching to different intervals.
 * 
 * \test CubicSplineTest::SinFunction
 * Tests interpolation accuracy on a smooth function (sine). Creates a spline
 * from sampled sine function data and verifies that interpolated values match
 * the true function values within acceptable tolerance.
 * 
 * \test CubicSplineTest::InvalidInput
 * Tests error handling for invalid input data. Verifies that exceptions are
 * thrown when: (1) too few data points (< 2) are provided, and (2) x-values
 * are not strictly increasing.
 * 
 * \test CubicSplineTest::GSLCompatibleInterface
 * Tests the GSL-compatible wrapper functions (gsl_spline_alloc, gsl_spline_init,
 * gsl_spline_eval, etc.). Verifies that the GSL API is correctly implemented
 * and that the accelerator reset functionality works properly.
 */

#include <gtest/gtest.h>
#include "Utilities/CubicSpline.h"
#include <cmath>
#include <vector>

class CubicSplineTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Test data: y = x^2
        x_data = {0.0, 1.0, 2.0, 3.0, 4.0};
        y_data = {0.0, 1.0, 4.0, 9.0, 16.0};
    }
    
    std::vector<double> x_data;
    std::vector<double> y_data;
};

TEST_F(CubicSplineTest, BasicInterpolation) {
    CubicSpline spline(x_data.data(), y_data.data(), x_data.size());
    
    // Test interpolation at data points
    EXPECT_NEAR(spline.eval(0.0), 0.0, 1e-10);
    EXPECT_NEAR(spline.eval(1.0), 1.0, 1e-10);
    EXPECT_NEAR(spline.eval(2.0), 4.0, 1e-10);
    EXPECT_NEAR(spline.eval(3.0), 9.0, 1e-10);
    EXPECT_NEAR(spline.eval(4.0), 16.0, 1e-10);
    
    // Test interpolation between points
    double val = spline.eval(1.5);
    EXPECT_GT(val, 1.0);
    EXPECT_LT(val, 4.0);
}

TEST_F(CubicSplineTest, Extrapolation) {
    CubicSpline spline(x_data.data(), y_data.data(), x_data.size());
    
    // Test left extrapolation
    double val_left = spline.eval(-1.0);
    EXPECT_LT(val_left, 0.0);
    
    // Test right extrapolation
    double val_right = spline.eval(5.0);
    EXPECT_GT(val_right, 16.0);
}

TEST_F(CubicSplineTest, Accelerator) {
    CubicSpline spline(x_data.data(), y_data.data(), x_data.size());
    CubicSpline::Accelerator accel;
    
    // First evaluation
    double val1 = spline.eval(1.5, accel);
    EXPECT_GT(val1, 1.0);
    EXPECT_LT(val1, 4.0);
    
    // Second evaluation in same interval (should use cached index)
    double val2 = spline.eval(1.6, accel);
    EXPECT_GT(val2, val1);
    
    // Evaluation in different interval
    double val3 = spline.eval(2.5, accel);
    EXPECT_GT(val3, 4.0);
    EXPECT_LT(val3, 9.0);
}

TEST_F(CubicSplineTest, SinFunction) {
    std::vector<double> x(10);
    std::vector<double> y(10);
    for (size_t i = 0; i < 10; ++i) {
        x[i] = i * 0.5;
        y[i] = std::sin(x[i]);
    }
    
    CubicSpline spline(x.data(), y.data(), x.size());
    
    // Test interpolation
    double val = spline.eval(1.25);
    double expected = std::sin(1.25);
    EXPECT_NEAR(val, expected, 0.01);
}

TEST_F(CubicSplineTest, InvalidInput) {
    // Too few points
    std::vector<double> x = {0.0};
    std::vector<double> y = {0.0};
    EXPECT_THROW(CubicSpline spline(x.data(), y.data(), 1), std::invalid_argument);
    
    // Non-increasing x
    x = {0.0, 1.0, 0.5, 2.0};
    y = {0.0, 1.0, 2.0, 3.0};
    EXPECT_THROW(CubicSpline spline(x.data(), y.data(), 4), std::invalid_argument);
}

TEST_F(CubicSplineTest, GSLCompatibleInterface) {
    gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, x_data.size());
    gsl_interp_accel* accel = gsl_interp_accel_alloc();
    
    gsl_spline_init(spline, x_data.data(), y_data.data(), x_data.size());
    
    // Test evaluation
    double val = gsl_spline_eval(spline, 1.5, accel);
    EXPECT_GT(val, 1.0);
    EXPECT_LT(val, 4.0);
    
    // Test at boundary
    val = gsl_spline_eval(spline, 0.0, accel);
    EXPECT_NEAR(val, 0.0, 1e-10);
    
    // Reset accelerator
    gsl_interp_accel_reset(accel);
    val = gsl_spline_eval(spline, 2.0, accel);
    EXPECT_NEAR(val, 4.0, 1e-10);
    
    gsl_spline_free(spline);
    gsl_interp_accel_free(accel);
}

