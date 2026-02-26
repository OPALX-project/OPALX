/**
 * \file TestGSLIntegration.cpp
 * \brief Unit tests for GSL numerical integration routines
 *
 * This file contains unit tests for GSL-compatible numerical integration
 * routines. The implementation uses adaptive Simpson's rule for computing
 * definite integrals with specified error tolerances.
 *
 * \test GSLIntegrationTest::LinearFunction
 * Tests integration of a linear function \f$f(x) = x\f$ over the interval
 * \f$[0, 1]\f$. Verifies that \f$\int_0^1 x \, dx = \frac{1}{2}\f$ is computed
 * correctly within the specified tolerance.
 *
 * \test GSLIntegrationTest::QuadraticFunction
 * Tests integration of a quadratic function \f$f(x) = x^2\f$ over the interval
 * \f$[0, 2]\f$. Verifies that \f$\int_0^2 x^2 \, dx = \frac{8}{3}\f$ is computed
 * correctly.
 *
 * \test GSLIntegrationTest::SinFunction
 * Tests integration of the sine function \f$f(x) = \sin(x)\f$ over the interval
 * \f$[0, \pi]\f$. Verifies that \f$\int_0^{\pi} \sin(x) \, dx = 2\f$ is computed
 * correctly, confirming proper handling of oscillatory functions.
 *
 * \test GSLIntegrationTest::ExpFunction
 * Tests integration of the exponential function \f$f(x) = e^x\f$ over the interval
 * \f$[0, 1]\f$. Verifies that \f$\int_0^1 e^x \, dx = e - 1\f$ is computed
 * correctly.
 *
 * \test GSLIntegrationTest::NegativeRange
 * Tests integration over a range that includes negative values. Verifies that
 * \f$\int_{-1}^1 x \, dx = 0\f$ (odd function over symmetric interval) is
 * computed correctly.
 *
 * \test GSLIntegrationTest::SmallRange
 * Tests integration over a small interval. Verifies that the adaptive algorithm
 * works correctly for narrow integration ranges and maintains accuracy.
 */

#include <gtest/gtest.h>
#include <cmath>
#include "Utilities/GSLIntegration.h"

class GSLIntegrationTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Test setup
    }
};

// Test function: f(x) = x
double linear_func(double x, void* /* params */) { return x; }

// Test function: f(x) = x^2
double quadratic_func(double x, void* /* params */) { return x * x; }

// Test function: f(x) = sin(x)
double sin_func(double x, void* /* params */) { return std::sin(x); }

// Test function: f(x) = exp(x)
double exp_func(double x, void* /* params */) { return std::exp(x); }

TEST_F(GSLIntegrationTest, LinearFunction) {
    gsl_function F;
    F.function = linear_func;
    F.params   = nullptr;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

    double result, error;
    int status = gsl_integration_qag(
        &F, 0.0, 1.0, 1e-10, 1e-10, 1000, GSL_INTEG_GAUSS15, w, &result, &error
    );

    EXPECT_EQ(status, 0);
    EXPECT_NEAR(result, 0.5, 1e-6);  // ∫x dx from 0 to 1 = 0.5

    gsl_integration_workspace_free(w);
}

TEST_F(GSLIntegrationTest, QuadraticFunction) {
    gsl_function F;
    F.function = quadratic_func;
    F.params   = nullptr;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

    double result, error;
    int status = gsl_integration_qag(
        &F, 0.0, 2.0, 1e-10, 1e-10, 1000, GSL_INTEG_GAUSS15, w, &result, &error
    );

    EXPECT_EQ(status, 0);
    EXPECT_NEAR(result, 8.0 / 3.0, 1e-6);  // ∫x^2 dx from 0 to 2 = 8/3

    gsl_integration_workspace_free(w);
}

TEST_F(GSLIntegrationTest, SinFunction) {
    gsl_function F;
    F.function = sin_func;
    F.params   = nullptr;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

    double result, error;
    int status = gsl_integration_qag(
        &F, 0.0, M_PI, 1e-10, 1e-10, 1000, GSL_INTEG_GAUSS15, w, &result, &error
    );

    EXPECT_EQ(status, 0);
    EXPECT_NEAR(result, 2.0, 1e-6);  // ∫sin(x) dx from 0 to π = 2

    gsl_integration_workspace_free(w);
}

TEST_F(GSLIntegrationTest, ExpFunction) {
    gsl_function F;
    F.function = exp_func;
    F.params   = nullptr;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

    double result, error;
    int status = gsl_integration_qag(
        &F, 0.0, 1.0, 1e-10, 1e-10, 1000, GSL_INTEG_GAUSS15, w, &result, &error
    );

    EXPECT_EQ(status, 0);
    EXPECT_NEAR(result, std::exp(1.0) - 1.0, 1e-6);  // ∫e^x dx from 0 to 1 = e - 1

    gsl_integration_workspace_free(w);
}

TEST_F(GSLIntegrationTest, NegativeRange) {
    gsl_function F;
    F.function = linear_func;
    F.params   = nullptr;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

    double result, error;
    int status = gsl_integration_qag(
        &F, -1.0, 1.0, 1e-10, 1e-10, 1000, GSL_INTEG_GAUSS15, w, &result, &error
    );

    EXPECT_EQ(status, 0);
    EXPECT_NEAR(result, 0.0, 1e-6);  // ∫x dx from -1 to 1 = 0

    gsl_integration_workspace_free(w);
}

TEST_F(GSLIntegrationTest, SmallRange) {
    gsl_function F;
    F.function = quadratic_func;
    F.params   = nullptr;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

    double result, error;
    int status = gsl_integration_qag(
        &F, 0.0, 0.1, 1e-10, 1e-10, 1000, GSL_INTEG_GAUSS15, w, &result, &error
    );

    EXPECT_EQ(status, 0);
    EXPECT_NEAR(result, 0.1 * 0.1 * 0.1 / 3.0, 1e-6);

    gsl_integration_workspace_free(w);
}
