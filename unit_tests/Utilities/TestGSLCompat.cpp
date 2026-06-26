/**
 * \file TestGSLCompat.cpp
 * \brief Unit tests for GSL special functions compatibility
 *
 * This file contains unit tests for GSL-compatible special function wrappers.
 * These functions provide compatibility with GSL's special function API while
 * using C++ standard library implementations where possible.
 *
 * \test GSLCompatTest::ErfFunction
 * Tests the error function \f$\text{erf}(x)\f$ implementation. Verifies that
 * \c gsl_sf_erf() matches \c std::erf() for various input values, including
 * zero, positive, negative, and large values approaching the asymptotic limit
 * of 1.
 *
 * \test GSLCompatTest::PowInt
 * Tests integer power function \f$x^n\f$ for both positive and negative
 * integer exponents. Verifies correct computation of powers including edge
 * cases: \f$x^0 = 1\f$, \f$x^1 = x\f$, and negative powers \f$x^{-n} = 1/x^n\f$.
 *
 * \test GSLCompatTest::GammaFunction
 * Tests the gamma function \f$\Gamma(x)\f$ implementation. Verifies that
 * \c gsl_sf_gamma() matches \c std::tgamma() for integer values (factorial
 * property: \f$\Gamma(n+1) = n!\f$) and half-integer values
 * (\f$\Gamma(1/2) = \sqrt{\pi}\f$).
 *
 * \test GSLCompatTest::Factorial
 * Tests factorial function \f$n!\f$ for small integer values. Verifies correct
 * computation of factorials from 0 to 10, including the special case \f$0! = 1\f$.
 *
 * \test GSLCompatTest::FactorialLarge
 * Tests factorial function for larger values (e.g., 20!). For large factorials,
 * the implementation uses the gamma function to avoid overflow. Verifies that
 * large factorials are computed correctly within numerical precision limits.
 *
 * \test GSLCompatTest::BinomialCoefficient
 * Tests binomial coefficient \f$\binom{n}{m} = \frac{n!}{m!(n-m)!}\f$. Verifies
 * correct computation including symmetry property \f$\binom{n}{m} = \binom{n}{n-m}\f$
 * and edge cases: \f$\binom{n}{0} = \binom{n}{n} = 1\f$.
 *
 * \test GSLCompatTest::BinomialCoefficientEdgeCases
 * Tests edge cases for binomial coefficients: when \f$m > n\f$ (should return 0),
 * and when \f$n = 0\f$ or \f$n = 1\f$ with various \f$m\f$ values.
 *
 * \test GSLCompatTest::Hypot
 * Tests the hypotenuse function \f$\text{hypot}(x, y) = \sqrt{x^2 + y^2}\f$.
 * Verifies correct computation for Pythagorean triples (e.g., 3-4-5, 5-12-13)
 * and edge cases with zero values.
 *
 * \test GSLCompatTest::ErrorHandler
 * Tests GSL error handler functions. Verifies that error handlers can be set
 * and disabled without causing crashes or exceptions.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include "Utilities/GSLCompat.h"

class GSLCompatTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Test setup
    }
};

TEST_F(GSLCompatTest, ErfFunction) {
    EXPECT_NEAR(gsl_sf_erf(0.0), 0.0, 1e-10);
    EXPECT_NEAR(gsl_sf_erf(1.0), std::erf(1.0), 1e-10);
    EXPECT_NEAR(gsl_sf_erf(-1.0), std::erf(-1.0), 1e-10);
    EXPECT_NEAR(gsl_sf_erf(2.0), std::erf(2.0), 1e-10);

    // erf(inf) should be 1
    EXPECT_NEAR(gsl_sf_erf(10.0), 1.0, 1e-6);
}

TEST_F(GSLCompatTest, PowInt) {
    EXPECT_DOUBLE_EQ(gsl_sf_pow_int(2.0, 0), 1.0);
    EXPECT_DOUBLE_EQ(gsl_sf_pow_int(2.0, 1), 2.0);
    EXPECT_DOUBLE_EQ(gsl_sf_pow_int(2.0, 2), 4.0);
    EXPECT_DOUBLE_EQ(gsl_sf_pow_int(2.0, 3), 8.0);
    EXPECT_DOUBLE_EQ(gsl_sf_pow_int(2.0, -1), 0.5);
    EXPECT_DOUBLE_EQ(gsl_sf_pow_int(2.0, -2), 0.25);

    EXPECT_DOUBLE_EQ(gsl_sf_pow_int(3.0, 4), 81.0);
    EXPECT_DOUBLE_EQ(gsl_sf_pow_int(5.0, 3), 125.0);
}

TEST_F(GSLCompatTest, GammaFunction) {
    EXPECT_NEAR(gsl_sf_gamma(1.0), 1.0, 1e-10);
    EXPECT_NEAR(gsl_sf_gamma(2.0), 1.0, 1e-10);
    EXPECT_NEAR(gsl_sf_gamma(3.0), 2.0, 1e-10);
    EXPECT_NEAR(gsl_sf_gamma(4.0), 6.0, 1e-10);
    EXPECT_NEAR(gsl_sf_gamma(0.5), std::sqrt(M_PI), 1e-10);
}

TEST_F(GSLCompatTest, Factorial) {
    EXPECT_DOUBLE_EQ(gsl_sf_fact(0), 1.0);
    EXPECT_DOUBLE_EQ(gsl_sf_fact(1), 1.0);
    EXPECT_DOUBLE_EQ(gsl_sf_fact(2), 2.0);
    EXPECT_DOUBLE_EQ(gsl_sf_fact(3), 6.0);
    EXPECT_DOUBLE_EQ(gsl_sf_fact(4), 24.0);
    EXPECT_DOUBLE_EQ(gsl_sf_fact(5), 120.0);
    EXPECT_DOUBLE_EQ(gsl_sf_fact(10), 3628800.0);
}

TEST_F(GSLCompatTest, FactorialLarge) {
    // Test large factorial (uses gamma function)
    double fact20   = gsl_sf_fact(20);
    double expected = 2432902008176640000.0;  // 20!
    EXPECT_NEAR(fact20, expected, 1e6);       // Allow some tolerance for large numbers
}

TEST_F(GSLCompatTest, BinomialCoefficient) {
    EXPECT_DOUBLE_EQ(gsl_sf_choose(5, 0), 1.0);
    EXPECT_DOUBLE_EQ(gsl_sf_choose(5, 1), 5.0);
    EXPECT_DOUBLE_EQ(gsl_sf_choose(5, 2), 10.0);
    EXPECT_DOUBLE_EQ(gsl_sf_choose(5, 3), 10.0);
    EXPECT_DOUBLE_EQ(gsl_sf_choose(5, 4), 5.0);
    EXPECT_DOUBLE_EQ(gsl_sf_choose(5, 5), 1.0);

    EXPECT_DOUBLE_EQ(gsl_sf_choose(10, 3), 120.0);
    EXPECT_DOUBLE_EQ(gsl_sf_choose(10, 7), 120.0);  // Symmetry
}

TEST_F(GSLCompatTest, BinomialCoefficientEdgeCases) {
    EXPECT_DOUBLE_EQ(gsl_sf_choose(5, 6), 0.0);  // m > n
    EXPECT_DOUBLE_EQ(gsl_sf_choose(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(gsl_sf_choose(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(gsl_sf_choose(1, 1), 1.0);
}

TEST_F(GSLCompatTest, Hypot) {
    EXPECT_NEAR(gsl_hypot(3.0, 4.0), 5.0, 1e-10);
    EXPECT_NEAR(gsl_hypot(5.0, 12.0), 13.0, 1e-10);
    EXPECT_NEAR(gsl_hypot(0.0, 0.0), 0.0, 1e-10);
    EXPECT_NEAR(gsl_hypot(1.0, 0.0), 1.0, 1e-10);
    EXPECT_NEAR(gsl_hypot(0.0, 1.0), 1.0, 1e-10);
}

TEST_F(GSLCompatTest, ErrorHandler) {
    gsl_error_handler_t old_handler = gsl_set_error_handler_off();
    // Should not crash
    gsl_set_error_handler(old_handler);
}
