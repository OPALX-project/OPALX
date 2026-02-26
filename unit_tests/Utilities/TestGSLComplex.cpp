/**
 * \file TestGSLComplex.cpp
 * \brief Unit tests for GSL complex number operations
 *
 * This file contains unit tests for GSL-compatible complex number operations.
 * The implementation provides a wrapper around \c std::complex<double> with
 * GSL-compatible function names and interfaces.
 *
 * \test GSLComplexTest::CreateRectangular
 * Tests creation of complex numbers in rectangular form \f$z = a + bi\f$.
 * Verifies that real and imaginary parts are correctly stored and can be
 * retrieved using \c GSL_REAL and \c GSL_IMAG macros.
 *
 * \test GSLComplexTest::CreatePolar
 * Tests creation of complex numbers in polar form \f$z = r e^{i\theta}\f$.
 * Verifies that polar coordinates are correctly converted to rectangular form:
 * \f$a = r\cos\theta\f$, \f$b = r\sin\theta\f$.
 *
 * \test GSLComplexTest::Addition
 * Tests complex number addition: \f$(a_1 + b_1i) + (a_2 + b_2i) = (a_1+a_2) + (b_1+b_2)i\f$.
 * Verifies correct computation of real and imaginary parts.
 *
 * \test GSLComplexTest::Subtraction
 * Tests complex number subtraction: \f$(a_1 + b_1i) - (a_2 + b_2i) = (a_1-a_2) + (b_1-b_2)i\f$.
 *
 * \test GSLComplexTest::Multiplication
 * Tests complex number multiplication using the identity:
 * \f$(a_1 + b_1i)(a_2 + b_2i) = (a_1a_2 - b_1b_2) + (a_1b_2 + a_2b_1)i\f$.
 *
 * \test GSLComplexTest::Division
 * Tests complex number division. Verifies special case: \f$1/i = -i\f$.
 *
 * \test GSLComplexTest::Conjugate
 * Tests complex conjugate: \f$\overline{a + bi} = a - bi\f$. Verifies that
 * the real part remains unchanged and the imaginary part is negated.
 *
 * \test GSLComplexTest::Absolute
 * Tests complex absolute value (modulus): \f$|a + bi| = \sqrt{a^2 + b^2}\f$.
 * Verifies correct computation for a 3-4-5 right triangle case.
 *
 * \test GSLComplexTest::Argument
 * Tests complex argument (phase): \f$\arg(a + bi) = \arctan(b/a)\f$.
 * Verifies correct computation for \f$1 + i\f$ (should be \f$\pi/4\f$).
 *
 * \test GSLComplexTest::Tanh
 * Tests complex hyperbolic tangent function. Verifies behavior for real
 * arguments (\f$\tanh(1)\f$) and purely imaginary arguments
 * (\f$\tanh(i) = i\tan(1)\f$).
 *
 * \test GSLComplexTest::Exp
 * Tests complex exponential function. Verifies Euler's identity:
 * \f$e^{i\pi} = -1\f$.
 *
 * \test GSLComplexTest::Log
 * Tests complex logarithm function. Verifies that \f$\ln(1) = 0\f$ for the
 * principal branch.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include "Utilities/GSLComplex.h"

class GSLComplexTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Test setup
    }
};

TEST_F(GSLComplexTest, CreateRectangular) {
    gsl_complex z = gsl_complex_rect(3.0, 4.0);
    EXPECT_DOUBLE_EQ(GSL_REAL(z), 3.0);
    EXPECT_DOUBLE_EQ(GSL_IMAG(z), 4.0);
}

TEST_F(GSLComplexTest, CreatePolar) {
    gsl_complex z = gsl_complex_polar(5.0, M_PI / 4.0);
    EXPECT_NEAR(GSL_REAL(z), 5.0 * std::cos(M_PI / 4.0), 1e-10);
    EXPECT_NEAR(GSL_IMAG(z), 5.0 * std::sin(M_PI / 4.0), 1e-10);
}

TEST_F(GSLComplexTest, Addition) {
    gsl_complex z1  = gsl_complex_rect(1.0, 2.0);
    gsl_complex z2  = gsl_complex_rect(3.0, 4.0);
    gsl_complex sum = gsl_complex_add(z1, z2);

    EXPECT_DOUBLE_EQ(GSL_REAL(sum), 4.0);
    EXPECT_DOUBLE_EQ(GSL_IMAG(sum), 6.0);
}

TEST_F(GSLComplexTest, Subtraction) {
    gsl_complex z1   = gsl_complex_rect(5.0, 6.0);
    gsl_complex z2   = gsl_complex_rect(2.0, 3.0);
    gsl_complex diff = gsl_complex_sub(z1, z2);

    EXPECT_DOUBLE_EQ(GSL_REAL(diff), 3.0);
    EXPECT_DOUBLE_EQ(GSL_IMAG(diff), 3.0);
}

TEST_F(GSLComplexTest, Multiplication) {
    gsl_complex z1   = gsl_complex_rect(1.0, 2.0);
    gsl_complex z2   = gsl_complex_rect(3.0, 4.0);
    gsl_complex prod = gsl_complex_mul(z1, z2);

    // (1+2i)(3+4i) = 3 + 4i + 6i + 8i^2 = 3 + 10i - 8 = -5 + 10i
    EXPECT_DOUBLE_EQ(GSL_REAL(prod), -5.0);
    EXPECT_DOUBLE_EQ(GSL_IMAG(prod), 10.0);
}

TEST_F(GSLComplexTest, Division) {
    gsl_complex z1   = gsl_complex_rect(1.0, 0.0);
    gsl_complex z2   = gsl_complex_rect(0.0, 1.0);
    gsl_complex quot = gsl_complex_div(z1, z2);

    // 1 / i = -i
    EXPECT_NEAR(GSL_REAL(quot), 0.0, 1e-10);
    EXPECT_NEAR(GSL_IMAG(quot), -1.0, 1e-10);
}

TEST_F(GSLComplexTest, Conjugate) {
    gsl_complex z    = gsl_complex_rect(3.0, 4.0);
    gsl_complex conj = gsl_complex_conjugate(z);

    EXPECT_DOUBLE_EQ(GSL_REAL(conj), 3.0);
    EXPECT_DOUBLE_EQ(GSL_IMAG(conj), -4.0);
}

TEST_F(GSLComplexTest, Absolute) {
    gsl_complex z = gsl_complex_rect(3.0, 4.0);
    double abs    = gsl_complex_abs(z);

    EXPECT_NEAR(abs, 5.0, 1e-10);
}

TEST_F(GSLComplexTest, Argument) {
    gsl_complex z = gsl_complex_rect(1.0, 1.0);
    double arg    = gsl_complex_arg(z);

    EXPECT_NEAR(arg, M_PI / 4.0, 1e-10);
}

TEST_F(GSLComplexTest, Tanh) {
    gsl_complex z      = gsl_complex_rect(1.0, 0.0);
    gsl_complex tanh_z = gsl_complex_tanh(z);

    // tanh(1) ≈ 0.7616
    EXPECT_NEAR(GSL_REAL(tanh_z), std::tanh(1.0), 1e-10);
    EXPECT_NEAR(GSL_IMAG(tanh_z), 0.0, 1e-10);

    // Test with imaginary part
    z      = gsl_complex_rect(0.0, 1.0);
    tanh_z = gsl_complex_tanh(z);
    // tanh(i) = i*tan(1)
    EXPECT_NEAR(GSL_REAL(tanh_z), 0.0, 1e-10);
    EXPECT_NEAR(GSL_IMAG(tanh_z), std::tan(1.0), 1e-10);
}

TEST_F(GSLComplexTest, Exp) {
    gsl_complex z     = gsl_complex_rect(0.0, M_PI);
    gsl_complex exp_z = gsl_complex_exp(z);

    // e^(i*pi) = -1
    EXPECT_NEAR(GSL_REAL(exp_z), -1.0, 1e-10);
    EXPECT_NEAR(GSL_IMAG(exp_z), 0.0, 1e-10);
}

TEST_F(GSLComplexTest, Log) {
    gsl_complex z     = gsl_complex_rect(1.0, 0.0);
    gsl_complex log_z = gsl_complex_log(z);

    // ln(1) = 0
    EXPECT_NEAR(GSL_REAL(log_z), 0.0, 1e-10);
    EXPECT_NEAR(GSL_IMAG(log_z), 0.0, 1e-10);
}
