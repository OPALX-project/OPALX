//
// GSL compatibility stubs
//
// Copyright (c) 2023, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#ifndef OPAL_GSL_COMPAT_HH
#define OPAL_GSL_COMPAT_HH

#include <cmath>
#include <limits>

/// \see https://www.gnu.org/software/gsl/doc/html/err.html
/// \brief Signature of a GSL-style error handler callback.
typedef void (*gsl_error_handler_t)(const char* reason, const char* file, int line, int gsl_errno);

/// \brief Set a global error handler and return the previous handler.
/// \param handler Input: new error handler.
/// \return Output: previous handler (may be null).
inline gsl_error_handler_t gsl_set_error_handler(gsl_error_handler_t handler) {
    // Store and return previous handler (simplified - no actual storage)
    static gsl_error_handler_t prev_handler = nullptr;
    gsl_error_handler_t old = prev_handler;
    prev_handler = handler;
    return old;
}

/// \brief Disable error handling and return a null handler.
/// \return Output: null handler.
inline gsl_error_handler_t gsl_set_error_handler_off() {
    // Return a null handler
    return nullptr;
}

/// \brief Error function \f$\mathrm{erf}(x)\f$.
/// \param x Input: argument.
/// \return Output: \f$\mathrm{erf}(x)\f$.
inline double gsl_sf_erf(double x) {
    return std::erf(x);
}

/// \see https://www.gnu.org/software/gsl/doc/html/specfunc.html
/// \brief Integer power \f$x^n\f$.
/// \param x Input: base.
/// \param n Input: integer exponent.
/// \return Output: \f$x^n\f$.
inline double gsl_sf_pow_int(double x, int n) {
    return std::pow(x, n);
}

/// \see https://www.gnu.org/software/gsl/doc/html/specfunc.html
/// \brief Gamma function \f$\Gamma(x)\f$.
/// \param x Input: argument.
/// \return Output: \f$\Gamma(x)\f$.
inline double gsl_sf_gamma(double x) {
    return std::tgamma(x);
}

/// \see https://www.gnu.org/software/gsl/doc/html/specfunc.html
/// \brief Factorial \f$n!\f$ (computed via \f$\Gamma(n+1)\f$ for large \f$n\f$).
/// \param n Input: non-negative integer.
/// \return Output: \f$n!\f$.
inline double gsl_sf_fact(unsigned int n) {
    if (n > 170) {
        // Factorial of numbers > 170 overflows double
        // Use gamma function for large numbers
        return std::tgamma(static_cast<double>(n + 1));
    }
    double result = 1.0;
    for (unsigned int i = 2; i <= n; ++i) {
        result *= static_cast<double>(i);
    }
    return result;
}

/// \see https://www.gnu.org/software/gsl/doc/html/specfunc.html
/// \brief Binomial coefficient \f$\binom{n}{m}\f$.
/// \param n Input: total count.
/// \param m Input: selected count.
/// \return Output: \f$\binom{n}{m}\f$.
inline double gsl_sf_choose(unsigned int n, unsigned int m) {
    if (m > n) return 0.0;
    if (m == 0 || m == n) return 1.0;
    
    // Use symmetry: choose(n, m) = choose(n, n-m)
    if (m > n / 2) m = n - m;
    
    // Compute using iterative method to avoid overflow
    double result = 1.0;
    for (unsigned int i = 0; i < m; ++i) {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

/// \brief GSL-compatible math constants.
namespace gsl {
    constexpr double GSL_POSINF = std::numeric_limits<double>::infinity();
    constexpr double GSL_NEGINF = -std::numeric_limits<double>::infinity();
    constexpr double GSL_NAN = std::numeric_limits<double>::quiet_NaN();
}

/// \see https://www.gnu.org/software/gsl/doc/html/specfunc.html
/// \brief Hypotenuse \f$\sqrt{x^2 + y^2}\f$ with overflow-safe std::hypot.
/// \param x Input: first leg.
/// \param y Input: second leg.
/// \return Output: \f$\sqrt{x^2 + y^2}\f$.
inline double gsl_hypot(double x, double y) {
    return std::hypot(x, y);
}

// GSL_REAL and GSL_IMAG macros for complex numbers
#define GSL_REAL(z) ((z).dat[0])
#define GSL_IMAG(z) ((z).dat[1])

#endif // OPAL_GSL_COMPAT_HH

