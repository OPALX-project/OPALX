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

// Error handling stubs
typedef void (*gsl_error_handler_t)(const char* reason, const char* file, int line, int gsl_errno);

inline gsl_error_handler_t gsl_set_error_handler(gsl_error_handler_t handler) {
    // Store and return previous handler (simplified - no actual storage)
    static gsl_error_handler_t prev_handler = nullptr;
    gsl_error_handler_t old = prev_handler;
    prev_handler = handler;
    return old;
}

inline gsl_error_handler_t gsl_set_error_handler_off() {
    // Return a null handler
    return nullptr;
}

// Special functions - use standard library
inline double gsl_sf_erf(double x) {
    return std::erf(x);
}

// Power function
inline double gsl_sf_pow_int(double x, int n) {
    return std::pow(x, n);
}

// Gamma function
inline double gsl_sf_gamma(double x) {
    return std::tgamma(x);
}

// Factorial function
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

// Binomial coefficient: choose(n, m) = n! / (m! * (n-m)!)
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

// Math constants and utilities
namespace gsl {
    constexpr double GSL_POSINF = std::numeric_limits<double>::infinity();
    constexpr double GSL_NEGINF = -std::numeric_limits<double>::infinity();
    constexpr double GSL_NAN = std::numeric_limits<double>::quiet_NaN();
}

// Hypotenuse function: sqrt(x^2 + y^2)
inline double gsl_hypot(double x, double y) {
    return std::hypot(x, y);
}

// GSL_REAL and GSL_IMAG macros for complex numbers
#define GSL_REAL(z) ((z).dat[0])
#define GSL_IMAG(z) ((z).dat[1])

#endif // OPAL_GSL_COMPAT_HH

