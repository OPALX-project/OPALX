//
// GSL Complex compatibility to replace gsl_complex
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

#ifndef OPAL_GSL_COMPLEX_HH
#define OPAL_GSL_COMPLEX_HH

#include <complex>
#include <cmath>

// GSL complex type compatibility
struct gsl_complex {
    double dat[2];  // [real, imag]
    
    gsl_complex() : dat{0.0, 0.0} {}
    gsl_complex(double r, double i) : dat{r, i} {}
    gsl_complex(const std::complex<double>& c) : dat{c.real(), c.imag()} {}
    
    operator std::complex<double>() const {
        return std::complex<double>(dat[0], dat[1]);
    }
    
    double real() const { return dat[0]; }
    double imag() const { return dat[1]; }
};

// GSL macros for accessing real and imaginary parts
#define GSL_REAL(z) ((z).dat[0])
#define GSL_IMAG(z) ((z).dat[1])

// GSL complex math functions
inline gsl_complex gsl_complex_rect(double x, double y) {
    return gsl_complex(x, y);
}

inline gsl_complex gsl_complex_polar(double r, double theta) {
    return gsl_complex(r * std::cos(theta), r * std::sin(theta));
}

inline double gsl_complex_arg(gsl_complex z) {
    return std::arg(std::complex<double>(z));
}

inline double gsl_complex_abs(gsl_complex z) {
    return std::abs(std::complex<double>(z));
}

inline double gsl_complex_abs2(gsl_complex z) {
    return z.dat[0] * z.dat[0] + z.dat[1] * z.dat[1];
}

inline gsl_complex gsl_complex_add(gsl_complex a, gsl_complex b) {
    return gsl_complex(a.dat[0] + b.dat[0], a.dat[1] + b.dat[1]);
}

inline gsl_complex gsl_complex_sub(gsl_complex a, gsl_complex b) {
    return gsl_complex(a.dat[0] - b.dat[0], a.dat[1] - b.dat[1]);
}

inline gsl_complex gsl_complex_mul(gsl_complex a, gsl_complex b) {
    std::complex<double> ca(a), cb(b);
    return gsl_complex(ca * cb);
}

inline gsl_complex gsl_complex_div(gsl_complex a, gsl_complex b) {
    std::complex<double> ca(a), cb(b);
    return gsl_complex(ca / cb);
}

inline gsl_complex gsl_complex_add_real(gsl_complex a, double x) {
    return gsl_complex(a.dat[0] + x, a.dat[1]);
}

inline gsl_complex gsl_complex_sub_real(gsl_complex a, double x) {
    return gsl_complex(a.dat[0] - x, a.dat[1]);
}

inline gsl_complex gsl_complex_mul_real(gsl_complex a, double x) {
    return gsl_complex(a.dat[0] * x, a.dat[1] * x);
}

inline gsl_complex gsl_complex_div_real(gsl_complex a, double x) {
    return gsl_complex(a.dat[0] / x, a.dat[1] / x);
}

inline gsl_complex gsl_complex_conjugate(gsl_complex a) {
    return gsl_complex(a.dat[0], -a.dat[1]);
}

inline gsl_complex gsl_complex_inverse(gsl_complex a) {
    std::complex<double> ca(a);
    return gsl_complex(1.0 / ca);
}

inline gsl_complex gsl_complex_negative(gsl_complex a) {
    return gsl_complex(-a.dat[0], -a.dat[1]);
}

inline gsl_complex gsl_complex_sqrt(gsl_complex a) {
    std::complex<double> ca(a);
    return gsl_complex(std::sqrt(ca));
}

inline gsl_complex gsl_complex_exp(gsl_complex a) {
    std::complex<double> ca(a);
    return gsl_complex(std::exp(ca));
}

inline gsl_complex gsl_complex_log(gsl_complex a) {
    std::complex<double> ca(a);
    return gsl_complex(std::log(ca));
}

inline gsl_complex gsl_complex_sin(gsl_complex a) {
    std::complex<double> ca(a);
    return gsl_complex(std::sin(ca));
}

inline gsl_complex gsl_complex_cos(gsl_complex a) {
    std::complex<double> ca(a);
    return gsl_complex(std::cos(ca));
}

inline gsl_complex gsl_complex_pow(gsl_complex a, gsl_complex b) {
    std::complex<double> ca(a), cb(b);
    return gsl_complex(std::pow(ca, cb));
}

inline gsl_complex gsl_complex_tanh(gsl_complex a) {
    std::complex<double> ca(a);
    return gsl_complex(std::tanh(ca));
}

#endif // OPAL_GSL_COMPLEX_HH

