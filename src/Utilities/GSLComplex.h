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

#include <cmath>
#include <complex>

/// \see Documentation on https://www.gnu.org/software/gsl/doc/html/complex.htm
/// \see Implementation on https://www.gnu.org/software/gsl/
/// \brief Complex number stored as \f$(\Re, \Im)\f$.
struct gsl_complex {
    double dat[2];  // [real, imag]

    gsl_complex() : dat{0.0, 0.0} {}
    gsl_complex(double r, double i) : dat{r, i} {}
    gsl_complex(const std::complex<double>& c) : dat{c.real(), c.imag()} {}

    operator std::complex<double>() const { return std::complex<double>(dat[0], dat[1]); }

    double real() const { return dat[0]; }
    double imag() const { return dat[1]; }
};

/// \brief Access real part of a \c gsl_complex.
#define GSL_REAL(z) ((z).dat[0])
/// \brief Access imaginary part of a \c gsl_complex.
#define GSL_IMAG(z) ((z).dat[1])

/// \brief Construct \f$z = x + i y\f$.
/// \param x Input: real part.
/// \param y Input: imaginary part.
/// \return Output: complex value \f$x + i y\f$.
inline gsl_complex gsl_complex_rect(double x, double y) { return gsl_complex(x, y); }

/// \brief Construct \f$z = r e^{i\theta}\f$.
/// \param r Input: magnitude.
/// \param theta Input: phase angle.
/// \return Output: complex value \f$r e^{i\theta}\f$.
inline gsl_complex gsl_complex_polar(double r, double theta) {
    return gsl_complex(r * std::cos(theta), r * std::sin(theta));
}

/// \brief Argument (phase) \f$\arg(z)\f$.
/// \param z Input: complex value.
/// \return Output: phase angle.
inline double gsl_complex_arg(gsl_complex z) { return std::arg(std::complex<double>(z)); }

/// \brief Magnitude \f$|z|\f$.
/// \param z Input: complex value.
/// \return Output: \f$|z|\f$.
inline double gsl_complex_abs(gsl_complex z) { return std::abs(std::complex<double>(z)); }

/// \brief Squared magnitude \f$|z|^2\f$.
/// \param z Input: complex value.
/// \return Output: \f$|z|^2\f$.
inline double gsl_complex_abs2(gsl_complex z) { return z.dat[0] * z.dat[0] + z.dat[1] * z.dat[1]; }

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Sum \f$a + b\f$.
/// \param a Input: first operand.
/// \param b Input: second operand.
/// \return Output: \f$a + b\f$.
inline gsl_complex gsl_complex_add(gsl_complex a, gsl_complex b) {
    return gsl_complex(a.dat[0] + b.dat[0], a.dat[1] + b.dat[1]);
}

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Difference \f$a - b\f$.
/// \param a Input: first operand.
/// \param b Input: second operand.
/// \return Output: \f$a - b\f$.
inline gsl_complex gsl_complex_sub(gsl_complex a, gsl_complex b) {
    return gsl_complex(a.dat[0] - b.dat[0], a.dat[1] - b.dat[1]);
}

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Product \f$a \cdot b\f$.
/// \param a Input: first operand.
/// \param b Input: second operand.
/// \return Output: \f$a b\f$.
inline gsl_complex gsl_complex_mul(gsl_complex a, gsl_complex b) {
    std::complex<double> ca(a), cb(b);
    return gsl_complex(ca * cb);
}

/// \brief Quotient \f$a / b\f$.
/// \param a Input: numerator.
/// \param b Input: denominator.
/// \return Output: \f$a / b\f$.
inline gsl_complex gsl_complex_div(gsl_complex a, gsl_complex b) {
    std::complex<double> ca(a), cb(b);
    return gsl_complex(ca / cb);
}

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Add real scalar \f$a + x\f$.
/// \param a Input: complex operand.
/// \param x Input: real scalar.
/// \return Output: \f$a + x\f$.
inline gsl_complex gsl_complex_add_real(gsl_complex a, double x) {
    return gsl_complex(a.dat[0] + x, a.dat[1]);
}

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Subtract real scalar \f$a - x\f$.
/// \param a Input: complex operand.
/// \param x Input: real scalar.
/// \return Output: \f$a - x\f$.
inline gsl_complex gsl_complex_sub_real(gsl_complex a, double x) {
    return gsl_complex(a.dat[0] - x, a.dat[1]);
}

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Multiply by real scalar \f$a \cdot x\f$.
/// \param a Input: complex operand.
/// \param x Input: real scalar.
/// \return Output: \f$a x\f$.
inline gsl_complex gsl_complex_mul_real(gsl_complex a, double x) {
    return gsl_complex(a.dat[0] * x, a.dat[1] * x);
}

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Divide by real scalar \f$a / x\f$.
/// \param a Input: complex operand.
/// \param x Input: real scalar.
/// \return Output: \f$a / x\f$.
inline gsl_complex gsl_complex_div_real(gsl_complex a, double x) {
    return gsl_complex(a.dat[0] / x, a.dat[1] / x);
}

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Complex conjugate \f$\overline{a}\f$.
/// \param a Input: complex operand.
/// \return Output: \f$\overline{a}\f$.
inline gsl_complex gsl_complex_conjugate(gsl_complex a) { return gsl_complex(a.dat[0], -a.dat[1]); }

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Multiplicative inverse \f$1/a\f$.
/// \param a Input: complex operand.
/// \return Output: \f$1/a\f$.
inline gsl_complex gsl_complex_inverse(gsl_complex a) {
    std::complex<double> ca(a);
    return gsl_complex(1.0 / ca);
}

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Negation \f$-a\f$.
/// \param a Input: complex operand.
/// \return Output: \f$-a\f$.
inline gsl_complex gsl_complex_negative(gsl_complex a) { return gsl_complex(-a.dat[0], -a.dat[1]); }

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Complex square root \f$\sqrt{a}\f$.
/// \param a Input: complex operand.
/// \return Output: \f$\sqrt{a}\f$.
inline gsl_complex gsl_complex_sqrt(gsl_complex a) {
    std::complex<double> ca(a);
    return gsl_complex(std::sqrt(ca));
}

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Complex exponential \f$e^{a}\f$.
/// \param a Input: complex operand.
/// \return Output: \f$e^{a}\f$.
inline gsl_complex gsl_complex_exp(gsl_complex a) {
    std::complex<double> ca(a);
    return gsl_complex(std::exp(ca));
}

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Complex natural logarithm \f$\log(a)\f$.
/// \param a Input: complex operand.
/// \return Output: \f$\log(a)\f$.
inline gsl_complex gsl_complex_log(gsl_complex a) {
    std::complex<double> ca(a);
    return gsl_complex(std::log(ca));
}

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Complex sine \f$\sin(a)\f$.
/// \param a Input: complex operand.
/// \return Output: \f$\sin(a)\f$.
inline gsl_complex gsl_complex_sin(gsl_complex a) {
    std::complex<double> ca(a);
    return gsl_complex(std::sin(ca));
}

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Complex cosine \f$\cos(a)\f$.
/// \param a Input: complex operand.
/// \return Output: \f$\cos(a)\f$.
inline gsl_complex gsl_complex_cos(gsl_complex a) {
    std::complex<double> ca(a);
    return gsl_complex(std::cos(ca));
}

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Complex power \f$a^b\f$.
/// \param a Input: base.
/// \param b Input: exponent.
/// \return Output: \f$a^b\f$.
inline gsl_complex gsl_complex_pow(gsl_complex a, gsl_complex b) {
    std::complex<double> ca(a), cb(b);
    return gsl_complex(std::pow(ca, cb));
}

/// \see https://www.gnu.org/software/gsl/doc/html/complex.html
/// \brief Complex hyperbolic tangent \f$\tanh(a)\f$.
/// \param a Input: complex operand.
/// \return Output: \f$\tanh(a)\f$.
inline gsl_complex gsl_complex_tanh(gsl_complex a) {
    std::complex<double> ca(a);
    return gsl_complex(std::tanh(ca));
}

#endif  // OPAL_GSL_COMPLEX_HH
