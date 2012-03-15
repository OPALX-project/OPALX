#ifndef CLASSIC_LinearMath_HH
#define CLASSIC_LinearMath_HH

// ------------------------------------------------------------------------
// $RCSfile: LinearMath.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Description:
//
// Declared template functions:
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:18:06 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/LinearFun.h"
#include "Utilities/DomainError.h"


// Class LinearFun; global functions acting on LinearFun objects.
// These functions all work in the linear approximation.
// ------------------------------------------------------------------------

/// Power (int y).
template <class T, int N>
LinearFun<T, N> pow(const LinearFun<T, N> &x, int y);

/// Square root.
template <class T, int N>
LinearFun<T, N> sqrt(const LinearFun<T, N> &x);

/// Sine.
template <class T, int N>
LinearFun<T, N> sin(const LinearFun<T, N> &x);

/// Cosine.
template <class T, int N>
LinearFun<T, N> cos(const LinearFun<T, N> &x);

/// Tangent.
template <class T, int N>
LinearFun<T, N> tan(const LinearFun<T, N> &x);

/// Cotangent.
template <class T, int N>
LinearFun<T, N> cot(const LinearFun<T, N> &x);

/// Secant.
template <class T, int N>
LinearFun<T, N> sec(const LinearFun<T, N> &x);

/// Cosecant.
template <class T, int N>
LinearFun<T, N> csc(const LinearFun<T, N> &x);

/// Exponential.
template <class T, int N>
LinearFun<T, N> exp(const LinearFun<T, N> &x);

/// Natural logarithm.
template <class T, int N>
LinearFun<T, N> log(const LinearFun<T, N> &x);

/// Hyperbolic sine.
template <class T, int N>
LinearFun<T, N> sinh(const LinearFun<T, N> &x);

/// Hyperbolic cosine.
template <class T, int N>
LinearFun<T, N> cosh(const LinearFun<T, N> &x);

/// Hyperbolic tangent.
template <class T, int N>
LinearFun<T, N> tanh(const LinearFun<T, N> &x);


// Implementation
// ------------------------------------------------------------------------

template <class T, int N>
LinearFun<T, N> pow(const LinearFun<T, N> &x, int y) {
    LinearFun<T, N> z;
    z[0] = std::pow(x[0], y);
    T fac = T(y) * std::pow(x[0], y - 1);
    for(int i = 1; i <= N; ++i) z[i] = fac * x[i];
    return z;
}


template <class T, int N>
LinearFun<T, N> sqrt(const LinearFun<T, N> &x) {
    T aZero = x[0];
    if(aZero <= 0.0) throw DomainError("sqrt(const LinearFun &)");
    T series[2];
    series[0] = sqrt(aZero);
    series[1] = series[0] / (2.0 * aZero);
    return x.taylor(series);
}


template <class T, int N>
LinearFun<T, N> sin(const LinearFun<T, N> &x) {
    T aZero = x[0];
    T series[2];
    series[0] = sin(aZero);
    series[1] = cos(aZero);
    return x.taylor(series);
}


template <class T, int N>
LinearFun<T, N> cos(const LinearFun<T, N> &x) {
    T aZero = x[0];
    T series[2];
    series[0] =   cos(aZero);
    series[1] = - sin(aZero);
    return x.taylor(series);
}


template <class T, int N>
LinearFun<T, N> tan(const LinearFun<T, N> &x) {
    return sin(x) / cos(x);
}


template <class T, int N>
LinearFun<T, N> sec(const LinearFun<T, N> &x) {
    return cos(x).inverse();
}


template <class T, int N>
LinearFun<T, N> csc(const LinearFun<T, N> &x) {
    return sin(x).inverse();
}


template <class T, int N>
LinearFun<T, N> exp(const LinearFun<T, N> &x) {
    T aZero = x[0];
    T series[2];
    series[0] = series[1] = exp(aZero);
    return x.taylor(series);
}


template <class T, int N>
LinearFun<T, N> log(const LinearFun<T, N> &x) {
    T aZero = x[0];
    if(aZero <= 0) throw DomainError("log(const LinearFun &)");

    T series[2];
    series[0] = log(aZero);
    series[1] = T(1) / aZero;
    return x.taylor(series);
}


template <class T, int N>
LinearFun<T, N> sinh(const LinearFun<T, N> &x) {
    T aZero = x[0];
    T series[2];
    series[0] = sinh(aZero);
    series[1] = cosh(aZero);
    return x.taylor(series);
}


template <class T, int N>
LinearFun<T, N> cosh(const LinearFun<T, N> &x) {
    T aZero = x[0];
    T series[2];
    series[0] = cosh(aZero);
    series[1] = sinh(aZero);
    return x.taylor(series);
}


template <class T, int N>
LinearFun<T, N> tanh(const LinearFun<T, N> &x) {
    return sinh(x) / cosh(x);
}

#endif // CLASSIC_LinearMath_HH
