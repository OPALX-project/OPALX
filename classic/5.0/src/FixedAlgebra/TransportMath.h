#ifndef CLASSIC_TransportMath_HH
#define CLASSIC_TransportMath_HH

// ------------------------------------------------------------------------
// $RCSfile: TransportMath.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Description:
//
// Declared template functions:
//   Elementary functions acting on TransportFun<T,N> objects.
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:50:59 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/TransportFun.h"
#include "Utilities/DomainError.h"


// Class TransportFun; global functions acting on TransportFun objects.
// These functions all work in the Transport approximation.
// ------------------------------------------------------------------------

/// Power (int y).
template <class T, int N>
TransportFun<T, N> pow(const TransportFun<T, N> &x, int y);

/// Square root.
template <class T, int N>
TransportFun<T, N> sqrt(const TransportFun<T, N> &x);

/// Sine.
template <class T, int N>
TransportFun<T, N> sin(const TransportFun<T, N> &x);

/// Cosine.
template <class T, int N>
TransportFun<T, N> cos(const TransportFun<T, N> &x);

/// Tangent.
template <class T, int N>
TransportFun<T, N> tan(const TransportFun<T, N> &x);

/// Cotangent.
template <class T, int N>
TransportFun<T, N> cot(const TransportFun<T, N> &x);

/// Secant.
template <class T, int N>
TransportFun<T, N> sec(const TransportFun<T, N> &x);

/// Cosecant.
template <class T, int N>
TransportFun<T, N> csc(const TransportFun<T, N> &x);

/// Exponential.
template <class T, int N>
TransportFun<T, N> exp(const TransportFun<T, N> &x);

/// Natural logarithm.
template <class T, int N>
TransportFun<T, N> log(const TransportFun<T, N> &x);

/// Hyperbolic sine.
template <class T, int N>
TransportFun<T, N> sinh(const TransportFun<T, N> &x);

/// Hyperbolic cosine.
template <class T, int N>
TransportFun<T, N> cosh(const TransportFun<T, N> &x);

/// Hyperbolic tangent.
template <class T, int N>
TransportFun<T, N> tanh(const TransportFun<T, N> &x);


// Implementation
// ------------------------------------------------------------------------

template <class T, int N>
TransportFun<T, N> pow(const TransportFun<T, N> &x, int y) {
    TransportFun<T, N> z;
    z[0] = std::pow(x[0], y);
    T fac = T(y) * std::pow(x[0], y - 1);
    for(int i = 1; i <= N; ++i) z[i] = fac * x[i];
    return z;
}


template <class T, int N>
TransportFun<T, N> sqrt(const TransportFun<T, N> &x) {
    T aZero = x[0];
    if(aZero <= 0.0) throw DomainError("sqrt(const TransportFun &)");
    T series[2];
    series[0] = sqrt(aZero);
    series[1] = series[0] / (2.0 * aZero);
    return x.taylor(series);
}


template <class T, int N>
TransportFun<T, N> sin(const TransportFun<T, N> &x) {
    T aZero = x[0];
    T series[2];
    series[0] = sin(aZero);
    series[1] = cos(aZero);
    return x.taylor(series);
}


template <class T, int N>
TransportFun<T, N> cos(const TransportFun<T, N> &x) {
    T aZero = x[0];
    T series[2];
    series[0] =   cos(aZero);
    series[1] = - sin(aZero);
    return x.taylor(series);
}


template <class T, int N>
TransportFun<T, N> tan(const TransportFun<T, N> &x) {
    return sin(x) / cos(x);
}


template <class T, int N>
TransportFun<T, N> sec(const TransportFun<T, N> &x) {
    return cos(x).inverse();
}


template <class T, int N>
TransportFun<T, N> csc(const TransportFun<T, N> &x) {
    return sin(x).inverse();
}


template <class T, int N>
TransportFun<T, N> exp(const TransportFun<T, N> &x) {
    T aZero = x[0];
    T series[2];
    series[0] = series[1] = exp(aZero);
    return x.taylor(series);
}


template <class T, int N>
TransportFun<T, N> log(const TransportFun<T, N> &x) {
    T aZero = x[0];
    if(aZero <= 0) throw DomainError("log(const TransportFun &)");

    T series[2];
    series[0] = log(aZero);
    series[1] = T(1) / aZero;
    return x.taylor(series);
}


template <class T, int N>
TransportFun<T, N> sinh(const TransportFun<T, N> &x) {
    T aZero = x[0];
    T series[2];
    series[0] = sinh(aZero);
    series[1] = cosh(aZero);
    return x.taylor(series);
}


template <class T, int N>
TransportFun<T, N> cosh(const TransportFun<T, N> &x) {
    T aZero = x[0];
    T series[2];
    series[0] = cosh(aZero);
    series[1] = sinh(aZero);
    return x.taylor(series);
}


template <class T, int N>
TransportFun<T, N> tanh(const TransportFun<T, N> &x) {
    return sinh(x) / cosh(x);
}

#endif // CLASSIC_TransportMath_HH
