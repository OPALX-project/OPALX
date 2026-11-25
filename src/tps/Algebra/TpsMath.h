#ifndef TPS_TpsMath_HH
#define TPS_TpsMath_HH
#include "Utilities/DomainError.h"
// Copyright: see LICENSE

// Overloaded template functions:
//   Elementary functions acting on Tps<T> objects.

#include <cmath>
#include "tps/Algebra/Tps.hpp"

namespace tps {

    // Class Tps; global functions acting on Tps objects.

    /// Integer power.
    template <class T>
    Tps<T> pow(const Tps<T>& x, int y);

    /// Square root.
    template <class T>
    Tps<T> sqrt(const Tps<T>& x);

    /// Sine.
    template <class T>
    Tps<T> sin(const Tps<T>& x);

    /// Cosine.
    template <class T>
    Tps<T> cos(const Tps<T>& x);

    /// Tangent.
    template <class T>
    Tps<T> tan(const Tps<T>& x);

    /// Secant.
    template <class T>
    Tps<T> sec(const Tps<T>& x);

    /// Cosecant.
    template <class T>
    Tps<T> csc(const Tps<T>& x);

    /// Exponential.
    template <class T>
    Tps<T> exp(const Tps<T>& x);

    /// Natural logarithm.
    template <class T>
    Tps<T> log(const Tps<T>& x);

    /// Hyperbolic sine.
    template <class T>
    Tps<T> sinh(const Tps<T>& x);

    /// Hyperbolic cosine.
    template <class T>
    Tps<T> cosh(const Tps<T>& x);

    /// Hyperbolic tangent.
    template <class T>
    Tps<T> tanh(const Tps<T>& x);

    // Implementation

    template <class T>
    Tps<T> pow(const Tps<T>& x, int y) {
        Tps<T> z(T(1));

        if (y > 0) {
            while (y-- > 0)
                z = z * x;
        } else if (y < 0) {
            Tps<T> t = x.inverse();
            while (y++ < 0)
                z = z * t;
        }

        return z;
    }

    template <class T>
    Tps<T> sqrt(const Tps<T>& x) {
        using std::sqrt;
        T aZero = x[0];
        if (aZero <= 0.0) throw ::DomainError("sqrt(const Tps &)");

        // Build Taylor series.
        int cut = x.getTaylorOrder();
        Array1D<T> series(cut + 1);
        series[0] = sqrt(aZero);

        for (int i = 1; i <= cut; i++) {
            series[i] = (series[i - 1] * double(3 - 2 * i)) / (aZero * double(2 * i));
        }

        Tps<T> z = x.taylor(series);
        return z;
    }

    template <class T>
    Tps<T> sin(const Tps<T>& x) {
        using std::cos;
        using std::sin;
        T aZero = x[0];
        int cut = x.getTaylorOrder();
        Array1D<T> series(cut + 1);
        series[0] = sin(aZero);
        if (cut >= 1) series[1] = cos(aZero);

        for (int i = 2; i <= cut; i++) {
            series[i] = -series[i - 2] / double(i * (i - 1));
        }

        Tps<T> z = x.taylor(series);
        return z;
    }

    template <class T>
    Tps<T> cos(const Tps<T>& x) {
        using std::cos;
        using std::sin;
        T aZero = x[0];
        int cut = x.getTaylorOrder();
        Array1D<T> series(cut + 1);
        series[0] = cos(aZero);
        if (cut >= 1) series[1] = -sin(aZero);

        for (int i = 2; i <= cut; i++) {
            series[i] = -series[i - 2] / double(i * (i - 1));
        }

        Tps<T> z = x.taylor(series);
        return z;
    }

    template <class T>
    Tps<T> tan(const Tps<T>& x) {
        return sin(x) / cos(x);
    }

    template <class T>
    Tps<T> sec(const Tps<T>& x) {
        return cos(x).inverse();
    }

    template <class T>
    Tps<T> csc(const Tps<T>& x) {
        return sin(x).inverse();
    }

    template <class T>
    Tps<T> exp(const Tps<T>& x) {
        using std::exp;
        T aZero = x[0];
        int cut = x.getTaylorOrder();
        Array1D<T> series(cut + 1);
        series[0] = exp(aZero);

        for (int i = 1; i <= cut; i++) {
            series[i] = series[i - 1] / double(i);
        }

        Tps<T> z = x.taylor(series);
        return z;
    }

    template <class T>
    Tps<T> log(const Tps<T>& x) {
        using std::log;
        T aZero = x[0];
        if (aZero <= 0) throw ::DomainError("log(const Tps &)");

        int cut = x.getTaylorOrder();
        double power;
        Array1D<T> series(cut + 1);
        series[0] = log(aZero);
        if (cut >= 1) series[1] = power = T(1) / aZero;

        for (int i = 2; i <= cut; i++) {
            power     = -power / aZero;
            series[i] = power / double(i);
        }

        Tps<T> z = x.taylor(series);
        return z;
    }

    template <class T>
    Tps<T> sinh(const Tps<T>& x) {
        using std::cosh;
        using std::sinh;
        T aZero = x[0];
        int cut = x.getTaylorOrder();
        Array1D<T> series(cut + 1);
        series[0] = sinh(aZero);
        if (cut >= 1) series[1] = cosh(aZero);

        for (int i = 2; i <= cut; i++) {
            series[i] = series[i - 2] / double(i * (i - 1));
        }

        Tps<T> z = x.taylor(series);
        return z;
    }

    template <class T>
    Tps<T> cosh(const Tps<T>& x) {
        using std::cosh;
        using std::sinh;
        T aZero = x[0];
        int cut = x.getTaylorOrder();
        Array1D<T> series(cut + 1);
        series[0] = cosh(aZero);
        if (cut >= 1) series[1] = sinh(aZero);

        for (int i = 2; i <= cut; i++) {
            series[i] = series[i - 2] / double(i * (i - 1));
        }

        Tps<T> z = x.taylor(series);
        return z;
    }

    template <class T>
    Tps<T> tanh(const Tps<T>& x) {
        return sinh(x) / cosh(x);
    }

}  // namespace tps

#endif  // TPS_TpsMath_HH
