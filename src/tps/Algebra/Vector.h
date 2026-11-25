#ifndef TPS_VECTOR_HPP
#define TPS_VECTOR_HPP
#include "Utilities/LogicalError.h"
// Copyright: see LICENSE (OPAL/CLASSIC API; Kokkos Kernels implementation).
#include <KokkosBatched_Axpy.hpp>
#include <KokkosBatched_Dot.hpp>
#include <KokkosBlas1_mult.hpp>
#include <KokkosBlas1_nrm2.hpp>
#include <KokkosBlas1_scal.hpp>
#include "tps/Algebra/Array1D.h"
namespace tps {
    /** @brief Checked host vector; numerical operations delegate to Kokkos Kernels. */
    template <class T>
    class Vector : public Array1D<T> {
    public:
        using Array1D<T>::Array1D;
        Vector() = default;
        Vector(const Array1D<T>& x) : Array1D<T>(x) {}
        Vector operator-() const {
            Vector r(*this);
            r *= T(-1);
            return r;
        }
        Vector& operator*=(const T& a) {
            KokkosBlas::SerialScale::invoke(a, this->view());
            return *this;
        }
        Vector& operator/=(const T& a) { return *this *= T(1) / a; }
        Vector& operator+=(const Vector& b) {
            check(b);
            KokkosBatched::SerialAxpy::invoke(T(1), b.view(), this->view());
            return *this;
        }
        Vector& operator-=(const Vector& b) {
            check(b);
            KokkosBatched::SerialAxpy::invoke(T(-1), b.view(), this->view());
            return *this;
        }

    private:
        void check(const Vector& b) const {
            if (this->size() != b.size()) throw ::LogicalError("Vector", "Dimension mismatch");
        }
    };
    template <class T>
    Vector<T> operator+(Vector<T> a, const Vector<T>& b) {
        return a += b;
    }
    template <class T>
    Vector<T> operator-(Vector<T> a, const Vector<T>& b) {
        return a -= b;
    }
    template <class T>
    Vector<T> operator*(Vector<T> a, const T& b) {
        return a *= b;
    }
    template <class T>
    Vector<T> operator*(const T& a, Vector<T> b) {
        return b *= a;
    }
    template <class T>
    Vector<T> operator/(Vector<T> a, const T& b) {
        return a /= b;
    }
    /** @brief Bilinear dot product; complex coefficients are NOT conjugated. */
    template <class T>
    T operator*(const Vector<T>& a, const Vector<T>& b) {
        if (a.size() != b.size()) throw ::LogicalError("Vector dot", "Dimension mismatch");
        T result{};
        Kokkos::View<T, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scalar(&result);
        KokkosBatched::SerialDot<KokkosBatched::Trans::Transpose, 0>::invoke(
                a.view(), b.view(), scalar);
        return result;
    }
    /** @brief Euclidean magnitude, with conjugation for complex input. */
    template <class T>
    auto euclidean_norm(const Vector<T>& a) {
        return KokkosBlas::nrm2(a.view());
    }
    template <class T>
    auto scaled_norm(const Array1D<T>& d, const Vector<T>& a) {
        if (d.size() != a.size()) throw ::LogicalError("scaled_norm", "Dimension mismatch");
        Vector<T> temporary(a.size());
        KokkosBlas::mult(T(0), temporary.view(), T(1), d.view(), a.view());
        Kokkos::fence();
        return euclidean_norm(temporary);
    }
}  // namespace tps
#endif
