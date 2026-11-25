#ifndef TPS_MATRIX_HPP
#define TPS_MATRIX_HPP
#include "Utilities/LogicalError.h"
// Copyright: see LICENSE (OPAL/CLASSIC API; Kokkos Kernels implementation).
#include <KokkosBatched_Copy_Decl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas2_serial_gemv.hpp>
#include "tps/Algebra/Array2D.h"
#include "tps/Algebra/Vector.h"
namespace tps {
    /** @brief Checked row-major host matrix backed by Kokkos storage.
     * @details Matrix(r,c,value) fills ALL entries with value. Use identity(n) for
     * an identity matrix. This replaces the old constructor that ignored value.
     */
    template <class T>
    class Matrix : public Array2D<T> {
    public:
        using Array2D<T>::Array2D;
        Matrix() = default;
        Matrix(const Array2D<T>& a) : Array2D<T>(a) {}
        using Array2D<T>::operator=;
        static Matrix identity(int n) {
            Matrix r(n, n, T(0));
            for (int i = 0; i < n; ++i)
                r(i, i) = T(1);
            return r;
        }
        Matrix& operator*=(const T& a) {
            KokkosBlas::SerialScale::invoke(a, this->view());
            return *this;
        }
        Matrix& operator/=(const T& a) { return *this *= T(1) / a; }
        Matrix& operator+=(const Matrix& b) {
            check(b);
            KokkosBatched::SerialAxpy::invoke(T(1), b.flatView(), this->flatView());
            return *this;
        }
        Matrix& operator-=(const Matrix& b) {
            check(b);
            KokkosBatched::SerialAxpy::invoke(T(-1), b.flatView(), this->flatView());
            return *this;
        }
        /// Add a scalar multiple of the identity; requires a square matrix.
        Matrix& operator+=(const T& a) {
            if (this->nrows() != this->ncols())
                throw ::LogicalError("Matrix diagonal", "Not square");
            for (int i = 0; i < this->nrows(); ++i)
                (*this)(i, i) += a;
            return *this;
        }
        Matrix& operator-=(const T& a) { return *this += -a; }
        Matrix& operator*=(const Matrix& b) { return *this = dotm(b); }
        Matrix transpose() const {
            Matrix r(this->ncols(), this->nrows());
            KokkosBatched::SerialCopy<KokkosBatched::Trans::Transpose>::invoke(
                    this->view(), r.view());
            return r;
        }
        Matrix dotm(const Matrix& b) const {
            if (this->ncols() != b.nrows())
                throw ::LogicalError("Matrix product", "Dimension mismatch");
            Matrix r(this->nrows(), b.ncols(), T(0));
            KokkosBatched::SerialGemm<
                    KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose,
                    KokkosBatched::Algo::Gemm::Unblocked>::
                    invoke(T(1), this->view(), b.view(), T(0), r.view());
            return r;
        }
        Vector<T> dotcv(const Vector<T>& b) const {
            if (this->ncols() != b.size())
                throw ::LogicalError("Matrix-vector", "Dimension mismatch");
            Vector<T> r(this->nrows(), T(0));
            KokkosBlas::Experimental::serial_gemv(
                    'N', T(1), this->view(), b.view(), T(0), r.view());
            return r;
        }
        Vector<T> dotrv(const Vector<T>& b) const {
            if (this->nrows() != b.size())
                throw ::LogicalError("Vector-matrix", "Dimension mismatch");
            Vector<T> r(this->ncols(), T(0));
            KokkosBlas::Experimental::serial_gemv(
                    'T', T(1), this->view(), b.view(), T(0), r.view());
            return r;
        }

    private:
        void check(const Matrix& b) const {
            if (this->nrows() != b.nrows() || this->ncols() != b.ncols())
                throw ::LogicalError("Matrix", "Dimension mismatch");
        }
    };
    template <class T>
    Matrix<T> operator+(Matrix<T> a, const Matrix<T>& b) {
        return a += b;
    }
    template <class T>
    Matrix<T> operator-(Matrix<T> a, const Matrix<T>& b) {
        return a -= b;
    }
    template <class T>
    Matrix<T> operator-(Matrix<T> a) {
        return a *= T(-1);
    }
    template <class T>
    Matrix<T> operator*(const Matrix<T>& a, const Matrix<T>& b) {
        return a.dotm(b);
    }
    template <class T>
    Vector<T> operator*(const Matrix<T>& a, const Vector<T>& b) {
        return a.dotcv(b);
    }
    template <class T>
    Vector<T> operator*(const Vector<T>& a, const Matrix<T>& b) {
        return b.dotrv(a);
    }
    template <class T>
    Matrix<T> operator*(Matrix<T> a, const T& b) {
        return a *= b;
    }
    template <class T>
    Matrix<T> operator*(const T& a, Matrix<T> b) {
        return b *= a;
    }
    template <class T>
    Matrix<T> operator/(Matrix<T> a, const T& b) {
        return a /= b;
    }
}  // namespace tps
#endif
