#ifndef TPS_LINEAR_SOLVE_HPP
#define TPS_LINEAR_SOLVE_HPP
#include "tps/Algebra/Matrix.h"
#include "tps/Error.hpp"
#include <KokkosBatched_Getrf.hpp>
#include <KokkosBatched_Getrs.hpp>
#include <algorithm>
#include <cmath>
namespace tps {
/** @brief Small host-side LU factorization using Kokkos Kernels partial pivoting.
 * @details No hand-written elimination or triangular solve is retained. INFO
 * is checked, and the original matrix is not modified. Not an eigensolver.
 */
template <class T> class LUFactorization {
  public:
    /// Reject exact singularity and U diagonals with magnitude strictly below minimumPivot.
    /// The absolute floor is in matrix-entry units; zero (default) adds no floor.
    /// LU rounding/pivot ties may differ from other elimination algorithms.
    explicit LUFactorization(const Matrix<T>& matrix, double minimumPivot = 0)
        : factors_(matrix), pivots_(matrix.nrows()) {
        if (matrix.nrows() != matrix.ncols())
            throw SizeError("LUFactorization", "Matrix must be square");
        if (!std::isfinite(minimumPivot) || minimumPivot < 0)
            throw std::invalid_argument("LU minimum pivot must be finite and nonnegative");
        for (int i = 0; i < matrix.nrows(); ++i)
            for (int j = 0; j < matrix.ncols(); ++j)
                if (!std::isfinite(std::abs(matrix(i, j))))
                    throw std::invalid_argument("LU matrix must be finite");
        int info = KokkosBatched::SerialGetrf<KokkosBatched::Algo::Getrf::Unblocked>::invoke(
            factors_.view(), pivots_.view());
        if (info > 0)
            throw SingularMatrixError("Kokkos Getrf");
        if (info < 0)
            throw std::runtime_error("Kokkos Getrf invalid argument");
        for (int i = 0; i < factors_.nrows(); ++i) {
            const auto magnitude = std::abs(factors_(i, i));
            if (!std::isfinite(magnitude))
                throw std::runtime_error("Nonfinite LU factor");
            if (magnitude == 0 || magnitude < minimumPivot)
                throw SingularMatrixError("LU minimum pivot");
        }
    }
    void solve(Vector<T>& rhs) const {
        if (rhs.size() != factors_.nrows())
            throw SizeError("LUFactorization solve", "Dimension mismatch");
        solveView<KokkosBatched::Trans::NoTranspose>(rhs.view());
    }
    void solve(Matrix<T>& rhs) const {
        if (rhs.nrows() != factors_.nrows())
            throw SizeError("LUFactorization solve", "Dimension mismatch");
        for (int i = 0; i < rhs.ncols(); ++i)
            solveView<KokkosBatched::Trans::NoTranspose>(
                Kokkos::subview(rhs.view(), Kokkos::ALL, i));
    }
    /// Solve A^T x=b in place with the original factorization/pivots (no conjugation).
    void solveTranspose(Vector<T>& rhs) const {
        if (rhs.size() != factors_.nrows())
            throw SizeError("LU transpose solve", "Dimension mismatch");
        solveView<KokkosBatched::Trans::Transpose>(rhs.view());
    }
    /// Solve A^T X=B in place, one RHS per column.
    void solveTranspose(Matrix<T>& rhs) const {
        if (rhs.nrows() != factors_.nrows())
            throw SizeError("LU transpose solve", "Dimension mismatch");
        for (int i = 0; i < rhs.ncols(); ++i)
            solveView<KokkosBatched::Trans::Transpose>(Kokkos::subview(rhs.view(), Kokkos::ALL, i));
    }
    /// Product of U diagonals and row-swap parity; empty determinant is one.
    /// Direct product, not a log determinant: overflow/underflow follows scalar arithmetic.
    T determinant() const {
        T result(1);
        for (int i = 0; i < factors_.nrows(); ++i) {
            if (pivots_[i] != i)
                result = -result; // Getrf uses zero-based swap indices.
            result *= factors_(i, i);
        }
        return result;
    }
    Matrix<T> inverse() const {
        auto r = Matrix<T>::identity(factors_.nrows());
        solve(r);
        return r;
    }

  private:
    template <class Trans, class View> void solveView(const View& rhs) const {
        int info = KokkosBatched::SerialGetrs<Trans, KokkosBatched::Algo::Getrs::Unblocked>::invoke(
            factors_.view(), pivots_.view(), rhs);
        if (info)
            throw std::runtime_error("Kokkos Getrs failed");
    }
    Matrix<T> factors_;
    Array1D<int> pivots_;
};
/** @brief Determinant through pivoted Kokkos LU; exact singularity returns zero.
 * @details No absolute pivot floor is applied. Other validation/backend failures
 * propagate, rather than being mistaken for singularity. Direct scalar product.
 */
template <class T> T determinant(const Matrix<T>& matrix) {
    try {
        return LUFactorization<T>(matrix).determinant();
    } catch (const SingularMatrixError&) {
        return T(0);
    }
}

/** @brief Maximum absolute row sum, including the original matrix's physical scaling. */
template <class T> auto infinityNorm(const Matrix<T>& matrix) {
    using Magnitude = decltype(std::abs(T{}));
    Magnitude maximum(0);
    for (int i = 0; i < matrix.nrows(); ++i) {
        Magnitude sum(0);
        for (int j = 0; j < matrix.ncols(); ++j)
            sum += std::abs(matrix(i, j));
        // Propagate NaN as well as overflow instead of reporting a smaller finite norm.
        if (!std::isfinite(sum))
            return sum;
        maximum = std::max(maximum, sum);
    }
    return maximum;
}
} // namespace tps
#endif
