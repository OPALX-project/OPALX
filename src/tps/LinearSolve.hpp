#ifndef TPS_LINEAR_SOLVE_HPP
#define TPS_LINEAR_SOLVE_HPP
#include "tps/Algebra/Matrix.h"
#include "tps/Error.hpp"
#include <KokkosBatched_Getrf.hpp>
#include <KokkosBatched_Getrs.hpp>
namespace tps {
/** @brief Small host-side LU factorization using Kokkos Kernels partial pivoting.
 * @details No hand-written elimination or triangular solve is retained. INFO
 * is checked, and the original matrix is not modified. Not an eigensolver.
 */
template <class T> class LUFactorization {
  public:
    explicit LUFactorization(const Matrix<T>& matrix) : factors_(matrix), pivots_(matrix.nrows()) {
        if (matrix.nrows() != matrix.ncols())
            throw SizeError("LUFactorization", "Matrix must be square");
        int info = KokkosBatched::SerialGetrf<KokkosBatched::Algo::Getrf::Unblocked>::invoke(
            factors_.view(), pivots_.view());
        if (info > 0)
            throw SingularMatrixError("Kokkos Getrf");
        if (info < 0)
            throw std::runtime_error("Kokkos Getrf invalid argument");
    }
    void solve(Vector<T>& rhs) const {
        if (rhs.size() != factors_.nrows())
            throw SizeError("LUFactorization solve", "Dimension mismatch");
        solveView(rhs.view());
    }
    void solve(Matrix<T>& rhs) const {
        if (rhs.nrows() != factors_.nrows())
            throw SizeError("LUFactorization solve", "Dimension mismatch");
        for (int i = 0; i < rhs.ncols(); ++i)
            solveView(Kokkos::subview(rhs.view(), Kokkos::ALL, i));
    }
    Matrix<T> inverse() const {
        auto r = Matrix<T>::identity(factors_.nrows());
        solve(r);
        return r;
    }

  private:
    template <class View> void solveView(const View& rhs) const {
        int info = KokkosBatched::SerialGetrs<
            KokkosBatched::Trans::NoTranspose,
            KokkosBatched::Algo::Getrs::Unblocked>::invoke(factors_.view(), pivots_.view(), rhs);
        if (info)
            throw std::runtime_error("Kokkos Getrs failed");
    }
    Matrix<T> factors_;
    Array1D<int> pivots_;
};
} // namespace tps
#endif
