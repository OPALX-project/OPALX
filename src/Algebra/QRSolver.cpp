// ------------------------------------------------------------------------
// $RCSfile: QRSolver.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: QRSolver
//   This class implements routines for the least-square solution of
//   systems of linear equations.
//   Given an  m  by  n  matrix  A, an  n  by  n  diagonal matrix  D,
//   and an m-vector  B, two problem can be solved:
//
//   1. Solve the the system
//           A*X = B,
//      in the least squares sense. The first step to solve this problem is
//           QRSolver solver(A, pivot);
//      The second step is then
//           solver.solveR(X);
//
//   2. Solve the the two systems
//           A*X = B,     D*X = 0,
//      in the least squares sense. The first step to solve this problem is
//           QRSolver solver(A, pivot);
//      The second step is
//           solver.solveS(D, X);
//      it can be repeated as many times as required for different diagonal
//      matrices D.
//
//   In both cases, the method
//           solver.getColNorm(C);
//   can be called to return the original column norms of A.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:35 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Algebra/QRSolver.h"
#include "Algebra/Matrix.h"
#include "Algebra/Vector.h"
#include <cassert>
#include <cmath>


// This method uses Householder transformations with optional column
// pivoting to compute a QR-factorization of the m-by-n matrix  A.
// QRSolver determines an orthogonal matrix  Q, a permutation matrix  P,
// and an upper trapezoidal matrix  R  with diagonal elements of
// nonincreasing magnitude, such that  A*P = Q*R.  The Householder
// transformation for column  k, k = 1, 2, ..., min(m,n),  is of the form
//
//                         T
//         I - (1/U(k))*U*U
//
// where  U  has zeros in the first  k-1  positions. The form of this
// transformation and the method of pivoting first appeared in the
// corresponding LINPACK subroutine.
//
QRSolver::QRSolver
(const Matrix<double> &A, const Vector<double> &B, bool pivot):
    n(A.ncols()), m(A.nrows()), Q(A), R(n, n, 0.0), QtB(B), S(n, n, 0.0),
    column_norm(n, 0.0), piv_col(n, 0) {
    assert(B.size() >= m);
    Vector<double> dot(n);

    // Compute the initial and internal column norms.
    for(int j = 0; j < n; ++j) {
        double temp = 0.0;
        for(int i = 0; i < m; ++i) temp += Q(i, j) * Q(i, j);
        dot(j) = column_norm(j) = sqrt(temp);
        piv_col(j) = j;
    }

    // Reduce  A  to  R  using Householder transformations.
    // Up to  min(m, n)  columns are used.
    for(int j = 0; j < std::min(m, n); ++j) {
        if(pivot) {
            // Bring column of largest norm into pivot position.
            int kmax = j;

            for(int k = j + 1; k < n; ++k) {
                if(dot(k) > dot(kmax)) kmax = k;
            }

            if(kmax != j) {
                std::swap_ranges(Q.col_begin(j), Q.col_end(j), Q.col_begin(kmax));
                dot(kmax) = dot(j);
                std::swap(piv_col(j), piv_col(kmax));
            }
        }

        // Compute the Householder transformation which reduces the j-th
        // column of  A  to a multiple of the j-th unit vector.
        double unorm = 0.0;
        for(int k = j; k < m; ++k) unorm += Q(k, j) * Q(k, j);
        unorm = (Q(j, j) > 0.0) ? sqrt(unorm) : - sqrt(unorm);

        // Build column   j  of  R.
        for(int i = 0; i < j; ++i) {
            R(i, j) = Q(i, j);
            Q(i, j) = 0.0;
        }
        R(j, j) = - unorm;

        if(unorm != 0.0) {
            // Modify column  j  of  Q.
            for(int i = j; i < m; ++i) {
                Q(i, j) /= unorm;
            }
            Q(j, j) += 1.0;

            // Transform remaining columns of  Q.
            for(int k = j + 1; k < n; ++k) {
                // Apply transformation to column  k.
                double sum = 0.0;
                for(int i = j; i < m; ++i) sum += Q(i, j) * Q(i, k);
                sum /= Q(j, j);
                for(int i = j; i < m; ++i) Q(i, k) -= sum * Q(i, j);

                if(pivot) {
                    // Recompute column norm.
                    double sum = 0.0;
                    for(int i = j + 1; i < m; ++i) sum += Q(i, k) * Q(i, k);
                    dot(k) = sqrt(sum);
                }
            }

            // Transform r.h.s.
            if(Q(j, j) != 0.0) {
                double sum = 0.0;
                for(int i = j; i < m; ++i) sum += Q(i, j) * QtB(i);
                sum /= Q(j, j);
                for(int i = j; i < m; ++i) QtB(i) -= Q(i, j) * sum;
            }
        }
    }
}


QRSolver::~QRSolver()
{}


void QRSolver::solveR(Vector<double> &X) const {
    solve(R, QtB, X);
}


void QRSolver::solveS(const Array1D<double> &D, double p, Vector<double> &X) {
    Vector<double> Z(m);
    Array1D<double> M(D);
    double root = sqrt(p);
    for(int i = 0; i < n; ++i) M(i) *= root;

    R_to_S(R, M, S, Z);
    solve(S, Z, X);
}


void QRSolver::solveRT(Vector<double> &V) const {
    solveT(R, V);
}


void QRSolver::solveST(Vector<double> &V) const {
    solveT(S, V);
}


void QRSolver::getColNorm(Array1D<double> &norms) const {
    norms = column_norm;
}


void QRSolver::R_to_S(const Matrix<double> &R, const Array1D<double> &D,
                      Matrix<double> &S, Vector<double> &Z) {
    // Eliminate the diagonal matrix  D using Givens rotations.
    assert(D.size() >= n);
    S = R;
    Z = QtB;

    for(int j = 0; j < n; ++j) {
        // Prepare the row of  D  to be eliminated, locating the diagonal
        // element using  P  from the QR-factorization.
        int l = piv_col(j);

        if(D(l) != 0.0) {
            Vector<double> row_D(n, 0.0);
            row_D(j) = D(l);

            // The transformations to eliminate the row of  D; modify only
            // a single element of  transpose(Q)*B  beyond the first n,
            // which is initially 0.
            double QtBpj = 0.0;

            for(int k = j; k < n; ++k) {
                // Determine a Givens rotation which eliminates the appropriate
                // element in the current row of D.
                if(row_D(k) != 0.0) {
                    double t = sqrt(S(k, k) * S(k, k) + row_D(k) * row_D(k));
                    double s = row_D(k) / t;
                    double c = S(k, k) / t;

                    // Compute the modified diagonal element of S and the modified
                    // element of ((Q transpose)*B, 0).
                    S(k, k) = t;
                    t =  c * Z(k) + s * QtBpj;
                    QtBpj = -s * Z(k) + c * QtBpj;
                    Z(k) = t;

                    // Accumulate the transformation in the row of S.
                    for(int i = k + 1; i < n; ++i) {
                        t        =  c * S(k, i) + s * row_D(i);
                        row_D(i) = -s * S(k, i) + c * row_D(i);
                        S(k, i)   = t;
                    }
                }
            }
        }
    }
}


void QRSolver::solve(const Matrix<double> &R, const Vector<double> &QtB,
                     Vector<double> &X) const {
    // Copy QtB to the solution.
    Vector<double> Z(QtB);

    // Solve the triangular system for Z.
    // If R is singular, then obtain a least squares solution.
    int rank = n;
    for(int j = 0; j < n; ++j) {
        if(R(j, j) == 0.0  &&  rank == n) rank = j;
        if(rank < n) Z(j) = 0.0;
    }

    for(int j = rank; j-- > 0;) {
        double sum = Z(j);
        for(int i = j + 1; i < rank; ++i) {
            sum -= R(j, i) * Z(i);
        }
        Z(j) = sum / R(j, j);
    }

    // Permute the components of Z back to components of X.
    X = Vector<double>(n, 0.0);
    for(int j = 0; j < n; ++j) {
        X(piv_col(j)) = Z(j);
    }
}


void QRSolver::solveT(const Matrix<double> &R, Vector<double> &V) const {
    assert(V.size() == n);
    Vector<double> X(n);
    for(int i = 0; i < n; ++i) X(i) = V(piv_col(i));

    for(int i = 0; i < n; ++i) {
        if(R(i, i) == 0.0) break;
        for(int j = 0; j < i; ++j) X(j) -= R(i, j) * V(i);
        X(i) /= R(i, i);
    }

    V = X;
}
