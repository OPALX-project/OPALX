#ifndef OPAL_QRSolver_HH
#define OPAL_QRSolver_HH 1

// ------------------------------------------------------------------------
// $RCSfile: QRSolver.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: QRSolver
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:35 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Algebra/Matrix.h"
#include "Algebra/Vector.h"
#include <cstddef>


// Class QRSolver
// ------------------------------------------------------------------------
/// Least-square solution of systems of linear equations.
//  Given an  m  by  n  matrix  A, an  n  by  n  diagonal matrix  D,
//  and an m-vector  B, two problem can be solved:
//  [ol][li]
//  Solve the the system $A*X = B$ in the least squares sense.
//  The first step to solve this problem is:
//  [pre]
//       QRSolver solver(A, pivot);
//  [/pre]
//  The second step is then
//  [pre]
//       solver.solveR(X);
//  [/pre]
//  [li]
//  Solve the the two systems $A*X = B, D*X = 0$ in the least squares sense.
//  The first step to solve this problem is
//  [pre]
//       QRSolver solver(A, pivot);
//  [/pre]
//      The second step is
//  [pre]
//       solver.solveS(D, X);
//  [/pre]
//  The second step can be repeated as many times as required for different
//  diagonal matrices $D$.
//  [/ol]
//  In both cases, the method
//  [pre]
//       solver.getColNorm(C);
//  [/pre]
//  can be called to return the original column norms of $A$.

class QRSolver {

public:

    /// Constructor.
    //  Determine the QR-factorisation $A = Q*R$ of the matrix $A$
    //  and transform the right-hand side $B$. If [b]pivot[/b] is
    //  true, then pivot search is done.
    //  [p]
    //  This method uses Householder transformations with optional column
    //  pivoting to compute a QR-factorization of the m-by-n matrix  $A$.
    //  It determines an orthogonal matrix $Q$, a permutation matrix $P$,
    //  and an upper trapezoidal matrix $R$ with diagonal elements of
    //  nonincreasing magnitude, such that $A*P = Q*R$. The Householder
    //  transformation for column $k$ is of the form $I - (1/U(k))*U*U^T$,
    //  where $U$ has zeros in the first $k-1$ positions.
    //  The form of this transformation and the method of pivoting first
    //  appeared in the corresponding LINPACK subroutine.
    QRSolver(const Matrix<double> &A, const Vector<double> &B, bool pivot);

    ~QRSolver();

    /// Solution of $A*X = B$ in the least-squares sense.
    void solveR(Vector<double> &X) const;

    /// Solution of $A*X = B, D*X = 0$ in the least-squares sense.
    void solveS(const Array1D<double> &D, double p, Vector<double> &X);

    /// Pre-multiply the vector  $V$  by  $R.transpose()^{-1}$.
    void solveRT(Vector<double> &V) const;

    /// Pre-multiply the vector  $V$  by  $S.transpose()^{-1}$.
    // Requires prior execution of solveS.
    void solveST(Vector<double> &V) const;

    /// Return the original column norms of the matrix  $A$.
    void getColNorm(Array1D<double> &) const;

private:

    // Transform  R  to  S,  given the diagonal matrix D.
    void R_to_S(const Matrix<double> &R,
                const Array1D<double> &D,
                Matrix<double> &S,
                Vector<double> &Z);

    // Solution of  R*X = QtB.
    void solve(const Matrix<double> &R,
               const Vector<double> &QtB,
               Vector<double> &X) const;

    // Solution of R.transpose()**(-1) * V
    void solveT(const Matrix<double> &R, Vector<double> &V) const;

    // The dimensions of the problem.
    int n;   // number of columns.
    int m;   // number of rows.

    // The lower trapezoidal part of the matrix  Q  contains the  n  vectors
    // defining a factored form of the orthogonal matrix  Q.
    Array2D<double> Q;

    // The upper-trapezoidal matrix  R  contains the matrix  R.
    Matrix<double> R;

    // The transformed r.h.s. Q.transpose() * B.
    Vector<double> QtB;

    // The upper-trapezoidal matrix  S  contains the matrix  S.
    Matrix<double> S;

    // The vector  column_norm  contains the original column norms of  M.
    Array1D<double> column_norm;

    // The vector  piv_col  defines the permutation matrix  P.
    // Column j of  P  is column  piv_col(j)  of the identity matrix.
    Array1D<int> piv_col;
};

#endif // OPAL_QRSolver_HH
