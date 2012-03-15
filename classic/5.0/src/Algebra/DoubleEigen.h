#ifndef CLASSIC_DoubleEigen_HH
#define CLASSIC_DoubleEigen_HH

// ------------------------------------------------------------------------
// $RCSfile: DoubleEigen.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class ComplexEigen:
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algebra/Matrix.h"
#include "Algebra/Vector.h"
#include <complex>

using std::complex;


// Class DoubleEigen
// ------------------------------------------------------------------------
/// Double eigenvector routines.
//  Representation of eigenvalues and eigenvectors
//  for a matrix of type Matrix<double>.

class DoubleEigen {

public:

  /// Constructor.
  //  Construct the vector of complex eigenvalues for the matrix [b]M[/b].
  //  If [b]vec[/b] is true, the matrix of eigenvectors is also built.
  DoubleEigen(const Matrix<double> &m, bool vec=false);

  DoubleEigen();
  DoubleEigen(const DoubleEigen &);
  ~DoubleEigen();

  /// Get eigenvalues.
  //  Return the eigenvalues as a complex vector.
  Vector<complex<double> > eigenValues() const;

  /// Get eigenvectors.
  //  Return the eigenvectors as the column vectors of a complex matrix.
  Matrix<complex<double> > eigenVectors() const;

  /// Get eigenvectors.
  //  Return the eigenvectors packed in a real matrix.
  //  Real eigenvectors appear a single column in the same position as
  //  the corresponding real eigenvalue.
  //  Complex eigenvalues occur in conjugate pairs, and the corresponding
  //  eigenvectors appear as two real columns containing the real and 
  //  imaginary parts respectively.
  Matrix<double> packedEigenVectors() const;

private:

  // Not implemented.
  void operator=(const DoubleEigen &);

  // Used by eigenvalue and eigenvector routines
  static void balance(Matrix<double> &, int &low, int &high,
		      Vector<double> &);

  static void exchange(Matrix<double> &, int j, int m,
		       int low, int high);

  static void elmhes(Matrix<double> &, int low, int high,
		     Array1D<int> &index);

  void elmtran(Matrix<double> &, int low, int high,
	       Array1D<int> &index);

  int hqr(Matrix<double> &, int low, int high);

  int hqr2(Matrix<double> &, int low, int high);

  void balbak(int low, int high, Vector<double> &scale);

  // Representation of the eigenvalues and eigenvectors.
  Vector<complex<double> > lambda;
  Matrix<double>           vectors;
};

#endif // CLASSIC_DoubleEigen_HH
