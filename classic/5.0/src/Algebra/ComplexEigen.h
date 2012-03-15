#ifndef CLASSIC_ComplexEigen_HH
#define CLASSIC_ComplexEigen_HH

// ------------------------------------------------------------------------
// $RCSfile: ComplexEigen.h,v $
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


// Class ComplexEigen
// ------------------------------------------------------------------------
/// Complex eigenvector routines.
//  Representation of eigenvalues and eigenvectors
//  for a matrix of type Matrix<complex<double> >.

class ComplexEigen {

public:

  /// Constructor.
  //  Construct the vector of eigenvalues for the matrix [b]M[/b].
  //  If [b]vec[/b] is true, the matrix of eigenvectors is also built.
  ComplexEigen(const Matrix<complex<double> > &m, bool vec=false);

  ComplexEigen();
  ComplexEigen(const ComplexEigen &);
  ~ComplexEigen();

  /// Get eigenvalues.
  //  Return the eigenvalues as a complex vector.
  const Vector<complex<double> > &eigenValues() const;

  /// Get eigenvectors.
  //  Return the eigenvectors as the column vectors of a complex matrix.
  const Matrix<complex<double> > &eigenVectors() const;

private:

  // Not implemented.
  void operator=(const ComplexEigen &);

  // Used by eigenvalue and eigenvector routines
  static void balance
  (Matrix<complex<double> > &, int &low, int &high, Array1D<double> &);

  void balbak(int low, int high, const Array1D<double> &scale);

  static void exchange(Matrix<complex<double> > &, int j, int m,
		       int low, int high);

  int hqr(Matrix<complex<double> > &, int low, int high);

  int hqr2(Matrix<complex<double> > &, int low, int high,
	      Array1D<complex<double> > &ort);

  static void orthes(Matrix<complex<double> > &, int low, int high,
		     Array1D<complex<double> > &ort);

  // Representation of the eigenvalues and eigenvectors.
  Vector<complex<double> > lambda;
  Matrix<complex<double> > vectors;
};

#endif // CLASSIC_ComplexEigen_HH
