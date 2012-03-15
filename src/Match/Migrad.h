#ifndef OPAL_Migrad_HH
#define OPAL_Migrad_HH 1

// ------------------------------------------------------------------------
// $RCSfile: Migrad.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Migrad
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"
#include "Match/MatchState.h"

template <class T> class Matrix;
template <class T> class Vector;


// Class Migrad
// ------------------------------------------------------------------------
/// The MIGRAD command.
//  This class encapsulates a minimisation according to the MIGRAD method,
//  Minimization by a gradient method due to Davidon/Fletcher/Powell
//  (Computer Journal 13, 317, 1970).  Also described in
//  [center]
//    J. F. Bonnans et al., Optimisation Numerique, pp. 44-46.
//    Springer, Berlin, 1997
//  [/center]
//  Algorithm rewritten following the MINUIT package.

class Migrad: public Action {

public:

  /// Exemplar constructor.
  Migrad();

  virtual ~Migrad();

  /// Make clone.
  virtual Migrad *clone(const string &name);
  
  /// Execute the command.
  virtual void execute();

private:

  // Not implemented.
  Migrad(const Migrad &);
  void operator=(const Migrad &);

  // Clone constructor.
  Migrad(const string &name, Migrad *parent);

  // Find covariance matrix.
  void hessenberg(Vector<double> &X, Vector<double> &F,
		  Matrix<double> &V, Vector<double> &G, Vector<double> &G2);

  // Find first and second derivatives.
  void derivatives(Vector<double> &X, Vector<double> &F,
		   Vector<double> &G, Vector<double> &G2);

  // Search for minimum along a line.
  MatchState lineSearch
  (Vector<double> &X, Vector<double> &dX, Vector<double> &F, double tol);

  // Force symmetric matrix positive definite.
  void forcePositiveDefinite(Matrix<double> &A);

  // Invert symmetric matrix.
  bool invertSymmetric(Matrix<double> &A);

  // Find real eigenvalues of a symmetric matrix.
  int symmetricEigen(const Matrix<double> &A, Vector<double> &eigen);

  // The sum of squares of the matching functions (local copy).
  double fmin;
};

#endif // OPAL_Migrad_HH
