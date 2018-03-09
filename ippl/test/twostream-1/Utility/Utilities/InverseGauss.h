// ------------------------------------------------------------------------
// $RCSfile: InverseGauss.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Inverse Gaussian density function.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------


//: Inverse Gaussian density function.
//  The inverse function of
//  <br>
//  y = Phi(x) = 2 * integral exp(-x^2) * dx / sqrt(pi)
//  <br>
//  [br]
//  Argument: A value $y$ between -1 and +1.
//  [br]
//  Return value: the value $x$.

extern double InverseGauss(double);
