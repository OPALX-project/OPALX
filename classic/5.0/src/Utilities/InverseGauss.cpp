// ------------------------------------------------------------------------
// $RCSfile: InverseGauss.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Inverse Gaussian.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2001/08/24 19:33:44 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Utilities/InverseGauss.h"
#include "Utilities/DomainError.h"
#include <cmath>


//@ Inverse Gaussian density function.
double InverseGauss(double p0)
{
  static const double a0 =  2.5050240;
  static const double a1 =  2.7724523;
  static const double a2 =  2.1635544;
  static const double a3 =  4.5585614e+01;
  static const double b0 =  2.2157257;
  static const double b1 = -1.4946512e+01;
  static const double b2 =  7.9731883e+01;
  static const double b3 = -2.7713713e+02;
  static const double b4 =  4.0314354e+02;
  static const double c0 =  2.4335936;
  static const double c1 = -2.4719139e+01;
  static const double c2 =  2.4648581e+02;
  static const double c3 = -1.5585873e+03;
  static const double c4 =  4.1394487e+03;
  static const double d0 =  2.6346872;
  static const double d1 = -4.1028898e+01;
  static const double d2 =  7.4942805e+02;
  static const double d3 = -8.5400893e+03;
  static const double d4 =  4.0895693e+04;
  static const double e0 =  2.8224654;
  static const double e1 = -6.8317697e+01;
  static const double e2 =  2.2566998e+03;
  static const double e3 = -4.6004775e+04;
  static const double e4 =  3.9399134e+05;
  static const double f0 = -8.1807613e-02;
  static const double f1 = -2.8358733;
  static const double f2 =  1.4902469;
  static const double pp1 = 0.334624883253;
  static const double qq2 = 0.090230446775;
  static const double qq3 = 0.049905685242;
  static const double qq4 = 0.027852994157;
  static const double qq5 = 0.015645650215;

  double p = p0 - 0.5;
  double p1 = std::abs(p);
  if (p1 < pp1) {
    double p2 = p*p;
    return (((a3 * p2 + a2) * p2 + a1) * p2 + a0) * p;
  } else {
    double gauinv;
    double q = 0.5 - p1;

    if (q > qq2) {
      gauinv = (((b4 * q + b3) * q + b2) * q + b1) * q + b0;
    } else if (q > qq3) {
      gauinv = (((c4 * q + c3) * q + c2) * q + c1) * q + c0;
    } else if (q > qq4) {
      gauinv = (((d4 * q + d3) * q + d2) * q + d1) * q + d0;
    } else if (q > qq5) {
      gauinv = (((e4 * q + e3) * q + e2) * q + e1) * q + e0;
    } else if (q > 0.0) {
      double t = sqrt(- 2.0 * log(q));
      gauinv = t + f0 + f1 / (f2 + t);
    } else {
      throw DomainError("InverseGauss()");
    }

    if (p < 0.0) gauinv = - gauinv;
    return gauinv;
  }
}
