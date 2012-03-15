// ------------------------------------------------------------------------
// $RCSfile: DynamicFixedPoint.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DynamicFixedPoint
//   Find dynamic fixed point of a truncated power series map.
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algebra/DynamicFixedPoint.h"
#include "Algebra/Matrix.h"
#include "Algebra/LUMatrix.h"
#include "Algebra/Vector.h"
#include "Algebra/VpsInvMap.h"
#include <cmath>


// Tolerance for fixed point search.
namespace {
  const double tol = 1.0e-10;
}


// Class DynamicFixedPoint
// ------------------------------------------------------------------------

DynamicFixedPoint::DynamicFixedPoint():
  fixedPoint(),
  fixedPointMap()

{}


DynamicFixedPoint::DynamicFixedPoint(const DynamicFixedPoint &rhs):
  fixedPoint(rhs.fixedPoint),
  fixedPointMap(rhs.fixedPointMap)
{}


DynamicFixedPoint::DynamicFixedPoint(const VpsInvMap<double> &map):
  fixedPoint(Vector<double>(map.getDimension(), 0.0)),
  fixedPointMap(map)
{
  int nDim = map.getDimension();
  VpsInvMap<double> ident = VpsInvMap<double>::identity(nDim);

  while (true) {
    // Get system of equations for fixed point.
    Matrix<double> A     = fixedPointMap.linearTerms();
    Vector<double> Error = fixedPointMap.constantTerm();
    double error = 0.0;

    for (int i = 0; i < nDim; i++) {
      A(i,i) -= 1.0;
      if (abs(Error(i)) > error) error = abs(Error(i));
    }

    // Test for convergence.
    if (error < tol) break;
    
    // Correction for fixed point.
    LUMatrix<double> lu(A);
    lu.backSubstitute(Error);
    fixedPoint -= Error;
    
    // Build map around fixed point found so far.
    fixedPointMap = map.substitute(ident + fixedPoint) - fixedPoint;
  }
}


DynamicFixedPoint::~DynamicFixedPoint()
{}


const Vector<double> &DynamicFixedPoint::getFixedPoint() const
{
  return fixedPoint;
}


const VpsInvMap<double> &DynamicFixedPoint::getFixedPointMap() const
{
  return fixedPointMap;
}

