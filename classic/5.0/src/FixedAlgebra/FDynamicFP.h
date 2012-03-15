#ifndef MAD_FDynamicFP_HH
#define MAD_FDynamicFP_HH

// ------------------------------------------------------------------------
// $RCSfile: FDynamicFP.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: FDynamicFP
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:50:59 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/FVps.h"

template <class T, int N> class FVector;


// Class FDynamicFP
// ------------------------------------------------------------------------
/// Dynamic fixed point of a Truncated power series map.

template <int N>
class FDynamicFP {

public:

  /// Constructor.
  //  Find fixed point of [b]map[/b].
  FDynamicFP(const FVps<double,N> &map);

  FDynamicFP();
  FDynamicFP(const FDynamicFP &);
  ~FDynamicFP();

  /// Get the transformation to the fixed point.
  const FVector<double,N> &getFixedPoint() const;

  /// Get the map around the fixed point.
  const FVps<double,N> &getFixedPointMap() const;

private:

  // Not implemented.
  void operator=(const FDynamicFP &);

  // Fixed point position.
  FVector<double,N> fixedPoint;

  // Map around the fixed point.
  FVps<double,N> fixedPointMap;
};


#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FLUMatrix.h"
#include "FixedAlgebra/FVector.h"
#include "FixedAlgebra/FVps.h"


// Class FDynamicFP
// ------------------------------------------------------------------------

template <int N>
FDynamicFP<N>::FDynamicFP():
  fixedPoint(), fixedPointMap()
{}


template <int N>
FDynamicFP<N>::FDynamicFP(const FDynamicFP &rhs):
  fixedPoint(rhs.fixedPoint), fixedPointMap(rhs.fixedPointMap)
{}


template <int N>
FDynamicFP<N>::FDynamicFP(const FVps<double,N> &map):
  fixedPoint(), fixedPointMap(map)
{
  FVps<double,N> identity;

  // Determine fixPoint of map T by iterating the contraction map
  //   z_{n+1} = z_n - (M(z_n) - I)^{-1} . (T(z_n) - z_n)
  // Here M(z_n) is the linear part of T(z_n + z), and z_0 = 0.
  if (map.minOrder() == 0) {
    while (true) {
      // Get system of equations for fixed point.
      FMatrix<double,N,N> A   = fixedPointMap.linearTerms();
      FVector<double,N> Error = fixedPointMap.constantTerm();
      double error = 0.0;
  
      for (int i = 0; i < N; i++) {
        A(i,i) -= 1.0;
        if (abs(Error(i)) > error) error = abs(Error(i));
      }
  
      // Test for convergence.
      static const double tol = 1.0e-10;
      if (error < tol) break;
      
      // Correction for fixed point.
      FLUMatrix<double,N> lu(A);
      lu.backSubstitute(Error);
      fixedPoint -= Error;
      
      // Build map around fixed point found so far.
      fixedPointMap = map.substitute(identity + fixedPoint) - fixedPoint;
    }
  }
}


template <int N>
FDynamicFP<N>::~FDynamicFP()
{}


template <int N>
const FVector<double,N> &FDynamicFP<N>::getFixedPoint() const
{
  return fixedPoint;
}


template <int N>
const FVps<double,N> &FDynamicFP<N>::getFixedPointMap() const
{
  return fixedPointMap;
}

#endif // MAD_FDynamicFP_HH
