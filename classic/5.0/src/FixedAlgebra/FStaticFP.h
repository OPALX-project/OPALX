#ifndef MAD_FStaticFP_HH
#define MAD_FStaticFP_HH

// ------------------------------------------------------------------------
// $RCSfile: FStaticFP.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: FStaticFP
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2003/11/07 18:03:01 $
// $Author: dbruhwil $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/FVps.h"
#include <cmath>


// Class FStaticFP
// ------------------------------------------------------------------------
/// Static fixed point of a Truncated power series map.

template <int N>
class FStaticFP {

public:

  /// Constructor.
  //  Find fixed point for [b]map[/b].
  FStaticFP(const FVps<double,N> &map);

  FStaticFP();
  FStaticFP(const FStaticFP &rhs);
  ~FStaticFP();

  /// Get the dispersion map.
  const FVps<double,N> &getDispersion() const;

  /// Get the transformation to the fixed point.
  const FVps<double,N> &getFixedPoint() const;

  /// Get the map around the fixed point.
  const FVps<double,N> &getFixedPointMap() const;

private:

  // Not implemented.
  void operator=(const FStaticFP &);

  // Fixed point position.
  FVps<double,N> fixedPoint;

  // Map around the fixed point.
  FVps<double,N> fixedPointMap;
};


#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FLUMatrix.h"
#include "FixedAlgebra/FMonomial.h"
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FVector.h"


// Class FStaticFP
// ------------------------------------------------------------------------

template <int N>
FStaticFP<N>::FStaticFP():
  fixedPoint(), fixedPointMap()
{}


template <int N>
FStaticFP<N>::FStaticFP(const FStaticFP &fp):
  fixedPoint(fp.fixedPoint), fixedPointMap(fp.fixedPointMap)
{}


template <int N>
FStaticFP<N>::FStaticFP(const FVps<double,N> &map):
  fixedPoint(), fixedPointMap()
{
  //std::cerr << "==> In FStaticFP<N>::FStaticFP(const FVps<double,N> &map)..." << std::endl;
  const int nFree = 4;
  FLUMatrix<double,N> lu;
  FVector<double,N> fixPoint;

  // Find fixed point;
  // tempMap will become the map around the fixed point.
  FVps<double,N> tempMap(map);

  // Determine fixPoint of map T by iterating the contraction map
  //   z_{n+1} = z_n - (M(z_n) - I)^{-1} . (T(z_n) - z_n)
  // Here M(z_n) is the matrix part of T(z_n + z); and the iterating atarts at z_0 = 0.
  if (tempMap.getMinOrder() == 0) {
    while (true) {
      // Get system of equations for fixed point.
      FMatrix<double,N,N>   A  = tempMap.linearTerms();
      FVector<double,N>  Error = tempMap.constantTerm();
      double error = 0.0;
        
      // Form (A - I).
      for (int i = 0; i < nFree; i++) {
        A(i,i) -= 1.0;
        if (abs(Error(i)) > error) error = abs(Error(i));
      }
      for (int i = nFree; i < N; i++) {
        for (int j = 0; j < N; j++) A(i,j) = A(j,i) = 0.0;
        A(i,i) = 1.0;
        Error(i) = 0.0;
      }
      
      // Correct the fixed point.
      lu = FLUMatrix<double,N>(A);
      lu.backSubstitute(Error);
      fixPoint -= Error;
  
      // Test for convergence.
      static const double tol = 1.0e-10;
      if (error < tol) break;
      
      // Build map around fixed point found so far.
      // (Note that inside this while loop fixedPoint == Identity.)
      tempMap = map.substitute(fixPoint + fixedPoint) - fixPoint;
    }
  } else {
    // Force computation of lu.
    FMatrix<double,N,N> A = tempMap.linearTerms();
    for (int i = 0; i < nFree; i++) A(i,i) -= 1.0;
    for (int i = nFree; i < N; i++) {
      for (int j = 0; j < N; j++) A(i,j) = A(j,i) = 0.0;
      A(i,i) = 1.0;
    }
    lu = FLUMatrix<double,N>(A);
  }

  // Fixed point map for orders one and higher in delta (dispersion).
  fixedPoint.setTruncOrder(map.getMaxOrder());
  for (int limit = 1; limit <= map.getMaxOrder(); limit++) {
    // Remove pure delta terms.
    FVps<double,N> Q = tempMap.substitute(fixedPoint) -
                         fixedPoint.substitute(fixedPointMap);
    FMonomial<N> delta;
    delta[nFree+1] = limit;
    FVector<double,N> Error;
    for (int i = 0; i < nFree; i++) {
      fixedPointMap[i] += Q[i].filter(limit, limit);
      Error(i) = Q[i].getCoefficient(delta);
      fixedPointMap[i].setCoefficient(delta, 0.0);
    }

    lu.backSubstitute(Error);

    for (int i = 0; i < nFree; i += 2) {
      // Effect on orbit.
      fixedPoint[i].setCoefficient(delta, - Error(i));
      fixedPoint[i+1].setCoefficient(delta, - Error(i+1));

      // Orbit terms acting on path length difference.
      FMonomial<N> powers;
      powers[nFree+1] = limit - 1;
      powers[i] = 1;
      fixedPoint[nFree].setCoefficient(powers, double(limit) * Error(i+1));
      powers[i] = 0;
      powers[i+1] = 1;
      fixedPoint[nFree].setCoefficient(powers, - double(limit) * Error(i));
      powers[i+1] = 0;
    }

    // Path length term in delta^limit.
    fixedPoint[nFree].setCoefficient(delta, Q[nFree].getCoefficient(delta));
  }

  // Fill in path length terms and identity for delta(p)/p.
  fixedPointMap[nFree+1] = FTps<double,N>::makeVariable(nFree + 1);
  fixedPointMap[nFree] =
    fixedPoint.inverse()[nFree].substitute(tempMap).substitute(fixedPoint);

  // Combine fixed point and dispersion.
  fixedPoint += fixPoint;
  //std::cerr << "Leaving FStaticFP<N>::FStaticFP(...)" << std::endl;
}


template <int N>
FStaticFP<N>::~FStaticFP()
{}


template <int N>
const FVps<double,N> &FStaticFP<N>::getFixedPoint() const
{
  return fixedPoint;
}


template <int N>
const FVps<double,N> &FStaticFP<N>::getFixedPointMap() const
{
  return fixedPointMap;
}

#endif // MAD_FStaticFP_HH
