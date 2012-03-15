#ifndef CLASSIC_LieMap_HH
#define CLASSIC_LieMap_HH

// ------------------------------------------------------------------------
// $RCSfile: LieMap.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: LieMap
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:18:05 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Algebra/Matrix.h"
#include "Algebra/Vector.h"
#include "Algebra/Tps.h"
#include "Algebra/TpsData.h"
#include "Algebra/VpsInvMap.h"
#include "Utilities/ConvergenceError.h"
#include "Utilities/LogicalError.h"
#include <iosfwd>


// Template class LieMap<T>
// ------------------------------------------------------------------------
/// Lie algebraic map.
//  LieMap<T> is a truncated Taylor series map with coefficients of type
//  [b]T[/b].  It implements operations available only for symplectic maps.

template <class T>
class LieMap: public VpsInvMap<T> {

public:

  /// Constructor.
  //  Construct identity with [b]nDim[/b] variables.
  LieMap(int nDim);

  /// Convert.
  //  Construct a map from the matrix [b]M[/b], which should be symplectic.
  LieMap(const Matrix<T> &M);

  /// Convert.
  //  Construct an identity map with a displacement vector [b]V[/b].
  LieMap(const Vector<T> &V);

  /// Conversion.
  //  Convert possibly non-symplectic map.
  LieMap(const Vps<T> &rhs);
  
  LieMap();
  LieMap(const LieMap<T> &rhs);
  ~LieMap();
  LieMap<T> &operator=(const LieMap<T> &);

  /// Check consistency.
  void check() const;

  /// Lie series.
  //  Build the exponential series exp(:H:)Z from the Tps<T> H,
  //  acting on the identity map [b]Z[/b].
  static LieMap<T> ExpMap(const Tps<T> &H);

  /// Lie series.
  //  Build the exponential series exp(:H:)M from the Tps<T> H,
  //  acting on the map M.
  static LieMap<T> ExpMap(const Tps<T> &H, const LieMap<T> &M);

private:

  // The iteration limit for Lie series convergence.
  static const int MAXITER = 100;
};


// Operations on Vps<T>.
// ------------------------------------------------------------------------

/// Poisson bracket.
template <class T>
Tps<T> PoissonBracket(const Tps<T> &x, const Tps<T> &y);

/// Poisson bracket.
template <class T>
Vps<T> PoissonBracket(const Tps<T> &x, const Vps<T> &y);

/// Extract LieMap<T> from stream.
template <class T>
std::istream &operator>>(std::istream &, LieMap<T> &x);

/// Insert LieMap<T> to stream.
template <class T>
std::ostream &operator<<(std::ostream &, const LieMap<T> &x);


// Implementation of global functions for class LieMap<T>.
// ------------------------------------------------------------------------

template <class T>
LieMap<T>::LieMap(): VpsInvMap<T>()
{}


template <class T>
LieMap<T>::LieMap(int nDim): VpsInvMap<T>(nDim)
{
  check();
}


template <class T>
LieMap<T>::LieMap(const LieMap<T> &rhs): VpsInvMap<T>(rhs)
{}


template <class T>
LieMap<T>::LieMap(const Matrix<T> &x): VpsInvMap<T>(x)
{}


template <class T>
LieMap<T>::LieMap(const Vector<T> &x): VpsInvMap<T>(x)
{}


template <class T>
LieMap<T>::LieMap(const Vps<T> &rhs):
  VpsInvMap<T>(rhs)
{
  check();
}


template <class T> inline
LieMap<T>::~LieMap()
{}


template <class T> inline
LieMap<T> &LieMap<T>::operator=(const LieMap<T> &rhs)
{
  VpsInvMap<T>::operator=(rhs);
  return *this;
}


template <class T>
Tps<T> PoissonBracket(const Tps<T> &x, const Tps<T> &y)
{
  Tps<T> z;
  int nDim = x.getDimension();
    
  for (int q = 0; q < nDim; q +=2) {
    int p = q + 1;
    z += x.derivative(q) * y.derivative(p) - x.derivative(p) * y.derivative(q);
  }

  return z;
}


template <class T>
Vps<T> PoissonBracket(const Tps<T> &x, const Vps<T> &y)
{
  Vps<T> z;
  int nDim = x.getDimension();
    
  for (int q = 0; q < nDim; q +=2) {
    int p = q + 1;
    z += x.derivative(q) * y.derivative(p) - x.derivative(p) * y.derivative(q);
  }

  return z;
}


template <class T>
LieMap<T> LieMap<T>::ExpMap(const Tps<T> &H)
{
  return ExpMap(H, VpsInvMap<T>::identity(H.getVariables()));
}


template <class T>
LieMap<T> LieMap<T>::ExpMap(const Tps<T> &H, const LieMap<T> &M)
{
  int nDim = H.getVariables();
  if (nDim % 2 != 0) {
    throw LogicalError("LieMap::ExpMap()", "Tps dimension is zero or odd.");
  }
  if (nDim == 0) return M;
  if (nDim != M.getVariables()) {
    throw LogicalError
      ("LieMap::ExpMap()",
       "Lie map and Hamiltonian have inconsistent dimensions.");
  }

  LieMap<T> temp(nDim);
  Vps<T> dH(nDim, nDim);

  for (int i = 0; i < nDim; i += 2) {
    dH[i]   = - H.derivative(i+1);
    dH[i+1] =   H.derivative(i);
  }

  for (int var = 0; var < nDim; var++) {
    Tps<T> old = 0.0;
    Tps<T> u = temp[var] = M[var];

    for (int count = 1; temp[var] != old; count++) {
      if (count >= MAXITER) {
	throw ConvergenceError("LieMap::ExpMap()",
			       "No convergence in ExpMap()");
      }

      Tps<T> w = 0.0;
      for (int v = 0; v < nDim; v++) w += dH[v] * u.derivative(v);
      u = w / T(count);
      old = temp[var];
      temp[var] += u;
    }
  }

  return temp;
}


template <class T>
std::istream &operator>>(std::istream &is, LieMap<T> &x)
{ return x.get(is); }


template <class T>
std::ostream &operator<<(std::ostream &os, const LieMap<T> &x)
{ return x.put(os); }


template <class T>
void LieMap<T>::check() const
{
  VpsInvMap<T>::check();

  if (this->getDimension() % 2 != 0) {
    throw LogicalError("LieMap::check()", "LieMap representation corrupted.");
  }
}

#endif // CLASSIC_LieMap_HH
