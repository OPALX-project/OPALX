#ifndef CLASSIC_VpsInvMap_HH
#define CLASSIC_VpsInvMap_HH

// ------------------------------------------------------------------------
// $RCSfile: VpsInvMap.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: VpsInvMap
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:18:06 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Algebra/LUMatrix.h"
#include "Algebra/Matrix.h"
#include "Algebra/Vector.h"
#include "Algebra/VpsInvMap.h"
#include "Algebra/Tps.h"
#include "Algebra/VpsMap.h"
#include "Utilities/LogicalError.h"
#include <iosfwd>


// Class VpsInvMap
// ------------------------------------------------------------------------
/// Invertible power series map
//  An invertible truncated power series map.
//  The number of variables and the dimension are equal.

template <class T>
class VpsInvMap: public VpsMap<T> {

public:

  /// Constructor.
  //  Construct identity in [b]nDim[/b] variables.
  VpsInvMap(int nDim);

  /// Convert.
  //  Throw SizeError if [b]rhs[/b] is not R**n --> R**n.
  VpsInvMap(const VpsMap<T> &rhs);

  /// Convert.
  //  Throw SizeError if [b]rhs[/b] is not R**n --> R**n.
  VpsInvMap(const Vps<T> &rhs);

  /// Constructor.
  //  Construct map of dimension [b]nDim[/b] and fill from [b]rhs[/b].
  VpsInvMap(int nDim, const Tps<T> rhs[]);

  /// Convert.
  //  The constant terms are zero.
  //  The linear terms are taken from [b]M[/b].
  //  Throw SizeError if [b]M[/b] is not square.
  VpsInvMap(const Matrix<T> &M);

  /// Convert.
  //  The constant terms are taken from [b]V[/b].
  //  The linear terms are the identity map.
  VpsInvMap(const Vector<T> &);

  VpsInvMap();
  VpsInvMap(const VpsInvMap<T> &rhs); 
  ~VpsInvMap();
  VpsInvMap<T> &operator=(const VpsInvMap<T> &y);

  /// Check consistency.
  void check() const;

  /// Inverse.
  VpsInvMap<T> inverse() const;

  /// Set to identity.
  static VpsInvMap<T> identity(int nVar);
};


// ------------------------------------------------------------------------
///{ Global Operators on VpsInvMap<T>

/// Multiply.
template <class T>
  VpsInvMap<T> operator*(const VpsInvMap<T> &x, const Tps<T> &y);

/// Multiply.
template <class T>
  VpsInvMap<T> operator*(const Tps<T> &x, const VpsInvMap<T> &y);

/// Multiply.
template <class T>
  VpsInvMap<T> operator*(const VpsInvMap<T> &x, const T &y);

/// Multiply.
template <class T>
  VpsInvMap<T> operator*(const T &x, const VpsInvMap<T> &y);

/// Divide.
//  Throw DivideError, if constant part of [b]y[/b] is zero.
template <class T>
  VpsInvMap<T> operator/(const VpsInvMap<T> &x, const Tps<T> &y);

/// Divide.
//  Throw DivideError, if [b]y[/b] is zero.
template <class T>
  VpsInvMap<T> operator/(const VpsInvMap<T> &x, const T &y);

/// Add.
template <class T>
  VpsInvMap<T> operator+(const VpsInvMap<T> &x, const VpsInvMap<T> &y);

/// Subtract.
template <class T>
  VpsInvMap<T> operator-(const VpsInvMap<T> &x, const VpsInvMap<T> &y);

/// Add.
template <class T>
  VpsInvMap<T> operator+(const VpsInvMap<T> &x, const Vector<T> &y);

/// Subtract.
template <class T> 
 VpsInvMap<T> operator-(const VpsInvMap<T> &x, const Vector<T> &y);

/// Add.
template <class T>
  VpsInvMap<T> operator+(const Vector<T> &x, const VpsInvMap<T> &y);

/// Subtract.
template <class T>
  VpsInvMap<T> operator-(const Vector<T> &x, const VpsInvMap<T> &y);

/// Extract from stream.
template <class T>
  std::istream &operator>>(std::istream &, VpsInvMap<T> &x);

/// Insert to stream.
template <class T>
  std::ostream &operator<<(std::ostream &, const VpsInvMap<T> &x);

/// Multiply.
template <class T> VpsInvMap<T>
  operator*(const Matrix<T> &x, const VpsInvMap<T> &y);


// Implementation of global functions for class Vps<T>.
// ------------------------------------------------------------------------

template <class T> inline
VpsInvMap<T> operator*(const VpsInvMap<T> &x, const Tps<T> &y)
{ VpsInvMap<T> z(x); return z *= y; }

template <class T> inline
VpsInvMap<T> operator*(const Tps<T> &x, const VpsInvMap<T> &y)
{ VpsInvMap<T> z(x); return z *= y; }

template <class T> inline
VpsInvMap<T> operator*(const VpsInvMap<T> &x, const T &y)
{ VpsInvMap<T> z(x); return z *= y; }

template <class T> inline
VpsInvMap<T> operator*(const T &x, const VpsInvMap<T> &y)
{ VpsInvMap<T> z(x); return z *= y; }

template <class T> inline
VpsInvMap<T> operator/(const VpsInvMap<T> &x, const Tps<T> &y)
{ VpsInvMap<T> z(x); return z /= y; }

template <class T> inline
VpsInvMap<T> operator/(const VpsInvMap<T> &x, const T &y)
{ VpsInvMap<T> z(x); return z /= y; }


template <class T> inline
VpsInvMap<T> operator+(const VpsInvMap<T> &x, const VpsInvMap<T> &y)
{ VpsInvMap<T> z(x); return z += y; }
  
template <class T> inline
VpsInvMap<T> operator-(const VpsInvMap<T> &x, const VpsInvMap<T> &y)
{ VpsInvMap<T> z(x); return z -= y; }

template <class T> inline
VpsInvMap<T> operator+(const VpsInvMap<T> &x, const Vector<T> &y)
{ VpsInvMap<T> z(x); return z += y; }
  
template <class T> inline
VpsInvMap<T> operator-(const VpsInvMap<T> &x, const Vector<T> &y)
{ VpsInvMap<T> z(x); return z -= y; }

template <class T> inline
VpsInvMap<T> operator+(const Vector<T> &x, const VpsInvMap<T> &y)
{ VpsInvMap<T> z(y); return z += x; }

template <class T> inline
VpsInvMap<T> operator-(const Vector<T> &x, const VpsInvMap<T> &y)
{ VpsInvMap<T> z(- y); return z += x; }

template <class T> inline
VpsInvMap<T> operator*(const Matrix<T> &x, const VpsInvMap<T> &y)
{ return y.substituteInto(x); }

template <class T> inline
std::istream &operator>>(std::istream &is, VpsInvMap<T> &x)
{ return x.get(is); }

template <class T> inline
std::ostream &operator<<(std::ostream &os, const VpsInvMap<T> &x)
{ return x.put(os); }


// Class VpsInvMap, implementation.
// ------------------------------------------------------------------------

template <class T>
VpsInvMap<T>::VpsInvMap():
  VpsMap<T>()
{}


template <class T>
VpsInvMap<T>::VpsInvMap(int nDim):
  VpsMap<T>(nDim, nDim)
{
  for(int i = 0; i < nDim; i++) this->data[i] = Tps<T>::makeVariable(nDim, i);
}


template <class T>
VpsInvMap<T>::VpsInvMap(const VpsInvMap<T> &rhs):
  VpsMap<T>(rhs)
{}


template <class T>
VpsInvMap<T>::VpsInvMap(const VpsMap<T> &rhs):
  VpsMap<T>(rhs)
{}


template <class T>
VpsInvMap<T>::VpsInvMap(const Vps<T> &rhs):
  VpsMap<T>(rhs)
{}


template <class T>
VpsInvMap<T>::VpsInvMap(int nDim, const Tps<T> rhs[]):
  VpsMap<T>(nDim, rhs)
{
  check();
}


template <class T>
VpsInvMap<T>::VpsInvMap(const Matrix<T> &x):
  VpsMap<T>(x)
{
  if (x.nrows() != x.ncols()) {
    throw LogicalError("VpsInvMap::VpsInvMap()", "Matrix is not n by n.");
  }
}


template <class T>
VpsInvMap<T>::VpsInvMap(const Vector<T> &x):
  VpsMap<T>(x)
{
}


template <class T>
VpsInvMap<T>::~VpsInvMap()
{
}


template <class T>
VpsInvMap<T> &VpsInvMap<T>::operator=(const VpsInvMap &y)
{
  VpsMap<T>::operator=(y);
  return *this;
}


template <class T>
void VpsInvMap<T>::check() const
{
  Vps<T>::check();
  if (this->getDimension() != this->variables) {
    throw LogicalError("VpsInvMap::VpsInvMap()", "Map is not n->n.");
  }
}

template <class T>
VpsInvMap<T> VpsInvMap<T>::inverse() const
{
  // Separate out linear terms and change sign of all other terms:
  //    M := A(1), the linear part of map A,
  //    map1 := - (A - A(1)).
  check();
  Vector<T> vec = this->constantTerm();
  LUMatrix<T> lu(this->linearTerms());
  Matrix<T> mat(lu.inverse());
  VpsInvMap<T> map1 = - filter(2, Tps<T>::EXACT);
  
  // Initialize map B = M**(-1) (with truncation at maxOrder).
  VpsInvMap<T> B(mat);
  B -= mat * vec;
  
  // Compute inverse order by order.
  for (int maxOrder = 2; maxOrder <= this->getTruncOrder(); maxOrder++) {
    VpsInvMap<T> map2 = map1.substitute(B, maxOrder);
    B = mat * (map2 + identity(this->variables) - vec);
  }
  
  return VpsInvMap<T>(B);
}


template <class T>
VpsInvMap<T> VpsInvMap<T>::identity(int nDim)
{
  return VpsInvMap<T>(nDim);
}

#endif // CLASSIC_VpsInvMap_HH
