// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

//////////////////////////////////////////////////////////////////////
#ifndef GEN_VEKTOR_H
#define GEN_VEKTOR_H

//////////////////////////////////////////////////////////////////////

template<unsigned Dim, unsigned IDim>
class UnitComponentVektor
{
private:
  // No actual data.
public:
  // Allow access rather like a regular Vektor.
  bool operator[](const unsigned int &d) const { return (d==IDim); }
};

//////////////////////////////////////////////////////////////////////

template<class T, unsigned Dim, unsigned IDim>
class ComponentVektor
{
private:
  T Value;

public:
  // Construct with or without a value.
  ComponentVektor() {}
  ComponentVektor(T value) : Value(value) {}

  // Allow access rather like a regular Vektor.
  T& operator[](const unsigned int &d) { return (d==IDim) ? Value : 0; }
  T  operator[](const unsigned int &d) const { return (d==IDim) ? Value : 0; }
};

//////////////////////////////////////////////////////////////////////

// Dot product of a UnitComponentVektor with a regular Vektor.
// The UnitComponentVektor just selects a value, no floating point ops.

// In the first two it takes in a Vektor& and returns a T&.
// That way this could be put on the left hand side like
// dot(vec,zhat) = 1.0;

template<class T, unsigned Dim, unsigned IDim>
inline
T& dot_ref(Vektor<T,Dim>& v, const UnitComponentVektor<Dim,IDim>&)
{
  return v[IDim];
}

template<class T, unsigned Dim, unsigned IDim>
inline
T& dot_ref(const UnitComponentVektor<Dim,IDim>& , Vektor<T,Dim>& v)
{
  return v[IDim];
}

// If the Vektor is const though, return by value so you
// can't assign to it.

template<class T, unsigned Dim, unsigned IDim>
inline
T dot(const Vektor<T,Dim>& v, const UnitComponentVektor<Dim,IDim>&)
{
  return v[IDim];
}

template<class T, unsigned Dim, unsigned IDim>
inline
T dot(const UnitComponentVektor<Dim,IDim>& , const Vektor<T,Dim>& v)
{
  return v[IDim];
}

//////////////////////////////////////////////////////////////////////

//
// These dot a regular Vektor with a ComponentVektor.
// They just do one multipply.
// The return by value of course.
//

template<class T, unsigned Dim, unsigned IDim>
inline
T dot(const Vektor<T,Dim>& v, const ComponentVektor<T,Dim,IDim>& c)
{
  return v[IDim]*c[IDim];
}

template<class T, unsigned Dim, unsigned IDim>
inline
T dot(const ComponentVektor<T,Dim,IDim>& c, const Vektor<T,Dim>& v)
{
  return v[IDim]*c[IDim];
}


#endif
/***************************************************************************
 * $RCSfile: GenVektor.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:24 $
 * IPPL_VERSION_ID: $Id: GenVektor.h,v 1.1.1.1 2003/01/23 07:40:24 adelmann Exp $ 
 ***************************************************************************/
