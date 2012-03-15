// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by PSI. 
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit www.amas.web.psi for more details
 *
 ***************************************************************************/

// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// include files
#include "SubParticle/SubParticleAttrib.h"
#include "Particle/ParticleAttrib.h"
#include "Utility/IpplInfo.h"
#include "Utility/PAssert.h"
#include "Profile/Profiler.h"


//////////////////////////////////////////////////////////////////////
// assignment of another SubParticleAttrib
template<class PA, class T, unsigned Dim>
SubParticleAttrib<PA,T,Dim>&
SubParticleAttrib<PA,T,Dim>::operator=(const SubParticleAttrib<PA,T,Dim> &t) {
  TAU_TYPE_STRING(taustr, CT(*this) + " (" + CT(t) + " )" );
  TAU_PROFILE("SubParticleAttrib::operator=()", taustr,
	      TAU_SUBFIELD | TAU_ASSIGN);

  // just do a regular attribute assignment
  //A = t.getAttrib();
  assign(*this, t);
  return *this;
}


//////////////////////////////////////////////////////////////////////
// assignment of a scalar
template<class PA, class T, unsigned Dim>
SubParticleAttrib<PA,T,Dim>&
SubParticleAttrib<PA,T,Dim>::operator=(T t) {
  TAU_TYPE_STRING(taustr, CT(*this) + " (" + CT(t) + " )" );
  TAU_PROFILE("SubParticleAttrib::operator=()", taustr,
	      TAU_SUBFIELD | TAU_ASSIGN);

  // make a PETE_Scalar so that we can do a regular assignment
  //A = t;
  PETE_Scalar<T> scalar(t);
  assign(*this, scalar);
  return *this;
}




//////////////////////////////////////////////////////////////////////
// Return the beginning and end iterators for this class.
template<class PA, class T, unsigned Dim>
typename SubParticleAttrib<PA,T,Dim>::iterator
SubParticleAttrib<PA,T,Dim>::begin() const {
  TAU_TYPE_STRING(taustr, CT(*this) + "::iterator ()" );
  TAU_PROFILE("SubParticleAttrib::begin()", taustr, TAU_SUBFIELD );

  PA &p = const_cast<PA &>(A);
  return iterator(p, MyDomain.begin_iv(), 0, MyDomain);
}

template<class PA, class T, unsigned Dim>
typename SubParticleAttrib<PA,T,Dim>::iterator
SubParticleAttrib<PA,T,Dim>::end() const {
  TAU_TYPE_STRING(taustr, CT(*this) + "::iterator ()" );
  TAU_PROFILE("SubParticleAttrib::end()", taustr, TAU_SUBFIELD );

  PA &p = const_cast<PA &>(A);
  return iterator(p, MyDomain.end_iv(), MyDomain.size(), MyDomain);
}

//////////////////////////////////////////////////////////////////////
// Make sure the LHS ParticleAttrib has the proper length.  It should
// have the same length as the number of LOCAL sparse index points.
// If it does not, we adjust the length.  Return size of result.
template<class PA, class T, unsigned Dim>
int SubParticleAttrib<PA,T,Dim>::adjustSize() {
  TAU_TYPE_STRING(taustr, CT(*this) + " ()" );
  TAU_PROFILE("SubParticleAttrib::adjustSize()", taustr, TAU_SUBFIELD);

  // get length of sparse index list local pieces
  int points = MyDomain.size();
  int currpoints = A.size();

  // adjust length of attrib
  if (points < currpoints)
    A.destroy((currpoints - points), points);
  else if (points > currpoints)
    A.create(points - currpoints);

  return points;
}


//////////////////////////////////////////////////////////////////////
// write out the contents of this SubParticleAttrib to the given ostream.
// This not only prints out the values of the attribute, but also prints
// out the local sindex points
template<class PA, class T, unsigned Dim>
void SubParticleAttrib<PA,T,Dim>::write(std::ostream &o) const {
  TAU_TYPE_STRING(taustr, CT(*this) + " (" + CT(o) + " )" );
  TAU_PROFILE("SubParticleAttrib::write()", taustr, TAU_SUBFIELD );

  // make sure the sizes match
  PInsist(A.size() == MyDomain.size(),
	  "SubParticleAttrib::write must have an attrib with enough elements");

  // loop over all the local lfields, printing out particles as we go
  int i = 0;
  typename Index_t::const_iterator_iv lfi = MyDomain.begin_iv();
  for ( ; lfi != MyDomain.end_iv(); ++lfi) {
    // for each lfield, print out particles
    typename Index_t::const_iterator_indx lsi = (*lfi)->begin();
    for ( ; lsi != (*lfi)->end(); ++lsi)
      o << *lsi + MyDomain.getOffset() << " ==> " << A[i++] << endl;
  }

}


/***************************************************************************
 * $RCSfile: SubParticleAttrib.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: SubParticleAttrib.cpp,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/
