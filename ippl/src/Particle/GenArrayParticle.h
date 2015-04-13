// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef GEN_ARRAY_PARTICLE_H
#define GEN_ARRAY_PARTICLE_H

/*
 * GenArrayParticle - Specialized Particle container for N data elements.
 *
 * GenArrayParticle is a special version of a ParticleBase class, specialized
 * to have D-dimensional coordinates of type T1, and N extra data items of
 * type T2.
 */

// include files
#include "Particle/ParticleBase.h"


// GenArrayParticle class definition
template<class PLayout, class T, unsigned N>
class GenArrayParticle : public ParticleBase<PLayout> {

public:
  // attributes for this class: an array with N attributes
  ParticleAttrib<T> data[N];

  // constructor: user-provided Layout object must be supplied
  GenArrayParticle(PLayout* L) : ParticleBase<PLayout>(L) {
    for (int i = 0; i < N; i++)
      addAttribute(data[i]);
  }

private:
  // disable default constructor
  GenArrayParticle();

};

#endif // GEN_ARRAY_PARTICLE_H

/***************************************************************************
 * $RCSfile: GenArrayParticle.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:28 $
 * IPPL_VERSION_ID: $Id: GenArrayParticle.h,v 1.1.1.1 2003/01/23 07:40:28 adelmann Exp $ 
 ***************************************************************************/
