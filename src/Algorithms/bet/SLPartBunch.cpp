#ifdef HAVE_ENVELOPE_SOLVER
// ------------------------------------------------------------------------
// $RCSfile: SLPartBunch.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class SLPartBunch
//   Interface to a particle bunch.
//   Can be used to avoid use of a template in user code.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/bet/SLPartBunch.h"
#include <iostream>
#include <cfloat>
#include <fstream>
#include "Physics/Physics.h"
#include "AbstractObjects/OpalData.h"

using Physics::c;
using Physics::pi;

// Class SLPartBunch
// ------------------------------------------------------------------------

SLPartBunch::SLPartBunch(const PartData *ref):
  reference(ref),
  PartBunch(ref)
{
}


SLPartBunch::SLPartBunch(const SLPartBunch &rhs):
  reference(rhs.reference),
  PartBunch(rhs.reference)
{}


SLPartBunch::SLPartBunch(const std::vector<Particle> &rhs, const PartData *ref):
  reference(ref),
  PartBunch(ref)
{}


SLPartBunch::~SLPartBunch()
{

}

Inform &SLPartBunch::print(Inform &os)
{

  return os;
}

#endif
