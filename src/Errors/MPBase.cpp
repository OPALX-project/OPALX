// ------------------------------------------------------------------------
// $RCSfile: MPBase.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MPBase
//   Abstract mixin class for all OPAL objects which handle multipole errors.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/MPBase.h"
#include "AbstractObjects/OpalData.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"


// Class MPBase
// ------------------------------------------------------------------------


MPBase::MPBase()
{}


MPBase::~MPBase()
{}


double MPBase::absFactor(int comp)
{
  double fact = OPAL.getP0() / Physics::c;
  while (comp > 1) fact /= double(comp--);
  return fact;
}


double MPBase::relFactor(int comp, int order, double radius)
{
  if (comp < order) {
    double fact = 1.0;
    while (comp++ < order) fact *= radius;
    return fact;
  } else {
    double fact = 1.0;
    while (comp-- > order) fact /= radius;
    return fact;
  }
}
