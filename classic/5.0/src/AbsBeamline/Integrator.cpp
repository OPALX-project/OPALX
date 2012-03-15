// ------------------------------------------------------------------------
// $RCSfile: Integrator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Integrator
//   Integrator is a pure abstract class.  It forms the base class for all
//   special propagators through an element or a beamline.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Integrator.h"


// Class Integrator
// ------------------------------------------------------------------------

Integrator::Integrator(ElementBase *elem):
  ElementBase(elem->getName()), itsElement(elem)
{}


Integrator::Integrator(const Integrator &rhs):
  ElementBase(rhs.getName()), itsElement(rhs.itsElement)
{}


Integrator::~Integrator()
{}


void Integrator::makeSharable()
{
  shareFlag = true;
  itsElement->makeSharable();
}
