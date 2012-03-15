// ------------------------------------------------------------------------
// $RCSfile: Corrector.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Corrector
//   Defines the abstract interface for a orbit corrector.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Corrector
// ------------------------------------------------------------------------

Corrector::Corrector():
  Component()
{}


Corrector::Corrector(const Corrector &right):
  Component(right)
{}


Corrector::Corrector(const string &name):
  Component(name)
{}


Corrector::~Corrector()
{}


void Corrector::accept(BeamlineVisitor &visitor) const
{
  visitor.visitCorrector(*this);
}
