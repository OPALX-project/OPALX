// ------------------------------------------------------------------------
// $RCSfile: Drift.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Drift
//   Defines the abstract interface for a drift space.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Drift.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Drift
// ------------------------------------------------------------------------

Drift::Drift():
  Component()
{}


Drift::Drift(const Drift &right):
  Component(right)
{}


Drift::Drift(const string &name):
  Component(name)
{

}

Drift::~Drift()
{}


void Drift::accept(BeamlineVisitor &visitor) const
{
  visitor.visitDrift(*this);
}
