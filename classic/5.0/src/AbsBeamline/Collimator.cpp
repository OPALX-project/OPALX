// ------------------------------------------------------------------------
// $RCSfile: Collimator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Collimator
//   Defines the abstract interface for a beam Collimator.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Collimator.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Collimator
// ------------------------------------------------------------------------

Collimator::Collimator():
  Component()
{}


Collimator::Collimator(const Collimator &rhs):
  Component(rhs)
{}


Collimator::Collimator(const string &name):
  Component(name)
{}


Collimator::~Collimator()
{}


void Collimator::accept(BeamlineVisitor &visitor) const
{
  visitor.visitCollimator(*this);
}
