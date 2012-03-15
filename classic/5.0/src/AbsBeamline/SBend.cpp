// ------------------------------------------------------------------------
// $RCSfile: SBend.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: SBend
//   Defines the abstract interface for a sector bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/SBend.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class SBend
// ------------------------------------------------------------------------

SBend::SBend():
  Component()
{}


SBend::SBend(const SBend &right):
  Component(right)
{}


SBend::SBend(const string &name):
  Component(name)
{}


SBend::~SBend()
{}


void SBend::accept(BeamlineVisitor &visitor) const
{
  visitor.visitSBend(*this);
}


double SBend::getNormalComponent(int n) const
{
  return getField().getNormalComponent(n);
}


double SBend::getSkewComponent(int n) const
{
  return getField().getSkewComponent(n);
}


void SBend::setNormalComponent(int n, double v)
{
  getField().setNormalComponent(n, v);
}


void SBend::setSkewComponent(int n, double v)
{
  getField().setSkewComponent(n, v);
}
