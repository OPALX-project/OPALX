// ------------------------------------------------------------------------
// $RCSfile: RBend.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RBend
//   Defines the abstract interface for a rectangular bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/RBend.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class RBend
// ------------------------------------------------------------------------

RBend::RBend():
  Component()
{}


RBend::RBend(const RBend &right):
  Component(right)
{}


RBend::RBend(const string &name):
  Component(name)
{}


RBend::~RBend()
{}


void RBend::accept(BeamlineVisitor &visitor) const
{
  visitor.visitRBend(*this);
}


double RBend::getNormalComponent(int n) const
{
  return getField().getNormalComponent(n);
}


double RBend::getSkewComponent(int n) const
{
  return getField().getSkewComponent(n);
}


void RBend::setNormalComponent(int n, double v)
{
  getField().setNormalComponent(n, v);
}


void RBend::setSkewComponent(int n, double v)
{
  getField().setSkewComponent(n, v);
}
