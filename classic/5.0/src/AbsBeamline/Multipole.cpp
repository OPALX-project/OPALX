// ------------------------------------------------------------------------
// $RCSfile: Multipole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Multipole
//   Defines the abstract interface for a Multipole magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Multipole
// ------------------------------------------------------------------------

Multipole::Multipole():
  Component()
{}


Multipole::Multipole(const Multipole &right):
  Component(right)
{}


Multipole::Multipole(const string &name):
  Component(name)
{}


Multipole::~Multipole()
{}


void Multipole::accept(BeamlineVisitor &visitor) const
{
  visitor.visitMultipole(*this);
}


double Multipole::getNormalComponent(int n) const
{
  return getField().getNormalComponent(n);
}


double Multipole::getSkewComponent(int n) const
{
  return getField().getSkewComponent(n);
}


void Multipole::setNormalComponent(int n, double v)
{
  getField().setNormalComponent(n, v);
}


void Multipole::setSkewComponent(int n, double v)
{
  getField().setSkewComponent(n, v);
}
