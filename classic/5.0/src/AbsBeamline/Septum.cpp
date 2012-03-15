// ------------------------------------------------------------------------
// $RCSfile: Septum.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Septum
//   Defines the abstract interface for a septum magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Septum.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Septum
// ------------------------------------------------------------------------

Septum::Septum():
  Component()
{}


Septum::Septum(const Septum &rhs):
  Component(rhs)
{}


Septum::Septum(const string &name):
  Component(name)
{}


Septum::~Septum()
{}


void Septum::accept(BeamlineVisitor &visitor) const
{
  visitor.visitSeptum(*this);
}
