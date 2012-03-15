// ------------------------------------------------------------------------
// $RCSfile: Separator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Separator
//   Defines the abstract interface for an  separator.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Separator.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Separator
// ------------------------------------------------------------------------

Separator::Separator():
  Component()
{}


Separator::Separator(const Separator &right):
  Component(right)
{}


Separator::Separator(const string &name):
  Component(name)
{}


Separator::~Separator()
{}


void Separator::accept(BeamlineVisitor &visitor) const
{
  visitor.visitSeparator(*this);
}
