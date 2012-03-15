// ------------------------------------------------------------------------
// $RCSfile: Diagnostic.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Diagnostic
//   Defines the abstract interface for a beam diagnostic.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Diagnostic
// ------------------------------------------------------------------------

Diagnostic::Diagnostic():
  Component()
{}


Diagnostic::Diagnostic(const Diagnostic &rhs):
  Component(rhs)
{}


Diagnostic::Diagnostic(const string &name):
  Component(name)
{}


Diagnostic::~Diagnostic()
{}


void Diagnostic::accept(BeamlineVisitor &visitor) const
{
  visitor.visitDiagnostic(*this);
}
