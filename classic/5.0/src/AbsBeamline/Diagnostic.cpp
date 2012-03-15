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
{ }


Diagnostic::Diagnostic(const Diagnostic &rhs):
    Component(rhs)
{ }


Diagnostic::Diagnostic(const string &name):
    Component(name)
{ }


Diagnostic::~Diagnostic()
{ }


void Diagnostic::accept(BeamlineVisitor &visitor) const
{
    visitor.visitDiagnostic(*this);
}

bool Diagnostic::apply(const int &i, const double &t, double E[], double B[])
{
    return false;
}

bool Diagnostic::apply(const int &i, const double &t, Vector_t &E, Vector_t &B)
{
    return false;
}
  
bool Diagnostic::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B)
{
    return false;
}

void Diagnostic::initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor)
{
    RefPartBunch_m = bunch;
}

void Diagnostic::finalise()
{ }

void Diagnostic::rescaleFieldMap(const double &scaleFactor)
{ }

bool Diagnostic::bends() const
{
    return false;
}

const string& Diagnostic::getType() const
{
    static const string type("Diagnostic");
    return type;
}

void Diagnostic::getDimensions(double &zBegin, double &zEnd) const
{ }

