// ------------------------------------------------------------------------
// $RCSfile: Marker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Marker
//   Defines the abstract interface for a marker element.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Marker.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Marker
// ------------------------------------------------------------------------


Marker::Marker():
  Component()
{}


Marker::Marker(const Marker &right):
  Component(right)
{}


Marker::Marker(const string &name):
  Component(name)
{}


Marker::~Marker()
{}


void Marker::accept(BeamlineVisitor &visitor) const
{
  visitor.visitMarker(*this);
}

bool Marker::apply(const int &i, const double &t, double E[], double B[])
{
  return false;
}

bool Marker::apply(const int &i, const double &t, Vector_t &E, Vector_t &B)
{
  return false;
}
  
bool Marker::apply(const Vector_t &R, const double &t, Vector_t &E, Vector_t &B)
{
  return false;
}

void Marker::initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor)
{
  RefPartBunch_m = bunch;
}

void Marker::finalise()
{}

void Marker::rescaleFieldMap(const double &scaleFactor)
{}

bool Marker::bends() const
{
  return false;
}

