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

bool Septum::apply(const int &i, const double &t, double E[], double B[])
{
  return false;
}

bool Septum::apply(const int &i, const double &t, Vector_t &E, Vector_t &B)
{
  return false;
}
  
bool Septum::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B)
{
  return false;
}

void Septum::initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor)
{
  RefPartBunch_m = bunch;
}

void Septum::finalise()
{}

void Septum::rescaleFieldMap(const double &scaleFactor)
{}

bool Septum::bends() const
{
  return false;
}


void Septum::getDimensions(double &zBegin, double &zEnd) const
{

}


const string& Septum::getType() const
{
    static const string type("Septum");
    return type;
}

