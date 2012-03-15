// ------------------------------------------------------------------------
// $RCSfile: Corrector.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Corrector
//   Defines the abstract interface for a orbit corrector.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Corrector
// ------------------------------------------------------------------------

Corrector::Corrector():
  Component()
{}


Corrector::Corrector(const Corrector &right):
  Component(right)
{}


Corrector::Corrector(const string &name):
  Component(name)
{}


Corrector::~Corrector()
{}


void Corrector::accept(BeamlineVisitor &visitor) const
{
  visitor.visitCorrector(*this);
}

bool Corrector::apply(const int &i, const double &t, double E[], double B[])
{
  return false;
}

bool Corrector::apply(const int &i, const double &t, Vector_t &E, Vector_t &B)
{
  return false;
}
  
bool Corrector::apply(const Vector_t &R, const double &t, Vector_t &E, Vector_t &B)
{
  return false;
}

void Corrector::initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor)
{
  RefPartBunch_m = bunch;
}

void Corrector::finalise()
{}

void Corrector::rescaleFieldMap(const double &scaleFactor)
{}

bool Corrector::bends() const
{
  return false;
}

void Corrector::getDimensions(double &zBegin, double &zEnd) const
{

}

