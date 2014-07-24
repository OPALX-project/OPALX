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
  Component(),
  startField_m(0.0),
  endField_m(0.0),
  kickX_m(0.0),
  kickY_m(0.0)
{ }


Corrector::Corrector(const Corrector &right):
  Component(right),
  startField_m(right.startField_m),
  endField_m(right.endField_m),
  kickX_m(right.kickX_m),
  kickY_m(right.kickY_m)
{ }


Corrector::Corrector(const string &name):
  Component(name),
  startField_m(0.0),
  endField_m(0.0),
  kickX_m(0.0),
  kickY_m(0.0)
{ }


Corrector::~Corrector()
{ }


void Corrector::accept(BeamlineVisitor &visitor) const {
    visitor.visitCorrector(*this);
}

bool Corrector::apply(const size_t &i, const double &t, double E[], double B[]) {
  Inform m("Corrector::apply 1" );
  const double xk = GetKickX();
  const double yk = GetKickY();
  B[0] = xk;
  B[1] = yk;
  return false;
}

bool Corrector::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {

  const double xk = GetKickX();
  const double yk = GetKickY();
  B = Vector_t(xk,yk,0.0);
  return false;
}

bool Corrector::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
  return false;
}

void Corrector::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
  endField_m = endField = startField + getElementLength();
  RefPartBunch_m = bunch;
  startField_m = startField;
}

void Corrector::finalise()
{ }

bool Corrector::bends() const {
    return false;
}

void Corrector::getDimensions(double &zBegin, double &zEnd) const
{ 
  zBegin = startField_m;
  zEnd = startField_m + getElementLength();
}

const string &Corrector::getType() const {
    static const string type("Corrector");
    return type;
}

void Corrector::SetKickX(double k) {kickX_m = k; }

void Corrector::SetKickY(double k) {kickY_m = k; }

double Corrector::GetKickX() const {return kickX_m; }

double Corrector::GetKickY() const {return kickY_m; }
