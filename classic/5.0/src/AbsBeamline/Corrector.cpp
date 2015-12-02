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
#include "Algorithms/PartBunch.h"
#include "Physics/Physics.h"
#include "Utilities/GeneralClassicException.h"

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


Corrector::Corrector(const std::string &name):
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
  B[0] = kickField_m(0);
  B[1] = kickField_m(1);

  return false;
}

bool Corrector::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
  B += kickField_m;
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

void Corrector::goOnline(const double &kineticEnergy) {
    if (kineticEnergy < 0.0) {
        throw GeneralClassicException("Corrector::goOnline", "given kinetic energy is negative");
    }


    const double pathLength = getGeometry().getElementLength();
    const double minLength = Physics::c * RefPartBunch_m->getdT();
    if (pathLength < minLength) throw GeneralClassicException("Corrector::goOnline", "length of corrector, L= " + std::to_string(pathLength) + ", shorter than distance covered during one time step, dS= " + std::to_string(minLength));

    const double momentum = sqrt(std::pow(kineticEnergy * 1e6, 2.0) + 2.0 * kineticEnergy * 1e6 * RefPartBunch_m->getM());
    const double magnitude = momentum / (Physics::c * pathLength);
    kickField_m = magnitude * RefPartBunch_m->getQ() * Vector_t(kickY_m, -kickX_m, 0.0);

    online_m = true;
}

bool Corrector::bends() const {
    return false;
}

void Corrector::getDimensions(double &zBegin, double &zEnd) const
{
  zBegin = startField_m;
  zEnd = startField_m + getElementLength();
}

ElementBase::ElementType Corrector::getType() const {
    return CORRECTOR;
}