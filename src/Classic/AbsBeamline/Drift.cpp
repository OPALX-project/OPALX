// ------------------------------------------------------------------------
// $RCSfile: Drift.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Drift
//   Defines the abstract interface for a drift space.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Drift.h"
#include "AbsBeamline/BeamlineVisitor.h"

extern Inform *gmsg;

// Class Drift
// ------------------------------------------------------------------------

Drift::Drift():
    Component()
{ }


Drift::Drift(const Drift &right):
    Component(right)
{ }


Drift::Drift(const std::string &name):
    Component(name) {

}

Drift::~Drift()
{ }


void Drift::accept(BeamlineVisitor &visitor) const {
    visitor.visitDrift(*this);
}

bool Drift::apply(const size_t &i, const double &t, double E[], double B[]) {
    return false;
}

bool Drift::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

bool Drift::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

void Drift::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    endField = startField + getElementLength();
    RefPartBunch_m = bunch;
    startField_m = startField;
}

void Drift::finalise() {
}

bool Drift::bends() const {
    return false;
}

void Drift::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = startField_m + getElementLength();
}

ElementBase::ElementType Drift::getType() const {
    return DRIFT;
}
