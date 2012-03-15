// ------------------------------------------------------------------------
// $RCSfile: Patch.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Patch
//   Defines the abstract interface for a geometry patch.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Patch.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Patch
// ------------------------------------------------------------------------

Patch::Patch():
    Component()
{}


Patch::Patch(const Patch &rhs):
    Component(rhs)
{}


Patch::Patch(const string &name):
    Component(name)
{}


Patch::~Patch()
{}


void Patch::accept(BeamlineVisitor &visitor) const {
    visitor.visitPatch(*this);
}

bool Patch::apply(const int &i, const double &t, double E[], double B[]) {
    return false;
}

bool Patch::apply(const int &i, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

bool Patch::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

void Patch::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    RefPartBunch_m = bunch;
}

void Patch::finalise()
{}

bool Patch::bends() const {
    return false;
}

void Patch::getDimensions(double &zBegin, double &zEnd) const {

}


const string &Patch::getType() const {
    static const string type("Patch");
    return type;
}

