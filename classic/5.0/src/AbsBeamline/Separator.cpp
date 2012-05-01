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


void Separator::accept(BeamlineVisitor &visitor) const {
    visitor.visitSeparator(*this);
}

bool Separator::apply(const size_t &i, const double &t, double E[], double B[]) {
    return false;
}

bool Separator::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

bool Separator::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

void Separator::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    RefPartBunch_m = bunch;
}

void Separator::finalise()
{}

bool Separator::bends() const {
    return false;
}

void Separator::getDimensions(double &zBegin, double &zEnd) const {

}


const string &Separator::getType() const {
    static const string type("Separator");
    return type;
}

