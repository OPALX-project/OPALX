// ------------------------------------------------------------------------
// $RCSfile: Lambertson.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Lambertson
//   Defines the abstract interface for a Lambertson septum magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Lambertson
// ------------------------------------------------------------------------

Lambertson::Lambertson():
    Component()
{ }


Lambertson::Lambertson(const Lambertson &rhs):
    Component(rhs)
{ }


Lambertson::Lambertson(const string &name):
    Component(name)
{ }


Lambertson::~Lambertson()
{ }


void Lambertson::accept(BeamlineVisitor &visitor) const {
    visitor.visitLambertson(*this);
}

bool Lambertson::apply(const int &i, const double &t, double E[], double B[]) {
    return false;
}

bool Lambertson::apply(const int &i, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

bool Lambertson::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

void Lambertson::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    RefPartBunch_m = bunch;
}

void Lambertson::finalise()
{ }

bool Lambertson::bends() const {
    return false;
}

void Lambertson::getDimensions(double &zBegin, double &zEnd) const {

}

const string &Lambertson::getType() const {
    static const string type("Lambertson");
    return type;
}

