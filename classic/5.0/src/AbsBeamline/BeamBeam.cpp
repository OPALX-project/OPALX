// ------------------------------------------------------------------------
// $RCSfile: BeamBeam.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamBeam
//   Defines the abstract interface for a beam-beam interaction.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/BeamBeam.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class BeamBeam
// ------------------------------------------------------------------------

BeamBeam::BeamBeam():
    Component()
{}


BeamBeam::BeamBeam(const BeamBeam &right):
    Component(right)
{}


BeamBeam::BeamBeam(const string &name):
    Component(name)
{}


BeamBeam::~BeamBeam()
{}


void BeamBeam::accept(BeamlineVisitor &visitor) const {
    visitor.visitBeamBeam(*this);
}

bool BeamBeam::apply(const size_t &i, const double &t, double E[], double B[]) {
    return false;
}

bool BeamBeam::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

bool BeamBeam::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

void BeamBeam::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    RefPartBunch_m = bunch;
}

void BeamBeam::finalise()
{}

bool BeamBeam::bends() const {
    return false;
}

void BeamBeam::getDimensions(double &zBegin, double &zEnd) const {

}


const string &BeamBeam::getType() const {
    static const string type("BeamBeam");
    return type;
}

