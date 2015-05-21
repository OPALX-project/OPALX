// ------------------------------------------------------------------------
// $RCSfile: RFQuadrupole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RFQuadrupole
//   Defines the abstract interface for a RF quadrupole.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/BeamlineVisitor.h"

extern Inform *gmsg;

// Class RFQuadrupole
// ------------------------------------------------------------------------

RFQuadrupole::RFQuadrupole():
    Component()
{}


RFQuadrupole::RFQuadrupole(const RFQuadrupole &rhs):
    Component(rhs)
{}


RFQuadrupole::RFQuadrupole(const std::string &name):
    Component(name)
{}


RFQuadrupole::~RFQuadrupole()
{}


void RFQuadrupole::accept(BeamlineVisitor &visitor) const {
    visitor.visitRFQuadrupole(*this);
}

bool RFQuadrupole::apply(const size_t &i, const double &t, double E[], double B[]) {
    return false;
}

bool RFQuadrupole::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

bool RFQuadrupole::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

void RFQuadrupole::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    RefPartBunch_m = bunch;
}

void RFQuadrupole::finalise()
{}

bool RFQuadrupole::bends() const {
    return false;
}


void RFQuadrupole::getDimensions(double &zBegin, double &zEnd) const {

}


ElementBase::ElementType RFQuadrupole::getType() const {
    return RFQUADRUPOLE;
}

