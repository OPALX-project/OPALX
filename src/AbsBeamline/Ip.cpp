// ------------------------------------------------------------------------
// $RCSfile: Ip.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Ip
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

#include "AbsBeamline/Ip.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "PartBunch/PartBunch.h"

extern Inform *gmsg;

// Class Ip
// ------------------------------------------------------------------------

Ip::Ip():
    Ip("")
{ }


Ip::Ip(const Ip &right):
    Component(right),
    nSlices_m(right.nSlices_m)
{ }


Ip::Ip(const std::string &name):
    Component(name),
    nSlices_m(1)
{ }


Ip::~Ip()
{ }


void Ip::accept(BeamlineVisitor &visitor) const {
    visitor.visitIp(*this);
}

void Ip::initialise(PartBunch_t *bunch, double &startField, double &endField) {
    endField = startField + getElementLength();
    RefPartBunch_m = bunch;
    startField_m = startField;
}


//set the number of slices for map tracking
void Ip::setNSlices(const std::size_t& nSlices) { 
    nSlices_m = nSlices;
}

//get the number of slices for map tracking
std::size_t Ip::getNSlices() const {
    return nSlices_m;
}

void Ip::finalise() {
}

bool Ip::bends() const {
    return false;
}

void Ip::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = startField_m + getElementLength();
}

ElementType Ip::getType() const {
    return ElementType::DRIFT;
}
