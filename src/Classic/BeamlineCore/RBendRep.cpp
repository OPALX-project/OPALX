// ------------------------------------------------------------------------
// $RCSfile: RBendRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RBendRep
//   Defines a concrete representation for rectangular (straight) bend.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/RBendRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndexedChannel.h"
#include "Channels/IndirectChannel.h"
#include "ComponentWrappers/RBendWrapper.h"
#include <cctype>

// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(RBendRep::*get)() const;
        void (RBendRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &RBendRep::getElementLength,
            &RBendRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


// Class RBendRep
// ------------------------------------------------------------------------

RBendRep::RBendRep():
    RBend(),
    geometry(0.0, 0.0),
    field() {
    rEntry = rExit = hEntry = hExit = 0.0;
}


RBendRep::RBendRep(const RBendRep &rhs):
    RBend(rhs),
    geometry(rhs.geometry),
    field(rhs.field) {
    rEntry = rhs.rEntry;
    rExit  = rhs.hExit;
    hEntry = rhs.rEntry;
    hExit  = rhs.hExit;
}


RBendRep::RBendRep(const std::string &name):
    RBend(name),
    geometry(0.0, 0.0),
    field() {
    rEntry = rExit = hEntry = hExit = 0.0;
}


RBendRep::~RBendRep()
{}


ElementBase *RBendRep::clone() const {
    return new RBendRep(*this);
}


Channel *RBendRep::getChannel(const std::string &aKey, bool create) {
    return ElementBase::getChannel(aKey, create);
}


BMultipoleField &RBendRep::getField() {
    return field;
}

const BMultipoleField &RBendRep::getField() const {
    return field;
}


RBendGeometry &RBendRep::getGeometry() {
    return geometry;
}

const RBendGeometry &RBendRep::getGeometry() const {
    return geometry;
}


ElementImage *RBendRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}


double RBendRep::getB() const {
    return field.getNormalComponent(1);
}

void RBendRep::setB(double B) {
    field.setNormalComponent(1, B);
}


double RBendRep::getEntryFaceRotation() const {
    return rEntry;
}


double RBendRep::getExitFaceRotation() const {
    return rExit;
}

void RBendRep::setEntryFaceRotation(double e1) {
    rEntry = e1;
}

void RBendRep::setExitFaceRotation(double e2) {
    rExit = e2;
}

double RBendRep::getEntryFaceCurvature() const {
    return hEntry;
}

double RBendRep::getExitFaceCurvature() const {
    return hExit;
}

void RBendRep::setEntryFaceCurvature(double h1) {
    hEntry = h1;
}

void RBendRep::setExitFaceCurvature(double h2) {
    hExit = h2;
}


double RBendRep::getSlices() const {
    return slices;
}

double RBendRep::getStepsize() const {
    return stepsize;
}

void RBendRep::setSlices(double sl) {
    slices = sl;
}

void RBendRep::setStepsize(double ds) {
    stepsize = ds;
}


void RBendRep::setField(const BMultipoleField &f) {
    field = f;
}


ElementBase *RBendRep::makeFieldWrapper() {
    ElementBase *wrap = new RBendWrapper(this);
    wrap->setName(getName());
    return wrap;
}