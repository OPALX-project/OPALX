// ------------------------------------------------------------------------
// $RCSfile: DriftRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DriftRep
//   Defines a concrete drift space representation.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/DriftRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(DriftRep::*get)() const;
        void (DriftRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &DriftRep::getElementLength,
            &DriftRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


// Class DriftRep
// ------------------------------------------------------------------------

DriftRep::DriftRep():
    Drift(),
    geometry(0.0)
{}


DriftRep::DriftRep(const DriftRep &right):
    Drift(right),
    geometry(right.geometry)
{}


DriftRep::DriftRep(const string &name):
    Drift(name),
    geometry()
{}


DriftRep::~DriftRep()
{}


ElementBase *DriftRep::clone() const {
    return new DriftRep(*this);
}


Channel *DriftRep::getChannel(const string &aKey) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<DriftRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey);
}


NullField &DriftRep::getField() {
    return field;
}

const NullField &DriftRep::getField() const {
    return field;
}


StraightGeometry &DriftRep::getGeometry() {
    return geometry;
}

const StraightGeometry &DriftRep::getGeometry() const {
    return geometry;
}


ElementImage *DriftRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}


