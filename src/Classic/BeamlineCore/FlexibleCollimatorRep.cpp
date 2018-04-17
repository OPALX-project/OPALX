// ------------------------------------------------------------------------
// $RCSfile: FlexibleCollimatorRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: FlexibleCollimatorRep
//   Defines a concrete collimator representation.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/FlexibleCollimatorRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(FlexibleCollimatorRep::*get)() const;
        void (FlexibleCollimatorRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &FlexibleCollimatorRep::getElementLength,
            &FlexibleCollimatorRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


// Class FlexibleCollimatorRep
// ------------------------------------------------------------------------

FlexibleCollimatorRep::FlexibleCollimatorRep():
    FlexibleCollimator(),
    geometry(0.0)
{}


FlexibleCollimatorRep::FlexibleCollimatorRep(const FlexibleCollimatorRep &right):
    FlexibleCollimator(right),
    geometry(right.geometry)
{}


FlexibleCollimatorRep::FlexibleCollimatorRep(const std::string &name):
    FlexibleCollimator(name),
    geometry()
{}


FlexibleCollimatorRep::~FlexibleCollimatorRep()
{}


ElementBase *FlexibleCollimatorRep::clone() const {
    return new FlexibleCollimatorRep(*this);
}


Channel *FlexibleCollimatorRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<FlexibleCollimatorRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


NullField &FlexibleCollimatorRep::getField() {
    return field;
}

const NullField &FlexibleCollimatorRep::getField() const {
    return field;
}


StraightGeometry &FlexibleCollimatorRep::getGeometry() {
    return geometry;
}

const StraightGeometry &FlexibleCollimatorRep::getGeometry() const {
    return geometry;
}


ElementImage *FlexibleCollimatorRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}