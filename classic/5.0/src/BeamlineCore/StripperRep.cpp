// ------------------------------------------------------------------------
// $RCSfile: StripperRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: StripperRep
//   Defines a representation for a septa.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2011/07/08 10:51:59 $
// $Author: Jianjun Yang $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/StripperRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(StripperRep::*get)() const;
        void (StripperRep::*set)(double);
    };

    static const Entry entries[] = {
        {
            "L",
            &StripperRep::getElementLength,
            &StripperRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


// Class StripperRep
// ------------------------------------------------------------------------

StripperRep::StripperRep():
    Stripper(), field(), geometry(), active(true)
{}


StripperRep::StripperRep(const StripperRep &right):
    Stripper(right), field(), geometry(right.geometry), active(true)
{}


StripperRep::StripperRep(const string &name):
    Stripper(name), field(), geometry(), active(true)
{}


StripperRep::~StripperRep()
{}


ElementBase *StripperRep::clone() const {
    return new StripperRep(*this);
}


Channel *StripperRep::getChannel(const string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<StripperRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

NullField &StripperRep::getField() {
    return field;
}

const NullField &StripperRep::getField() const {
    return field;
}

StraightGeometry &StripperRep::getGeometry() {
    return geometry;
}

const StraightGeometry &StripperRep::getGeometry() const {
    return geometry;
}


ElementImage *StripperRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}

void StripperRep::setActive(bool flag) {
    active = flag;
}

