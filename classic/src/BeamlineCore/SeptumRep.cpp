// ------------------------------------------------------------------------
// $RCSfile: SeptumRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SeptumRep
//   Defines a representation for a septa.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2009/09/21 10:26:06 $
// $Author: bi $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/SeptumRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(SeptumRep::*get)() const;
        void (SeptumRep::*set)(double);
    };

    static const Entry entries[] = {
        {
            "L",
            &SeptumRep::getElementLength,
            &SeptumRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


// Class SeptumRep
// ------------------------------------------------------------------------

SeptumRep::SeptumRep():
    Septum(), field(), geometry(), active(true)
{}


SeptumRep::SeptumRep(const SeptumRep &right):
    Septum(right), field(), geometry(right.geometry), active(true)
{}


SeptumRep::SeptumRep(const std::string &name):
    Septum(name), field(), geometry(), active(true)
{}


SeptumRep::~SeptumRep()
{}


ElementBase *SeptumRep::clone() const {
    return new SeptumRep(*this);
}


Channel *SeptumRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<SeptumRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

NullField &SeptumRep::getField() {
    return field;
}

const NullField &SeptumRep::getField() const {
    return field;
}

StraightGeometry &SeptumRep::getGeometry() {
    return geometry;
}

const StraightGeometry &SeptumRep::getGeometry() const {
    return geometry;
}


ElementImage *SeptumRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}

void SeptumRep::setActive(bool flag) {
    active = flag;
}
