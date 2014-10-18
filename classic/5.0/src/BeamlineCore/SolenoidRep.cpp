// ------------------------------------------------------------------------
// $RCSfile: SolenoidRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SolenoidRep
//   Defines a concrete representation for a solenoid.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/SolenoidRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(SolenoidRep::*get)() const;
        void (SolenoidRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &SolenoidRep::getElementLength,
            &SolenoidRep::setElementLength
        },
        {
            "BZ",
            &SolenoidRep::getBz,
            &SolenoidRep::setBz
        },
        { 0, 0, 0 }
    };
}


// Class SolenoidRep
// ------------------------------------------------------------------------

SolenoidRep::SolenoidRep():
    Solenoid(),
    geometry(),
    field()
{}


SolenoidRep::SolenoidRep(const SolenoidRep &right):
    Solenoid(right),
    geometry(right.geometry),
    field(right.field)
{}


SolenoidRep::SolenoidRep(const std::string &name):
    Solenoid(name), geometry(), field()
{}


SolenoidRep::~SolenoidRep()
{}


ElementBase *SolenoidRep::clone() const {
    return new SolenoidRep(*this);
}


Channel *SolenoidRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<SolenoidRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


ConstBzField &SolenoidRep::getField() {
    return field;
}

const ConstBzField &SolenoidRep::getField() const {
    return field;
}


StraightGeometry &SolenoidRep::getGeometry() {
    return geometry;
}

const StraightGeometry &SolenoidRep::getGeometry() const {
    return geometry;
}


ElementImage *SolenoidRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}


double SolenoidRep::getBz() const {
    return field.getBz();
}


void SolenoidRep::setBz(double Bz) {
    field.setBz(Bz);
}
