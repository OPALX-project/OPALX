// ------------------------------------------------------------------------
// $RCSfile: RFCavityRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RFCavityRep
//   Defines a representation for a RF cavity.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/RFCavityRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

bool RFCavityRep::ignoreCavities = false;

namespace {
    struct Entry {
        const char *name;
        double(RFCavityRep::*get)() const;
        void (RFCavityRep::*set)(double);
    };

    static const Entry entries[] = {
        {
            "L",
            &RFCavityRep::getElementLength,
            &RFCavityRep::setElementLength
        },
        {
            "AMPLITUDE",
            &RFCavityRep::getAmplitude,
            &RFCavityRep::setAmplitude
        },
        {
            "FREQUENCY",
            &RFCavityRep::getFrequency,
            &RFCavityRep::setFrequency
        },
        {
            "PHASE",
            &RFCavityRep::getPhase,
            &RFCavityRep::setPhase
        },
        { 0, 0, 0 }
    };
}


// Class RFCavityRep
// ------------------------------------------------------------------------

RFCavityRep::RFCavityRep():
    RFCavity()
{}


RFCavityRep::RFCavityRep(const RFCavityRep &right):
    RFCavity(right)
{}


RFCavityRep::RFCavityRep(const string &name):
    RFCavity(name)
{}


RFCavityRep::~RFCavityRep()
{}


ElementBase *RFCavityRep::clone() const {
    return new RFCavityRep(*this);
}


Channel *RFCavityRep::getChannel(const string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<RFCavityRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


AcceleratingField &RFCavityRep::getField() {
    return field;
}

const AcceleratingField &RFCavityRep::getField() const {
    return field;
}


StraightGeometry &RFCavityRep::getGeometry() {
    return geometry;
}

const StraightGeometry &RFCavityRep::getGeometry() const {
    return geometry;
}


ElementImage *RFCavityRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}


double RFCavityRep::getAmplitude() const {
    return ignoreCavities ? 0.0 : field.getEz();
}


double RFCavityRep::getFrequency() const {
    return field.getFrequency();
}


double RFCavityRep::getPhase() const {
    return field.getPhase();
}


void RFCavityRep::setAmplitude(double amplitude) {
    field.setEz(amplitude);
}


void RFCavityRep::setFrequency(double frequency) {
    field.setFrequency(frequency);
}


void RFCavityRep::setPhase(double phase) {
    field.setPhase(phase);
}

void RFCavityRep::setIgnore(bool ignore) {
    ignoreCavities = ignore;
}
