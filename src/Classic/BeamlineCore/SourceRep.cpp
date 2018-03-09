#include "BeamlineCore/SourceRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(SourceRep::*get)() const;
        void (SourceRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &SourceRep::getElementLength,
            &SourceRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


// Class SourceRep
// ------------------------------------------------------------------------

SourceRep::SourceRep():
    Source(),
    geometry()
{}


SourceRep::SourceRep(const SourceRep &right):
    Source(right),
    geometry(right.geometry)
{}


SourceRep::SourceRep(const std::string &name):
    Source(name), geometry()
{}


SourceRep::~SourceRep()
{}


ElementBase *SourceRep::clone() const {
    return new SourceRep(*this);
}


Channel *SourceRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<SourceRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

NullField &SourceRep::getField() {
    return field;
}

const NullField &SourceRep::getField() const {
    return field;
}

StraightGeometry &SourceRep::getGeometry() {
    return geometry;
}

const StraightGeometry &SourceRep::getGeometry() const {
    return geometry;
}


ElementImage *SourceRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}