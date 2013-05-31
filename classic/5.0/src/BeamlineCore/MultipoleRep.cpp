// ------------------------------------------------------------------------
// $RCSfile: MultipoleRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MultipoleRep
//   Defines a concrete representation for a general multipole.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2001/08/24 19:31:20 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/MultipoleRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndexedChannel.h"
#include "Channels/IndirectChannel.h"
#include "ComponentWrappers/MultipoleWrapper.h"
#include <cctype>

// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(MultipoleRep::*get)() const;
        void (MultipoleRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &MultipoleRep::getElementLength,
            &MultipoleRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


// Class MultipoleRep
// ------------------------------------------------------------------------

MultipoleRep::MultipoleRep():
    Multipole(),
    geometry(),
    field()
{}


MultipoleRep::MultipoleRep(const MultipoleRep &multipole):
    Multipole(multipole),
    geometry(multipole.geometry),
    field(multipole.field)
{}


MultipoleRep::MultipoleRep(const string &name):
    Multipole(name),
    geometry(),
    field()
{}


MultipoleRep::~MultipoleRep()
{}


ElementBase *MultipoleRep::clone() const {
    return new MultipoleRep(*this);
}


Channel *MultipoleRep::getChannel(const string &aKey, bool create) {
    if(aKey[0] == 'A'  ||  aKey[0] == 'B') {
        int n = 0;

        for(string::size_type k = 1; k < aKey.length(); k++) {
            if(isdigit(aKey[k])) {
                n = 10 * n + aKey[k] - '0';
            } else {
                return 0;
            }
        }

        if(aKey[0] == 'B') {
            return new IndexedChannel<MultipoleRep>
                   (*this, &MultipoleRep::getNormalComponent,
                    &MultipoleRep::setNormalComponent, n);
        } else {
            return new IndexedChannel<MultipoleRep>
                   (*this, &MultipoleRep::getSkewComponent,
                    &MultipoleRep::setSkewComponent, n);
        }
    } else {
        for(const Entry *entry = entries; entry->name != 0; ++entry) {
            if(aKey == entry->name) {
                return new IndirectChannel<MultipoleRep>
                       (*this, entry->get, entry->set);
            }
        }

        return ElementBase::getChannel(aKey, create);
    }
}


StraightGeometry &MultipoleRep::getGeometry() {
    return geometry;
}

const StraightGeometry &MultipoleRep::getGeometry() const {
    return geometry;
}


ElementImage *MultipoleRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    for(int n = 1; n <= field.order(); n++) {
        char buffer[20];
        char *p = buffer;
        int k = n;

        while(k != 0) {
            *p++ = k % 10 + '0';
            k /= 10;
        }

        string name(" ");
        while(p > buffer) name += *--p;

        double b = getNormalComponent(n);
        if(b != 0.0) {
            name[0] = 'b';
            image->setAttribute(name, b);
        }

        double a = getSkewComponent(n);
        if(a != 0.0) {
            name[0] = 'a';
            image->setAttribute(name, a);
        }
    }

    return image;
}


BMultipoleField &MultipoleRep::getField() {
    return field;
}


const BMultipoleField &MultipoleRep::getField() const {
    return field;
}


void MultipoleRep::setField(const BMultipoleField &aField) {
    field = aField;
}


ElementBase *MultipoleRep::makeFieldWrapper() {
    ElementBase *wrap = new MultipoleWrapper(this);
    wrap->setName(getName());
    return wrap;
}
