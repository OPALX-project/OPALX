//
// Class IpRep
//   Representation for a drift space.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "BeamlineCore/IpRep.h"
#include "Channels/IndirectChannel.h"


namespace {
    struct Entry {
        const char *name;
        double(IpRep::*get)() const;
        void (IpRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &IpRep::getElementLength,
            &IpRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


IpRep::IpRep():
    Ip(),
    geometry(0.0)
{}


IpRep::IpRep(const IpRep &right):
    Ip(right),
    geometry(right.geometry)
{}


IpRep::IpRep(const std::string &name):
    Ip(name),
    geometry()
{}


IpRep::~IpRep()
{}


ElementBase *IpRep::clone() const {
    return new IpRep(*this);
}


Channel *IpRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<IpRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


NullField &IpRep::getField() {
    return field;
}

const NullField &IpRep::getField() const {
    return field;
}


StraightGeometry &IpRep::getGeometry() {
    return geometry;
}

const StraightGeometry &IpRep::getGeometry() const {
    return geometry;
}
