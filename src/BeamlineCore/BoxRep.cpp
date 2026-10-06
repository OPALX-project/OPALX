//
// Class BoxRep
//   Representation for a rectangular absorbing block.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "BeamlineCore/BoxRep.h"
#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char* name;
        double (BoxRep::*get)() const;
        void (BoxRep::*set)(double);
    };

    const Entry entries[] = {{0, 0, 0}};
}  // namespace

BoxRep::BoxRep() : Box(), geometry() {}

BoxRep::BoxRep(const BoxRep& right) : Box(right), geometry(right.geometry) {}

BoxRep::BoxRep(const std::string& name) : Box(name), geometry() {}

BoxRep::~BoxRep() {}

ElementBase* BoxRep::clone() const { return new BoxRep(*this); }

Channel* BoxRep::getChannel(const std::string& aKey, bool create) {
    for (const Entry* entry = entries; entry->name != 0; ++entry) {
        if (aKey == entry->name) {
            return new IndirectChannel<BoxRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

Geometry& BoxRep::getGeometry() { return geometry; }

const Geometry& BoxRep::getGeometry() const { return geometry; }
