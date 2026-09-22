//
// Class FieldmapElementRep
//   Representation for a field-map-driven element.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#include "BeamlineCore/FieldmapElementRep.h"
#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char* name;
        double (FieldmapElementRep::*get)() const;
        void (FieldmapElementRep::*set)(double);
    };

    const Entry entries[] = {
            {"SCALE", &FieldmapElementRep::getScale, &FieldmapElementRep::setScale},
            {"ESCALE", &FieldmapElementRep::getEScale, &FieldmapElementRep::setEScale},
            {0, 0, 0}};
}  // namespace

FieldmapElementRep::FieldmapElementRep() : FieldmapElement(), geometry() {}

FieldmapElementRep::FieldmapElementRep(const FieldmapElementRep& right)
    : FieldmapElement(right), geometry(right.geometry) {}

FieldmapElementRep::FieldmapElementRep(const std::string& name)
    : FieldmapElement(name), geometry() {}

FieldmapElementRep::~FieldmapElementRep() {}

ElementBase* FieldmapElementRep::clone() const { return new FieldmapElementRep(*this); }

Channel* FieldmapElementRep::getChannel(const std::string& aKey, bool create) {
    for (const Entry* entry = entries; entry->name != 0; ++entry) {
        if (aKey == entry->name) {
            return new IndirectChannel<FieldmapElementRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

Geometry& FieldmapElementRep::getGeometry() { return geometry; }

const Geometry& FieldmapElementRep::getGeometry() const { return geometry; }
