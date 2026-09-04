//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "BeamlineCore/MultipoleTRep.h"

#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char* name;
        double (MultipoleTRep::*get)() const;
        void (MultipoleTRep::*set)(double);
    };

    const Entry entries[] = {
            {"BY", &MultipoleTRep::getB, &MultipoleTRep::setB}, {nullptr, nullptr, nullptr}};
}  // namespace

MultipoleTRep::MultipoleTRep() : MultipoleT(), geometry_m(Geometry::makeStraight(0.0)) {}

MultipoleTRep::MultipoleTRep(const MultipoleTRep& right)
    : MultipoleT(right), geometry_m(right.geometry_m) {}

MultipoleTRep::MultipoleTRep(const std::string& name)
    : MultipoleT(name), geometry_m(Geometry::makeStraight(0.0)) {}

MultipoleTRep::~MultipoleTRep() = default;

ElementBase* MultipoleTRep::clone() const {
    return new MultipoleTRep(*this);
}

Channel* MultipoleTRep::getChannel(const std::string& aKey, bool create) {
    for (const Entry* entry = entries; entry->name != nullptr; ++entry) {
        if (aKey == entry->name) {
            return new IndirectChannel<MultipoleTRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

Geometry& MultipoleTRep::getGeometry() {
    return geometry_m;
}

const Geometry& MultipoleTRep::getGeometry() const {
    return geometry_m;
}
