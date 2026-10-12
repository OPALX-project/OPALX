//
// Class BeamBeamRep
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
#include "BeamlineCore/BeamBeamRep.h"
#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char* name;
        double (BeamBeamRep::*get)() const;
        void (BeamBeamRep::*set)(double);
    };

    const Entry entries[] = {{0, 0, 0}};
}  // namespace

BeamBeamRep::BeamBeamRep() : BeamBeam(), geometry() {}

BeamBeamRep::BeamBeamRep(const BeamBeamRep& right) : BeamBeam(right), geometry(right.geometry) {}

BeamBeamRep::BeamBeamRep(const std::string& name) : BeamBeam(name), geometry() {}

BeamBeamRep::~BeamBeamRep() {}

ElementBase* BeamBeamRep::clone() const { return new BeamBeamRep(*this); }

Channel* BeamBeamRep::getChannel(const std::string& aKey, bool create) {
    for (const Entry* entry = entries; entry->name != 0; ++entry) {
        if (aKey == entry->name) {
            return new IndirectChannel<BeamBeamRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

Geometry& BeamBeamRep::getGeometry() { return geometry; }

const Geometry& BeamBeamRep::getGeometry() const { return geometry; }
