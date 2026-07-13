//
// Class MonitorRep
//   Representation for an orbit position monitor.
//   The base class observes both planes.
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
#include "BeamlineCore/MonitorRep.h"
#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char* name;
        double (MonitorRep::*get)() const;
        void (MonitorRep::*set)(double);
    };

    const Entry entries[] = {{0, 0, 0}};
}  // namespace

MonitorRep::MonitorRep() : Monitor(), geometry(), active(true) {}

MonitorRep::MonitorRep(const MonitorRep& right)
    : Monitor(right), geometry(right.geometry), active(true) {}

MonitorRep::MonitorRep(const std::string& name) : Monitor(name), geometry(), active(true) {}

MonitorRep::~MonitorRep() {}

ElementBase* MonitorRep::clone() const { return new MonitorRep(*this); }

Channel* MonitorRep::getChannel(const std::string& aKey, bool create) {
    for (const Entry* entry = entries; entry->name != 0; ++entry) {
        if (aKey == entry->name) {
            return new IndirectChannel<MonitorRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

Geometry& MonitorRep::getGeometry() { return geometry; }

const Geometry& MonitorRep::getGeometry() const { return geometry; }

Monitor::Plane MonitorRep::getPlane() const { return active ? XY : OFF; }

void MonitorRep::setActive(bool flag) { active = flag; }
