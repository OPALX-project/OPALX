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

MonitorRep::MonitorRep() : Monitor(), geometry() {}

MonitorRep::MonitorRep(const MonitorRep& right) : Monitor(right), geometry(right.geometry) {}

MonitorRep::MonitorRep(const std::string& name) : Monitor(name), geometry() {}

MonitorRep::~MonitorRep() {}

ElementBase* MonitorRep::clone() const { return new MonitorRep(*this); }

Geometry& MonitorRep::getGeometry() { return geometry; }

const Geometry& MonitorRep::getGeometry() const { return geometry; }
