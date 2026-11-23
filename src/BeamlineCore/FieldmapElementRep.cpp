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

FieldmapElementRep::FieldmapElementRep() : FieldmapElement(), geometry() {}

FieldmapElementRep::FieldmapElementRep(const FieldmapElementRep& right)
    : FieldmapElement(right), geometry(right.geometry) {}

FieldmapElementRep::FieldmapElementRep(const std::string& name)
    : FieldmapElement(name), geometry() {}

FieldmapElementRep::~FieldmapElementRep() {}

ElementBase* FieldmapElementRep::clone() const { return new FieldmapElementRep(*this); }

Geometry& FieldmapElementRep::getGeometry() { return geometry; }

const Geometry& FieldmapElementRep::getGeometry() const { return geometry; }
