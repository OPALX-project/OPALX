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

BoxRep::BoxRep() : Box(), geometry() {}

BoxRep::BoxRep(const BoxRep& right) : Box(right), geometry(right.geometry) {}

BoxRep::BoxRep(const std::string& name) : Box(name), geometry() {}

BoxRep::~BoxRep() {}

ElementBase* BoxRep::clone() const { return new BoxRep(*this); }

Geometry& BoxRep::getGeometry() { return geometry; }

const Geometry& BoxRep::getGeometry() const { return geometry; }
