//
// Class MultipoleRep
//   Representation for a general multipole.
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
#include "BeamlineCore/MultipoleRep.h"

MultipoleRep::MultipoleRep() : Multipole() {}

MultipoleRep::MultipoleRep(const MultipoleRep& multipole) : Multipole(multipole) {}

MultipoleRep::MultipoleRep(const std::string& name) : Multipole(name) {}

MultipoleRep::~MultipoleRep() {}

ElementBase* MultipoleRep::clone() const { return new MultipoleRep(*this); }
