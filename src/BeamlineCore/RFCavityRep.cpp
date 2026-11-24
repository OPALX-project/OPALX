//
// Class RFCavityRep
//   Representation for a RF cavity.
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
#include "BeamlineCore/RFCavityRep.h"

RFCavityRep::RFCavityRep() : RFCavity() {}

RFCavityRep::RFCavityRep(const RFCavityRep& right) : RFCavity(right) {}

RFCavityRep::RFCavityRep(const std::string& name) : RFCavity(name) {}

RFCavityRep::~RFCavityRep() {}

ElementBase* RFCavityRep::clone() const { return new RFCavityRep(*this); }
