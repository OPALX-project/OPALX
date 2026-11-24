//
// Class TravelingWaveRep
//   Representation for a traveling wave.
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
#include "BeamlineCore/TravelingWaveRep.h"

TravelingWaveRep::TravelingWaveRep() : TravelingWave() {}

TravelingWaveRep::TravelingWaveRep(const TravelingWaveRep& right) : TravelingWave(right) {}

TravelingWaveRep::TravelingWaveRep(const std::string& name) : TravelingWave(name) {}

TravelingWaveRep::~TravelingWaveRep() {}

ElementBase* TravelingWaveRep::clone() const { return new TravelingWaveRep(*this); }
