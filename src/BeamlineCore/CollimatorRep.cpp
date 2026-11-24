//
// Class CollimatorRep
//   Representation for a transverse-aperture collimator.
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
#include "BeamlineCore/CollimatorRep.h"

CollimatorRep::CollimatorRep() : Collimator() {}

CollimatorRep::CollimatorRep(const CollimatorRep& right) : Collimator(right) {}

CollimatorRep::CollimatorRep(const std::string& name) : Collimator(name) {}

CollimatorRep::~CollimatorRep() {}
