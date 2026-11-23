//
// Class SolenoidRep
//   Representation for a solenoid magnet.
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
#include "BeamlineCore/SolenoidRep.h"

SolenoidRep::SolenoidRep() : Solenoid(), geometry() {}

SolenoidRep::SolenoidRep(const SolenoidRep& right)
    : Solenoid(right), geometry(right.geometry), Bz_m(right.Bz_m) {}

SolenoidRep::SolenoidRep(const std::string& name) : Solenoid(name), geometry() {}

SolenoidRep::~SolenoidRep() {}

ElementBase* SolenoidRep::clone() const { return new SolenoidRep(*this); }

Geometry& SolenoidRep::getGeometry() { return geometry; }

const Geometry& SolenoidRep::getGeometry() const { return geometry; }

double SolenoidRep::getBz() const { return Bz_m; }

void SolenoidRep::setBz(double Bz) { Bz_m = Bz; }
