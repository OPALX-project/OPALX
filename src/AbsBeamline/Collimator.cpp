//
// Class Collimator
//   Interface for a transverse-aperture collimator.
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
#include "AbsBeamline/Collimator.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "PartBunch/PartBunch.h"

Collimator::Collimator() : Collimator("") {}

Collimator::Collimator(const Collimator& right) : ElementBase(right) {}

Collimator::Collimator(const std::string& name) : ElementBase(name) {}

Collimator::~Collimator() {}

void Collimator::accept(BeamlineVisitor& visitor) const { visitor.visitCollimator(*this); }

void Collimator::initialise(PartBunch_t* bunch) { RefPartBunch_m = bunch; }

void Collimator::finalise() {}

void Collimator::getFieldExtent(double& zBegin, double& zEnd) const {
    // Local-chart support interval (field-free; report the body span).
    zBegin = 0.0;
    zEnd   = getGeometry().getElementLength();
}

ElementType Collimator::getType() const { return ElementType::COLLIMATOR; }
