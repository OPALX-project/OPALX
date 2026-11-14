//
// Class OpalDrift
//   The class of OPAL drift spaces.
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
#include "Elements/OpalDrift.h"

#include <memory>
#include "Attributes/Attributes.h"
#include "BeamlineCore/DriftRep.h"

OpalDrift::OpalDrift()
    : OpalElement(SIZE, "DRIFT", "The \"DRIFT\" element defines a drift space.") {
    itsAttr[GEOMETRY] = Attributes::makeString("GEOMETRY", "BoundaryGeometry for Drifts");

    itsAttr[NSLICES] = Attributes::makeReal(
            "NSLICES", "The number of slices/ steps for this element in Map Tracking", 1);

    registerOwnership();
}

OpalDrift::~OpalDrift() {}

bool OpalDrift::isDrift() const { return true; }

std::shared_ptr<ElementBase> OpalDrift::build(const ElementValues& values) {
    auto element = std::make_shared<DriftRep>(values.getName());
    fill(*element, values);
    return element;
}

void OpalDrift::fill(ElementBase& element, const ElementValues& values) {
    setCommonFields(element, values);
    element.getGeometry().setElementLength(values.getReal("L"));
}
