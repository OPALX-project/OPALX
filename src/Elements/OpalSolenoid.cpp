//
// Class OpalSolenoid
//   The SOLENOID element.
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
#include "Elements/OpalSolenoid.h"

#include <memory>
#include "AbsBeamline/Solenoid.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"

OpalSolenoid::OpalSolenoid()
    : OpalElement(SIZE, "SOLENOID", "The \"SOLENOID\" element defines a Solenoid.") {
    itsAttr[KS]     = Attributes::makeReal("KS", "Normalised solenoid strength in m**(-1)");
    itsAttr[DKS]    = Attributes::makeReal("DKS", "Normalised solenoid strength error in m**(-1)");
    itsAttr[FMAPFN] = Attributes::makeString("FMAPFN", "Solenoid field map filename ");
    itsAttr[FAST]   = Attributes::makeBool("FAST", "Faster but less accurate", true);
    itsAttr[ZREVERSE] = Attributes::makeBool(
            "ZREVERSE",
            "Read the field map back to front: mirror it in z and negate Bz, which is what "
            "turning the magnet around does. Only G4beamline cylinder maps support it.",
            false);

    registerOwnership();
}

OpalSolenoid::~OpalSolenoid() {}

std::shared_ptr<ElementBase> OpalSolenoid::build(const ObjectValues& values, double /*p0*/) {
    auto element = std::make_shared<Solenoid>(values.getName());
    fill(*element, values);
    return element;
}

void OpalSolenoid::fill(Solenoid& element, const ObjectValues& values) {
    setCommonFields(element, values);

    double length = values.getReal("L");
    bool fast     = values.getBool("FAST");

    element.getGeometry().setElementLength(length);
    element.setFieldMapFN(values.getString("FMAPFN"));
    element.setFast(fast);
    element.setIsZReversed(values.getBool("ZREVERSE"));
    element.setKS(values.getReal("KS"));
    element.setDKS(values.getReal("DKS"));
}
