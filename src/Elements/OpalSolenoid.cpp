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
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/SolenoidRep.h"
#include "Physics/Physics.h"

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

std::shared_ptr<ElementBase> OpalSolenoid::build(const ElementValues& values) {
    auto element = std::make_shared<SolenoidRep>(values.getName());
    fill(*element, values);
    return element;
}

void OpalSolenoid::fill(ElementBase& element, const ElementValues& values) {
    setCommonFields(element, values);

    SolenoidRep* sol = dynamic_cast<SolenoidRep*>(&element);
    double length    = values.getReal("L");
    double Bz        = values.getReal("KS") * OpalData::getInstance()->getP0() / Physics::c;
    bool fast        = values.getBool("FAST");

    sol->getGeometry().setElementLength(length);
    sol->setFieldMapFN(values.getString("FMAPFN"));
    sol->setFast(fast);
    sol->setIsZReversed(values.getBool("ZREVERSE"));
    sol->setBz(Bz);
    sol->setKS(values.getReal("KS"));
    sol->setDKS(values.getReal("DKS"));
}
