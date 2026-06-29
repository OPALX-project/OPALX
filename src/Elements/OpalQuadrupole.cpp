//
// Class OpalQuadrupole
//   The QUADRUPOLE element.
//
//   A convenience wrapper around MULTIPOLE that accepts a single
//   quadrupole strength K1 instead of the full KN vector.
//   Internally the element is represented as a MultipoleRep with
//   KN = {0, K1}.
//
// Copyright (c) 2024 – 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Elements/OpalQuadrupole.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MultipoleRep.h"
#include "Physics/Physics.h"

OpalQuadrupole::OpalQuadrupole()
    : OpalElement(
              SIZE, "QUADRUPOLE",
              "The \"QUADRUPOLE\" element defines a thick quadrupole.\n"
              "It is a convenience wrapper around MULTIPOLE with a single\n"
              "quadrupole strength K1.\n"
              "* If the length is non-zero, K1 is per unit length.\n"
              "* If the length is zero, K1 is the integrated strength.") {
    itsAttr[K1] = Attributes::makeReal("K1", "Normalised quadrupole strength in m^(-2)", 0.0);
    itsAttr[DK1] =
            Attributes::makeReal("DK1", "Normalised quadrupole strength error in m^(-2)", 0.0);
    itsAttr[K1S] =
            Attributes::makeReal("K1S", "Normalised skew quadrupole strength in m^(-2)", 0.0);
    itsAttr[DK1S] = Attributes::makeReal(
            "DK1S", "Normalised skew quadrupole strength error in m^(-2)", 0.0);

    registerOwnership();

    setElement(new MultipoleRep("QUADRUPOLE"));
}

OpalQuadrupole::OpalQuadrupole(const std::string& name, OpalQuadrupole* parent)
    : OpalElement(name, parent) {
    setElement(new MultipoleRep(name));
}

OpalQuadrupole::~OpalQuadrupole() {}

OpalQuadrupole* OpalQuadrupole::clone(const std::string& name) {
    return new OpalQuadrupole(name, this);
}

void OpalQuadrupole::print(std::ostream& os) const { OpalElement::print(os); }

void OpalQuadrupole::update() {
    OpalElement::update();

    MultipoleRep* mult = dynamic_cast<MultipoleRep*>(getElement());
    double length      = getLength();
    mult->setElementLength(length);

    // The KN vector convention is: index 0 = dipole, index 1 = quadrupole, ...
    // Populate the device coefficient views read by Multipole::apply().
    double k1   = Attributes::getReal(itsAttr[K1]);
    double dk1  = Attributes::getReal(itsAttr[DK1]);
    double k1s  = Attributes::getReal(itsAttr[K1S]);
    double dk1s = Attributes::getReal(itsAttr[DK1S]);

    // comp 0 (dipole): strength = 0
    mult->setNormalComponent(0, 0.0, 0.0);
    mult->setSkewComponent(0, 0.0, 0.0);

    // comp 1 (quadrupole)
    mult->setNormalComponent(1, k1, dk1);
    mult->setSkewComponent(1, k1s, dk1s);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(mult);
}
