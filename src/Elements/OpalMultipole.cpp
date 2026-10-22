//
// Class OpalMultipole
//   The MULTIPOLE element.
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
#include "Elements/OpalMultipole.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MultipoleRep.h"
#include "Expressions/SRefExpr.h"
#include "Expressions/SValue.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

OpalMultipole::OpalMultipole()
    : OpalElement(
              SIZE, "MULTIPOLE",
              "The \"MULTIPOLE\" element defines a straight thick multipole, from "
              "quadrupole up to decapole.\n"
              "* If the length is non-zero, the strengths are per unit "
              "length.\n* If the length is zero, the strengths are the "
              "values integrated over the length.\n"
              "* With zero length no synchrotron radiation can be calculated.") {
    itsAttr[KN] = Attributes::makeRealArray(
            "KN",
            "Normalised normal multipole strengths in m^(-k): {dipole, quadrupole, sextupole, "
            "octupole, decapole}. The first (dipole) entry must be 0; use SBEND or RBEND for a "
            "dipole field.");
    itsAttr[DKN] = Attributes::makeRealArray(
            "DKN", "Errors added to KN, in m^(-k). Same entries and rules as KN.");
    itsAttr[KS] = Attributes::makeRealArray(
            "KS", "Normalised skew multipole strengths in m^(-k). Same entries and rules as KN.");
    itsAttr[DKS] = Attributes::makeRealArray(
            "DKS", "Errors added to KS, in m^(-k). Same entries and rules as KN.");

    registerOwnership();

    setElement(new MultipoleRep("MULTIPOLE"));
}

OpalMultipole::OpalMultipole(const std::string& name, OpalMultipole* parent)
    : OpalElement(name, parent) {
    setElement(new MultipoleRep(name));
}

OpalMultipole::~OpalMultipole() {}

OpalMultipole* OpalMultipole::clone(const std::string& name) {
    return new OpalMultipole(name, this);
}

void OpalMultipole::print(std::ostream& os) const { OpalElement::print(os); }

void OpalMultipole::update() {
    OpalElement::update();

    // Get pointer to MultipoleRep object to set length
    MultipoleRep* mult = dynamic_cast<MultipoleRep*>(getElement());
    double length      = getLength();
    mult->getGeometry().setElementLength(length);

    // Get the vector with the multipole expansion components
    const std::vector<double> norm = Attributes::getRealArray(itsAttr[KN]);
    std::vector<double> normErrors = Attributes::getRealArray(itsAttr[DKN]);
    const std::vector<double> skew = Attributes::getRealArray(itsAttr[KS]);
    std::vector<double> skewErrors = Attributes::getRealArray(itsAttr[DKS]);

    // OpalData::update() calls update() on every object in the directory, including the
    // unmodified builtin prototype. Only check real (cloned) element instances.
    if (!isBuiltin()) {
        // MULTIPOLE is placed as a straight element, so a dipole entry would steer the beam away
        // from the line the element sits on. The errors are added to the strengths, so their
        // first entry would make a dipole too.
        const auto check = [this](const std::string& name, const std::vector<double>& values) {
            if (values.size() > static_cast<std::size_t>(MAX_MP_ORDER)) {
                throw OpalException(
                        "OpalMultipole::update",
                        "MULTIPOLE \"" + getOpalName() + "\": " + name + " has "
                                + std::to_string(values.size())
                                + " entries. MULTIPOLE goes up to decapole: at most 5 entries "
                                  "(dipole, quadrupole, sextupole, octupole, decapole).");
            }
            if (!values.empty() && values[0] != 0.0) {
                throw OpalException(
                        "OpalMultipole::update",
                        "MULTIPOLE \"" + getOpalName() + "\": the first " + name
                                + " entry is the dipole and must be 0. MULTIPOLE is a straight "
                                  "element; use SBEND or RBEND for a dipole field.");
            }
        };
        check("KN", norm);
        check("DKN", normErrors);
        check("KS", skew);
        check("DKS", skewErrors);
    }

    // Resize error arrays
    unsigned int normSize = norm.size();
    unsigned int skewSize = skew.size();
    normErrors.resize(normSize, 0.0);
    skewErrors.resize(skewSize, 0.0);

    unsigned int top = (normSize > skewSize) ? normSize : skewSize;

    // KN/KS are NORMALISED strengths [m^-(n+1)]; the components stored/applied by
    // Multipole are physical field coefficients [T/m^n]. Convert with Brho = P0/c,
    // one factor for every order, exactly as OpalSBend does (OpalSBend.cpp:71-76).
    double factor = OpalData::getInstance()->getP0() / Physics::c;

    // Loop over components (0=Dipole, 1=Quadrupole, ...) populating the device
    // coefficient views read by Multipole::apply().
    for (unsigned int comp = 0; comp < top; ++comp) {
        if (comp < normSize) {
            mult->setNormalComponent(comp, factor * norm[comp], factor * normErrors[comp]);
        }
        if (comp < skewSize) {
            mult->setSkewComponent(comp, factor * skew[comp], factor * skewErrors[comp]);
        }
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(mult);
}