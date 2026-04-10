//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "Elements/OpalPolynomialTimeDependence.h"
#include <string>
#include "Algorithms/PolynomialTimeDependence.h"
#include "Attributes/Attributes.h"

const std::string OpalPolynomialTimeDependence::doc_string =
        "The \"POLYNOMIAL_TIME_DEPENDENCE\" element defines "
        "polynomial coefficients for time dependent RF phase, "
        "frequency, amplitude, etc, given by "
        "f(t) = P0+P1*t+P2*t^2+P3*t^3 where t is the time in ns";

// Use either P0..P3 or COEFFICIENTS.  If COEFFICIENTS are present, P0..P3 are ignored.
OpalPolynomialTimeDependence::OpalPolynomialTimeDependence()
    : OpalElement(static_cast<int>(SIZE), "POLYNOMIAL_TIME_DEPENDENCE", doc_string.c_str()) {
    itsAttr[P0] = Attributes::makeReal("P0", "Constant term in the polynomial expansion.");
    itsAttr[P1] =
            Attributes::makeReal("P1", "First order (linear) term in the polynomial expansion.");
    itsAttr[P2] = Attributes::makeReal(
            "P2", "Second order (quadratic) term in the polynomial expansion.");
    itsAttr[P3] =
            Attributes::makeReal("P3", "Third order (cubic) term in the polynomial expansion.");
    itsAttr[COEFFICIENTS] = Attributes::makeRealArray(
            "COEFFICIENTS", "Polynomial coefficients as an array with arbitrary length.");

    registerOwnership();
}

OpalPolynomialTimeDependence* OpalPolynomialTimeDependence::clone(const std::string& name) {
    return new OpalPolynomialTimeDependence(name, this);
}

void OpalPolynomialTimeDependence::print(std::ostream& out) const { OpalElement::print(out); }

OpalPolynomialTimeDependence::OpalPolynomialTimeDependence(
        const std::string& name, OpalPolynomialTimeDependence* parent)
    : OpalElement(name, parent) {}

void OpalPolynomialTimeDependence::update() {
    // getOpalName() comes from AbstractObjects/Object.h
    if (itsAttr[COEFFICIENTS] && (itsAttr[P0] || itsAttr[P1] || itsAttr[P2] || itsAttr[P3])) {
        // Only use P0..P3 or COEFFICIENTS, not both
        throw std::invalid_argument(
                "OpalPolynomialTimeDependence::Update: "
                "Use P0..P3 or COEFFICIENTS to specify the coefficients, not both.");
    }
    std::vector<double> polynomial_coefficients = Attributes::getRealArray(itsAttr[COEFFICIENTS]);
    if (polynomial_coefficients.empty()) {
        polynomial_coefficients.push_back(Attributes::getReal(itsAttr[P0]));
        polynomial_coefficients.push_back(Attributes::getReal(itsAttr[P1]));
        polynomial_coefficients.push_back(Attributes::getReal(itsAttr[P2]));
        polynomial_coefficients.push_back(Attributes::getReal(itsAttr[P3]));
    }
    AbstractTimeDependence::setTimeDependence(
            getOpalName(), std::make_shared<PolynomialTimeDependence>(polynomial_coefficients));
}
