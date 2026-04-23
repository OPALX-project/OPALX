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

#include "Elements/OpalSplineTimeDependence.h"
#include <string>
#include "Algorithms/SplineTimeDependence.h"
#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"

const std::string OpalSplineTimeDependence::doc_string =
        "The \"SPLINE_TIME_DEPENDENCE\" element defines "
        "an array of times and corresponding values for time lookup, "
        "for use in time-dependent elements. Lookup is supported at "
        "first order or third order with quadratic smoothing.";

// I investigated using a StringArray or RealArray here;
// Don't seem to have capacity to handle variables, so for now not implemented
OpalSplineTimeDependence::OpalSplineTimeDependence()
    : OpalElement(static_cast<int>(SIZE), "SPLINE_TIME_DEPENDENCE", doc_string.c_str()) {
    itsAttr[ORDER] = Attributes::makeReal(
            "ORDER",
            "Order of the lookup - either 1 for linear interpolation, "
            "or 3 for cubic interpolation with quadratic smoothing. "
            "Other values make an error.",
            1);

    itsAttr[TIMES] = Attributes::makeRealArray(
            "TIMES",
            "Array of real times in ns. There must be at least \"ORDER\"+1 "
            "elements in the array and they must be strictly monotically "
            "increasing");

    itsAttr[VALUES] = Attributes::makeRealArray(
            "VALUES",
            "Array of real values. The length of \"VALUES\" must be the "
            "same as the length of \"TIMES\".");

    registerOwnership();
}

OpalSplineTimeDependence* OpalSplineTimeDependence::clone(const std::string& name) {
    return new OpalSplineTimeDependence(name, this);
}

void OpalSplineTimeDependence::print(std::ostream& out) const { OpalElement::print(out); }

OpalSplineTimeDependence::OpalSplineTimeDependence(
        const std::string& name, OpalSplineTimeDependence* parent)
    : OpalElement(name, parent) {}

void OpalSplineTimeDependence::update() {
    const double order = Attributes::getReal(itsAttr[ORDER]);
    if (order != SplineTimeDependence::LinearInterpolation
        && order != SplineTimeDependence::CubicInterpolation) {
        throw OpalException(
                "OpalSplineTimeDependence::update",
                "SPLINE_TIME_DEPENDENCE \"ORDER\" should be 1 or 3.");
    }
    // Note we set array defaults as it seems that the default object must be valid
    std::vector<double> times;
    if (itsAttr[TIMES]) {
        times = Attributes::getRealArray(itsAttr[TIMES]);
    } else {
        times = {0.0, 1.0};  // A default time array
    }
    std::vector<double> values;
    if (itsAttr[VALUES]) {
        values = Attributes::getRealArray(itsAttr[VALUES]);
    } else {
        values = {0.0, 0.0};  // A default value array
    }
    const auto spline =
            std::make_shared<SplineTimeDependence>(static_cast<size_t>(order), times, values);
    AbstractTimeDependence::setTimeDependence(getOpalName(), spline);
}
