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

#include "Elements/OpalMultipoleT.h"

#include <memory>
#include <vector>
#include "Attributes/Attributes.h"
#include "BeamlineCore/MultipoleTRep.h"

OpalMultipoleT::OpalMultipoleT()
    : OpalElement(
              SIZE, "MULTIPOLET",
              "The \"MULTIPOLET\" element defines a combined function multipole.") {
    // Attributes for a straight multipole
    itsAttr[TP] = Attributes::makeRealArray(
            "TP",
            "Array of multipolar field strengths b_k. The field generated in the "
            "flat top is B = b_k x^k [T m^(-k)]");
    itsAttr[LFRINGE]   = Attributes::makeReal("LFRINGE", "The length of the left end field [m]");
    itsAttr[RFRINGE]   = Attributes::makeReal("RFRINGE", "The length of the right end field [m]");
    itsAttr[HAPERT]    = Attributes::makeReal("HAPERT", "The aperture width [m]");
    itsAttr[VAPERT]    = Attributes::makeReal("VAPERT", "The aperture height [m]");
    itsAttr[MAXFORDER] = Attributes::makeReal(
            "MAXFORDER", "Number of terms used in each fringe component", DefaultMAXFORDER);
    itsAttr[ROTATION] = Attributes::makeReal(
            "ROTATION", "Rotation angle about its axis for skew elements [rad]");
    itsAttr[EANGLE] = Attributes::makeReal("EANGLE", "The entrance angle [rad]");
    itsAttr[BBLENGTH] =
            Attributes::makeReal("BBLENGTH", "Distance between centre of magnet and entrance [m]");

    // Further attributes for a constant radius curved multipole
    itsAttr[ANGLE] =
            Attributes::makeReal("ANGLE", "The azimuthal angle of the magnet in ring [rad]", 0.0);
    itsAttr[MAXXORDER] = Attributes::makeReal(
            "MAXXORDER", "Number of terms used in polynomial expansions", DefaultMAXXORDER);

    // Further attributes for a variable radius multipole
    itsAttr[VARRADIUS] =
            Attributes::makeBool("VARRADIUS", "Set true if radius of magnet is variable", false);
    itsAttr[ENTRYOFFSET] = Attributes::makeReal(
            "ENTRYOFFSET", "Longitudinal offset from standard entrance point [m]", 0.0);

    // Time dependence attributes
    itsAttr[SCALING_MODEL] = Attributes::makeUpperCaseString(
            "SCALING_MODEL",
            "The name of the time dependence model, which should give a scaling factor.");

    registerOwnership();
}

OpalMultipoleT::OpalMultipoleT(const std::string& name, OpalMultipoleT* parent)
    : OpalElement(name, parent) {}

OpalMultipoleT* OpalMultipoleT::clone(const std::string& name) {
    return new OpalMultipoleT(name, this);
}

void OpalMultipoleT::print(std::ostream& os) const { OpalElement::print(os); }

std::shared_ptr<ElementBase> OpalMultipoleT::build(const ElementValues& values) {
    auto element = std::make_shared<MultipoleT>(values.getName());
    fill(*element, values);
    return element;
}

void OpalMultipoleT::fill(ElementBase& element, const ElementValues& values) {
    // Base class first
    setCommonFields(element, values);
    // Make some sanity checks
    const auto maxFOrder = values.getReal("MAXFORDER");
    if (maxFOrder < MinimumMAXFORDER) {
        throw OpalException("OpalMultipoleT::build", "Attribute MAXFORDER must be >= 1.0");
    }
    if (maxFOrder > MaximumMAXFORDER) {
        throw OpalException("OpalMultipoleT::build", "Attribute MAXFORDER must be <= 9.0");
    }
    const double rotation  = values.getReal("ROTATION");
    const double bendAngle = values.getReal("ANGLE");
    if (bendAngle != 0.0 && rotation != 0.0) {
        throw OpalException(
                "OpalMultipoleT::build",
                "Non-zero ROTATION (a skew multipole) is only "
                "supported for straight magnets");
    }
    const bool varRadius = values.getBool("VARRADIUS");
    if (varRadius && bendAngle != 0.0) {
        throw OpalException("OpalMultipoleT::build", "Variable radius is not currently supported");
    }
    const double entryOffset = values.getReal("ENTRYOFFSET");
    if (!varRadius && entryOffset != 0.0) {
        throw OpalException(
                "OpalMultipoleT::build",
                "The ENTRYOFFSET is only supported for variable radius curved magnets");
    }
    auto tp = values.getRealArray("TP");
    // Set the attributes
    const auto length = values.getReal("L");
    if (bendAngle != 0.0 && length == 0.0) {
        throw OpalException(
                "OpalMultipoleT::build",
                "A curved magnet (non-zero ANGLE) needs a non-zero L: the curvature is "
                "ANGLE / L, with L the arc length");
    }
    auto* multT = dynamic_cast<MultipoleTRep*>(&element);
    multT->setElementLength(length);
    multT->setBendAngle(bendAngle, varRadius);
    multT->setAperture(values.getReal("VAPERT"), values.getReal("HAPERT"));
    multT->setFringeField(length * 0.5, values.getReal("LFRINGE"), values.getReal("RFRINGE"));
    multT->setBoundingBoxLength(values.getReal("BBLENGTH"));
    multT->setTransProfile(tp);
    multT->setMaxOrder(
            static_cast<size_t>(maxFOrder), static_cast<size_t>(values.getReal("MAXXORDER")));
    multT->setRotation(rotation);
    multT->setEntranceAngle(values.getReal("EANGLE"));
    multT->setEntryOffset(entryOffset);
    multT->setScalingName(values.getString("SCALING_MODEL"));
}
