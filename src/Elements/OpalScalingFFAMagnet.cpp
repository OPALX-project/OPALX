//
// Class OpalScalingFFAMagnet
//   The class provides the user interface for the SCALINGFFAMAGNET object.
//
// Copyright (c) 2017 - 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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
#include "Elements/OpalScalingFFAMagnet.h"

#include "AbsBeamline/EndFieldModel/EndFieldModelManager.h"
#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "AbsBeamline/ScalingFFAMagnet.h"
#include "Attributes/Attributes.h"
#include "Physics/Units.h"
#include "Utilities/OpalException.h"

OpalScalingFFAMagnet::OpalScalingFFAMagnet()
    : OpalElement(
              SIZE, "SCALINGFFAMAGNET",
              "The \"ScalingFFAMagnet\" element defines a FFA scaling magnet. Placement note: "
              "Because FFAs are normally circular, ScalingFFAMagnet elements that are placed "
              "using ELEMEDGE will be placed on a circle of radius R0 with ELEMEDGE = 0.0 "
              "corresponding to (x, y, z) = (0, 0, 0). If R0 is positive, the placements "
              "will have an anticlockwise arrangement; the ring curves towards the negative x "
              "direction. If R0 is negative placements will have a clockwise arrangement.") {
    itsAttr[B0] = Attributes::makeReal("B0", "The nominal dipole field of the magnet [T].");

    itsAttr[R0] = Attributes::makeReal("R0", "Radial scale [m].");

    itsAttr[FIELD_INDEX] = Attributes::makeReal("FIELD_INDEX", "The scaling magnet field index.");

    itsAttr[TAN_DELTA] = Attributes::makeReal(
            "TAN_DELTA", "Tangent of the spiral angle; set to 0 to make a radial sector magnet.");

    itsAttr[MAX_Y_POWER] = Attributes::makeReal(
            "MAX_Y_POWER",
            "The maximum power in y that will be considered in the field expansion (default 3).",
            3);

    itsAttr[END_FIELD_MODEL] = Attributes::makeString(
            "END_FIELD_MODEL",
            "Names the end field model of the magnet, giving the field magnitude along a line of "
            "constant radius. If blank, uses the 'FRINGE_LENGTH' and 'L' parameters to construct a  "
            "tanh model. If 'END_FIELD_MODEL' is not blank, OpalX will seek "
            "an END_FIELD_MODEL corresponding to the name defined in this string.");

    itsAttr[FRINGE_LENGTH] = Attributes::makeReal(
            "FRINGE_LENGTH",
            "The fringe field e-fold length of the spiral FFA, if END_FIELD_MODEL is not defined [m]. This "
            "determines the fringe field taper if no "
            "END_FIELD_MODEL is defined. Ignored if END_FIELD_MODEL is defined.");

    itsAttr[LENGTH] = Attributes::makeReal(
            "L",
            "The centre length of the spiral FFA, if END_FIELD_MODEL is not defined [m]. If"
            "END_FIELD_MODEL is defined, LENGTH is taken from the END_FIELD_MODEL.");

    itsAttr[RADIAL_NEG_EXTENT] = Attributes::makeReal(
            "RADIAL_NEG_EXTENT",
            "Particles are considered outside the tracking region if "
            "radius is less than R0-RADIAL_NEG_EXTENT relative to the FFA centre [m].",
            1);

    itsAttr[RADIAL_POS_EXTENT] = Attributes::makeReal(
            "RADIAL_POS_EXTENT",
            "Particles are considered outside the tracking region if "
            "radius is greater than R0+RADIAL_POS_EXTENT relative to the FFA centre [m].",
            1);

    itsAttr[HEIGHT] = Attributes::makeReal(
            "HEIGHT",
            "Full height of the magnet. Particles moving more than height/2. "
            "off the midplane (either above or below) are out of the aperture [m].", 1);

    itsAttr[AZIMUTHAL_EXTENT] = Attributes::makeReal(
            "AZIMUTHAL_EXTENT",
            "The field will be assumed zero if particles have S more than AZIMUTHAL_EXTENT "
            "from the magnet centre. Default is CENTRE_LENGTH/2.+5.*FRINGE_LENGTH [m].");

    registerOwnership();

    ScalingFFAMagnet* magnet = new ScalingFFAMagnet("ScalingFFAMagnet");
    magnet->setEndField(std::make_shared<endfieldmodel::Tanh>(1., 1., 1));
    setElement(magnet);
}

OpalScalingFFAMagnet::OpalScalingFFAMagnet(const std::string& name, OpalScalingFFAMagnet* parent)
    : OpalElement(name, parent) {
    ScalingFFAMagnet* magnet = new ScalingFFAMagnet(name);
    magnet->setEndField(std::make_shared<endfieldmodel::Tanh>(1., 1., 1));
    setElement(magnet);
}

OpalScalingFFAMagnet::~OpalScalingFFAMagnet() {}

OpalScalingFFAMagnet* OpalScalingFFAMagnet::clone(const std::string& name) {
    return new OpalScalingFFAMagnet(name, this);
}

void OpalScalingFFAMagnet::setupDefaultEndField() {
    ScalingFFAMagnet* magnet = dynamic_cast<ScalingFFAMagnet*>(getElement());
    // get centre length and end length in metres
    double end_length    = Attributes::getReal(itsAttr[FRINGE_LENGTH]);
    double x0 = Attributes::getReal(itsAttr[LENGTH]) / 2.;
    auto endField = std::make_shared<endfieldmodel::Tanh>();
    endField->setLambda(end_length);
    // x0 is the distance between B=0.5*B0 and B=B0 i.e. half the centre length
    endField->setX0(x0);
    magnet->setEndField(endField);
    std::string endName = "__opal_internal__" + getOpalName();
    magnet->setEndFieldName(endName);
    magnet->setPhiStart(0.0);
    endfieldmodel::EndFieldModelManager::getEFMManager()->setEndFieldModel(endName, endField);
}

void OpalScalingFFAMagnet::setupNamedEndField() {
    if (!itsAttr[END_FIELD_MODEL]) {
        return;
    }
    std::string name         = Attributes::getString(itsAttr[END_FIELD_MODEL]);
    ScalingFFAMagnet* magnet = dynamic_cast<ScalingFFAMagnet*>(getElement());
    magnet->setEndFieldName(name);
    magnet->setPhiStart(0.0);
}

void OpalScalingFFAMagnet::update() {
    OpalElement::update();
    ScalingFFAMagnet* magnet = dynamic_cast<ScalingFFAMagnet*>(getElement());
    // use L = r0*theta; we define the magnet into length for UI but into angles
    // internally; and use m as external default unit
    double r0Abs    = std::abs(Attributes::getReal(itsAttr[R0]));
    double r0Signed = Attributes::getReal(itsAttr[R0]);
    magnet->setR0(r0Signed);
    magnet->setDipoleConstant(Attributes::getReal(itsAttr[B0]));
    if (itsAttr[APERT]) {
        throw OpalException(
                "OpalScalingFFAMagnet::update()", "SCALINGFFAMAGNET does not use APERTURE command");
    }

    // dimensionless quantities
    magnet->setFieldIndex(Attributes::getReal(itsAttr[FIELD_INDEX]));
    magnet->setTanDelta(Attributes::getReal(itsAttr[TAN_DELTA]));
    int maxOrder = std::floor(Attributes::getReal(itsAttr[MAX_Y_POWER]));
    magnet->setMaxOrder(maxOrder);

    // get rmin and rmax bounding box edge
    double rmin = r0Abs - Attributes::getReal(itsAttr[RADIAL_NEG_EXTENT]);
    double rmax = r0Abs + Attributes::getReal(itsAttr[RADIAL_POS_EXTENT]);
    magnet->setRMin(rmin);
    magnet->setRMax(rmax);

    // we store maximum vertical displacement (which is half the height)
    double height = Attributes::getReal(itsAttr[HEIGHT]);
    magnet->setVerticalExtent(height / 2.);

    // get azimuthal extent in radians; this is just the bounding box
    if (itsAttr[AZIMUTHAL_EXTENT]) {
        if (Attributes::getReal(itsAttr[AZIMUTHAL_EXTENT]) < 0.0) {
            throw OpalException("OpalScalingFFAMagnet::update()", "AZIMUTHAL_EXTENT must be > 0.0");
        }
        magnet->setAzimuthalExtent(Attributes::getReal(itsAttr[AZIMUTHAL_EXTENT]) / r0Abs);
    } else {
        magnet->setAzimuthalExtent(-1);  // leave it for setupEndField
    }
    if (itsAttr[END_FIELD_MODEL]) {
        setupNamedEndField();
    } else {
        setupDefaultEndField();
    }
    magnet->initialise();
    setElement(magnet);
}
