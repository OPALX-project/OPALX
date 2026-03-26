/*
 *  Copyright (c) 2017, Titus Dascalu
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include "Elements/OpalMultipoleT.h"
#include <iostream>
#include <vector>
#include "AbsBeamline/MultipoleT.h"
#include "Attributes/Attributes.h"

OpalMultipoleT::OpalMultipoleT()
    : OpalElement(
          SIZE, "MULTIPOLET", "The \"MULTIPOLET\" element defines a combined function multipole.") {
    // Attributes for a straight multipole
    itsAttr[TP] = Attributes::makeRealArray(
        "TP",
        "Array of multipolar field strengths b_k. The field generated in the "
        "flat top is B = b_k x^k [T m^(-k)]");
    itsAttr[LFRINGE] = Attributes::makeReal("LFRINGE", "The length of the left end field [m]");
    itsAttr[RFRINGE] = Attributes::makeReal("RFRINGE", "The length of the right end field [m]");
    itsAttr[HAPERT]  = Attributes::makeReal("HAPERT", "The aperture width [m]");
    itsAttr[VAPERT]  = Attributes::makeReal("VAPERT", "The aperture height [m]");
    itsAttr[MAXFORDER] = Attributes::makeReal(
        "MAXFORDER", "Number of terms used in each fringe component", DefaultMAXFORDER);
    itsAttr[ROTATION] = Attributes::makeReal(
        "ROTATION", "Rotation angle about its axis for skew elements [rad]");
    itsAttr[EANGLE] = Attributes::makeReal("EANGLE", "The entrance angle [rad]");
    itsAttr[BBLENGTH] = Attributes::makeReal(
        "BBLENGTH", "Distance between centre of magnet and entrance [m]");

    // Further attributes for a constant radius curved multipole
    itsAttr[ANGLE] = Attributes::makeReal(
        "ANGLE", "The azimuthal angle of the magnet in ring [rad]", 0.0);
    itsAttr[MAXXORDER] = Attributes::makeReal(
        "MAXXORDER", "Number of terms used in polynomial expansions", DefaultMAXXORDER);

    // Further attributes for a variable radius multipole
    itsAttr[VARRADIUS] = Attributes::makeBool(
        "VARRADIUS", "Set true if radius of magnet is variable", false);
    itsAttr[ENTRYOFFSET] = Attributes::makeReal(
        "ENTRYOFFSET", "Longitudinal offset from standard entrance point [m]", 0.0);

    // Time dependence attributes
    itsAttr[SCALING_MODEL] = Attributes::makeUpperCaseString("SCALING_MODEL",
         "The name of the time dependence model, which should give a scaling factor.");

    // Misalignment attributes
    itsAttr[MISALIGN_H] = Attributes::makeReal(
        "MISALIGN_H", "Horizontal misalignment [m]", 0.0);
    itsAttr[MISALIGN_V] = Attributes::makeReal(
        "MISALIGN_V", "Vertical misalignment [m]", 0.0);
    itsAttr[MISALIGN_S] = Attributes::makeReal(
        "MISALIGN_S", "Longitudinal misalignment [m]", 0.0);
    itsAttr[MISALIGN_ROLL] = Attributes::makeReal(
        "MISALIGN_ROLL", "Roll misalignment [rad] about the longitudinal axis", 0.0);
    itsAttr[MISALIGN_PITCH] = Attributes::makeReal(
        "MISALIGN_PITCH", "Pitch misalignment [rad] about the horizontal axis", 0.0);
    itsAttr[MISALIGN_YAW] = Attributes::makeReal(
        "MISALIGN_YAW", "Yaw misalignment [rad] about the vertical axis", 0.0);
    itsAttr[MISALIGN_AXISOFFSET] = Attributes::makeReal(
        "MISALIGN_AXISOFFSET", "Vertical offset of roll axis [m]", 0.0);

    registerOwnership();
    setElement(new MultipoleT("MULTIPOLET"));
}

OpalMultipoleT::OpalMultipoleT(const std::string& name, OpalMultipoleT* parent)
    : OpalElement(name, parent) {
    setElement(new MultipoleT(name));
}

OpalMultipoleT* OpalMultipoleT::clone(const std::string& name) {
    return new OpalMultipoleT(name, this);
}

void OpalMultipoleT::print(std::ostream& os) const {
    OpalElement::print(os);
}

void OpalMultipoleT::update() {
    // Base class first
    OpalElement::update();
    // Make some sanity checks
    const auto maxFOrder = Attributes::getReal(itsAttr[MAXFORDER]);
    if(maxFOrder < MinimumMAXFORDER) {
        throw OpalException("OpalMultipoleT::Update",
                            "Attribute MAXFORDER must be >= 1.0");
    }
    if(maxFOrder > MaximumMAXFORDER) {
        *gmsg << "OpalMultipoleT::Update, a value of "
                << maxFOrder << " for MAXFORDER may lead to excessive run time";
    }
    const double rotation = Attributes::getReal(itsAttr[ROTATION]);
    const double bendAngle = Attributes::getReal(itsAttr[ANGLE]);
    if(bendAngle != 0.0 && rotation != 0.0) {
        throw OpalException("OpalMultipoleT::Update",
                            "Non-zero ROTATION (a skew multipole) is only "
                            "supported for straight magnets");
    }
    const bool varRadius = Attributes::getBool(itsAttr[VARRADIUS]);
    if(varRadius && bendAngle != 0.0) {
        *gmsg << "OpalMultipoleT::Update, the variable radius "
                 "multipole magnet implementation is very slow";
    }
    const double entryOffset = Attributes::getReal(itsAttr[ENTRYOFFSET]);
    if((!varRadius || bendAngle == 0.0) && entryOffset != 0.0) {
        throw OpalException("OpalMultipoleT::Update",
                            "The ENTRYOFFSET is only supported for variable radius curved magnets");
    }
    // Convert pole strengths from Tesla to internal units which are kGauss
    auto tp = Attributes::getRealArray(itsAttr[TP]);
    for(auto& i : tp) {
        i *= Units::T2kG;
    }
    // Set the attributes
    const auto length = Attributes::getReal(itsAttr[LENGTH]);
    auto* multT = dynamic_cast<MultipoleT*>(getElement());
    multT->setElementLength(length);
    multT->setBendAngle(bendAngle, varRadius);
    multT->setAperture(Attributes::getReal(itsAttr[VAPERT]), Attributes::getReal(itsAttr[HAPERT]));
    multT->setFringeField(
        length * 0.5, Attributes::getReal(itsAttr[LFRINGE]), Attributes::getReal(itsAttr[RFRINGE]));
    multT->setBoundingBoxLength(Attributes::getReal(itsAttr[BBLENGTH]));
    multT->setTransProfile(tp);
    multT->setMaxOrder(static_cast<size_t>(maxFOrder),
        static_cast<size_t>(Attributes::getReal(itsAttr[MAXXORDER])));
    multT->setRotation(rotation);
    multT->setEntranceAngle(Attributes::getReal(itsAttr[EANGLE]));
    multT->setEntryOffset(entryOffset);
    multT->setScalingName(Attributes::getString(itsAttr[SCALING_MODEL]));
    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(multT);
}
