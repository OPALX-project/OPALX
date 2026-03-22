//
// Class OpalIp
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
#include "Elements/OpalIp.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/IpRep.h"

OpalIp::OpalIp()
    : OpalElement(
          SIZE, "IP", "The \"IP\" element defines an interaction point for colliding beams.") {
    itsAttr[GEOMETRY] = Attributes::makeString("GEOMETRY", "BoundaryGeometry for Ips");

    itsAttr[COLWINLEN] =
        Attributes::makeReal("COLWINLEN", "The collision windlow lenght in [m]", 0.001);
    itsAttr[COPY] = Attributes::makeBool(
        "COPY",
        "If true, enable the future mirrored-bunch copy model for the interaction window.",
        false);
    itsAttr[VISUALIZE] = Attributes::makeBool(
        "VISUALIZE",
        "If true, emit the ASCII interaction-window visualization during tracking.",
        false);

    registerOwnership();

    setElement(new IpRep("DRIFT"));
}

OpalIp::OpalIp(const std::string& name, OpalIp* parent) : OpalElement(name, parent) {
    setElement(new IpRep(name));
}

OpalIp::~OpalIp() {
}

OpalIp* OpalIp::clone(const std::string& name) {
    return new OpalIp(name, this);
}

bool OpalIp::isIp() const {
    return true;
}

void OpalIp::update() {
    OpalElement::update();

    IpRep* ip = static_cast<IpRep*>(getElement());
    ip->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    ip->setAttribute("COLWINLEN", Attributes::getReal(itsAttr[COLWINLEN]));
    ip->setAttribute("COPY", Attributes::getBool(itsAttr[COPY]) ? 1.0 : 0.0);
    ip->setAttribute("VISUALIZE", Attributes::getBool(itsAttr[VISUALIZE]) ? 1.0 : 0.0);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(ip);
}
