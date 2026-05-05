//
// Class OpalBeamBeam
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
#include "Elements/OpalBeamBeam.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/BeamBeamRep.h"

OpalBeamBeam::OpalBeamBeam()
    : OpalElement(
              SIZE, "BEAMBEAM",
              "The \"BEAMBEAM\" element defines a beam-beam interaction point for colliding "
              "beams.") {
    itsAttr[GEOMETRY] =
            Attributes::makeString("GEOMETRY", "BoundaryGeometry for beam-beam elements");
    itsAttr[COPY] = Attributes::makeBool(
            "COPY",
            "If true, enable the future mirrored-bunch copy model for the beam-beam window.",
            false);
    itsAttr[VISUALIZE] = Attributes::makeBool(
            "VISUALIZE", "If true, emit the ASCII beam-beam-window visualization during tracking.",
            false);

    registerOwnership();

    setElement(new BeamBeamRep("DRIFT"));
}

OpalBeamBeam::OpalBeamBeam(const std::string& name, OpalBeamBeam* parent)
    : OpalElement(name, parent) {
    setElement(new BeamBeamRep(name));
}

OpalBeamBeam::~OpalBeamBeam() {}

OpalBeamBeam* OpalBeamBeam::clone(const std::string& name) { return new OpalBeamBeam(name, this); }

bool OpalBeamBeam::isBeamBeam() const { return true; }

void OpalBeamBeam::update() {
    OpalElement::update();

    BeamBeamRep* beamBeam = static_cast<BeamBeamRep*>(getElement());
    beamBeam->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    beamBeam->setAttribute("COPY", Attributes::getBool(itsAttr[COPY]) ? 1.0 : 0.0);
    beamBeam->setAttribute("VISUALIZE", Attributes::getBool(itsAttr[VISUALIZE]) ? 1.0 : 0.0);
    beamBeam->setAttribute("APERTURE_SET", itsAttr[APERT] ? 1.0 : 0.0);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(beamBeam);
}
