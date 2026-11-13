//
// Class OpalProbe
//   The Probe element.
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
#include "Elements/OpalProbe.h"

#include <memory>
#include "AbstractObjects/Attribute.h"
#include "Algorithms/CyclotronProbeRecorder.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/ProbeRep.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"

OpalProbe::OpalProbe() : OpalElement(SIZE, "PROBE", "The \"PROBE\" element defines a Probe.") {
    itsAttr[XSTART] = Attributes::makeReal("XSTART", "Start of laboratory X coordinate [mm]");
    itsAttr[XEND]   = Attributes::makeReal("XEND", "End of laboratory X coordinate [mm]");
    itsAttr[YSTART] =
            Attributes::makeReal("YSTART", "Start of laboratory Z coordinate (legacy Y name) [mm]");
    itsAttr[YEND] =
            Attributes::makeReal("YEND", "End of laboratory Z coordinate (legacy Y name) [mm]");
    itsAttr[WIDTH] = Attributes::makeReal("WIDTH", "Width of the probe, not used.");
    itsAttr[STEP]  = Attributes::makeReal("STEP", "Cyclotron radial histogram bin width [mm]", 1.0);

    registerOwnership();
}

OpalProbe::OpalProbe(const std::string& name, OpalProbe* parent) : OpalElement(name, parent) {}

OpalProbe::~OpalProbe() {}

OpalProbe* OpalProbe::clone(const std::string& name) { return new OpalProbe(name, this); }

std::shared_ptr<ElementBase> OpalProbe::build(const ElementValues& values) {
    auto element = std::make_shared<ProbeRep>(values.getName());
    fill(*element, values);
    return element;
}

void OpalProbe::fill(ElementBase& element, const ElementValues& values) {
    // Endpoint coordinates are already global. Ordinary element placement
    // validation would incorrectly require ELEMEDGE or a separate 6D pose.

    ProbeRep* prob = dynamic_cast<ProbeRep*>(&element);

    double xstart = Units::mm2m * values.getReal("XSTART");
    double xend   = Units::mm2m * values.getReal("XEND");
    double ystart = Units::mm2m * values.getReal("YSTART");
    double yend   = Units::mm2m * values.getReal("YEND");
    double step   = Units::mm2m * values.getReal("STEP");

    prob->getGeometry().setElementLength(0.0);
    prob->setDimensions(xstart, xend, ystart, yend);
    prob->setStep(step);
    prob->setOutputFN(values.getString("OUTFN"));
    CyclotronProbeRecorder::configuration(*prob);
}
