//
// Class OpalFieldmapElement
//   The FIELDMAP element.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#include "Elements/OpalFieldmapElement.h"

#include "Attributes/Attributes.h"
#include "BeamlineCore/FieldmapElementRep.h"
#include "Utilities/OpalException.h"

OpalFieldmapElement::OpalFieldmapElement()
    : OpalElement(
              SIZE, "FIELDMAP",
              "The \"FIELDMAP\" element defines an element whose field is a tabulated map.") {
    itsAttr[FMAPFN] = Attributes::makeString("FMAPFN", "Field map filename");
    itsAttr[SCALE]  = Attributes::makeReal(
            "SCALE",
            "Multiplier applied to the tabulated magnetic field. This is a plain factor, not "
             "a normalised strength: the G4beamline readers store absolute Tesla and do not "
             "normalise, so 1 reproduces the map as written. Same quantity as the current= "
             "given where the map is placed in a G4beamline input.",
            1.0);
    itsAttr[ESCALE] = Attributes::makeReal(
            "ESCALE",
            "Multiplier applied to the tabulated electric field, for a map that carries one. "
            "Separate from SCALE because G4beamline scales the two fields independently: "
            "current and normB for the magnetic field, gradient and normE for the electric "
            "one. Same quantity as the gradient= given where the map is placed. No effect on "
            "a map with no electric field.",
            1.0);
    itsAttr[ZREVERSE] = Attributes::makeBool(
            "ZREVERSE",
            "Read the field map back to front: mirror it in z and negate the longitudinal "
            "component, which is what turning the magnet around does. Only G4beamline cylinder "
            "maps support it.",
            false);

    registerOwnership();

    setElement(new FieldmapElementRep("FIELDMAP"));
}

OpalFieldmapElement::OpalFieldmapElement(const std::string& name, OpalFieldmapElement* parent)
    : OpalElement(name, parent) {
    setElement(new FieldmapElementRep(name));
}

OpalFieldmapElement::~OpalFieldmapElement() {}

OpalFieldmapElement* OpalFieldmapElement::clone(const std::string& name) {
    return new OpalFieldmapElement(name, this);
}

void OpalFieldmapElement::update() {
    // Checked before OpalElement::update(), so that the specific message wins over the
    // generic "placement is over-specified" one. OpalData::update() calls update() on every
    // object in the directory including the untouched builtin prototype, whose attributes are
    // all unset, so only real (cloned) instances are validated -- the same guard
    // OpalElement::validatePlacement() uses.
    if (!isBuiltin()) {
        if (itsAttr[ELEMEDGE]) {
            throw OpalException(
                    "OpalFieldmapElement::update()",
                    getOpalName()
                            + ": a FIELDMAP element cannot be placed with ELEMEDGE. ELEMEDGE "
                              "positions an element by path length s along the reference orbit, "
                              "which needs the element's length along that orbit; a FIELDMAP "
                              "element's extent comes from its field map and is not a path "
                              "length. Use an absolute 6D lab pose (X, Y, Z, THETA, PHI, PSI) "
                              "instead. OPALX requires one placement convention for the whole "
                              "beamline, so every other element in the line must use a 6D pose "
                              "too.");
        }

        if (itsAttr[LENGTH]) {
            throw OpalException(
                    "OpalFieldmapElement::update()",
                    getOpalName()
                            + ": a FIELDMAP element takes no L. Its length is the longitudinal "
                              "extent of the field map named by FMAPFN, and is read from the map "
                              "header.");
        }

        if (!itsAttr[FMAPFN] || Attributes::getString(itsAttr[FMAPFN]).empty()) {
            throw OpalException(
                    "OpalFieldmapElement::update()",
                    getOpalName() + ": a FIELDMAP element requires FMAPFN.");
        }
    }

    OpalElement::update();

    FieldmapElementRep* fm = dynamic_cast<FieldmapElementRep*>(getElement());

    // The geometry length is deliberately NOT set here: it comes from the map, which is only
    // loaded in FieldmapElement::initialise(). That runs on the tracked clone before the
    // element list is sorted and before PlacementResolver, so every consumer still sees it.
    fm->setFieldMapFN(Attributes::getString(itsAttr[FMAPFN]));
    fm->setScale(Attributes::getReal(itsAttr[SCALE]));
    fm->setEScale(Attributes::getReal(itsAttr[ESCALE]));
    fm->setIsZReversed(Attributes::getBool(itsAttr[ZREVERSE]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(fm);
}
