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

#include <memory>

#include "Attributes/Attributes.h"
#include "BeamlineCore/FieldmapElementRep.h"
#include "Utilities/OpalException.h"

OpalFieldmapElement::OpalFieldmapElement()
    : OpalElement(
              SIZE, "FIELDMAP",
              "The \"FIELDMAP\" element defines an element whose field is a tabulated map.") {
    itsAttr[FMAPFN] = Attributes::makeString("FMAPFN", "Field map filename");
    itsAttr[BSCALE] = Attributes::makeReal(
            "BSCALE",
            "Multiplier applied to the tabulated magnetic field. G4beamline maps are read in "
            "absolute Tesla, so 1 reproduces the map as written; this is the current= given "
            "where the map is placed in a G4beamline input. OPAL and ASTRA maps are normalised "
            "to an on-axis peak of 1 T unless their header turns that off; for them BSCALE is "
            "the peak field in Tesla.",
            1.0);
    itsAttr[ESCALE] = Attributes::makeReal(
            "ESCALE",
            "Multiplier applied to the tabulated electric field, for a map that carries one. "
            "Separate from BSCALE because G4beamline scales the two fields independently: "
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
}

OpalFieldmapElement::OpalFieldmapElement(const std::string& name, OpalFieldmapElement* parent)
    : OpalElement(name, parent) {}

OpalFieldmapElement::~OpalFieldmapElement() {}

OpalFieldmapElement* OpalFieldmapElement::clone(const std::string& name) {
    return new OpalFieldmapElement(name, this);
}

std::shared_ptr<ElementBase> OpalFieldmapElement::build(const ElementValues& values) {
    auto element = std::make_shared<FieldmapElementRep>(values.getName());
    fill(*element, values);
    return element;
}

void OpalFieldmapElement::fill(ElementBase& element, const ElementValues& values) {
    // Checked before setCommonFields(), so that the specific message wins over the generic
    // "placement is over-specified" one.
    if (values.hasValue("ELEMEDGE")) {
        throw OpalException(
                "OpalFieldmapElement::build",
                values.getName()
                        + ": a FIELDMAP element cannot be placed with ELEMEDGE. ELEMEDGE "
                          "positions an element by path length s along the reference orbit, "
                          "which needs the element's length along that orbit; a FIELDMAP "
                          "element's extent comes from its field map and is not a path "
                          "length. Use an absolute 6D lab pose (X, Y, Z, THETA, PHI, PSI) "
                          "instead. OPALX requires one placement convention for the whole "
                          "beamline, so every other element in the line must use a 6D pose "
                          "too.");
    }

    if (values.hasValue("L")) {
        throw OpalException(
                "OpalFieldmapElement::build",
                values.getName()
                        + ": a FIELDMAP element takes no L. Its length is the longitudinal "
                          "extent of the field map named by FMAPFN, and is read from the map "
                          "header.");
    }

    if (!values.hasValue("FMAPFN") || values.getString("FMAPFN").empty()) {
        throw OpalException(
                "OpalFieldmapElement::build",
                values.getName() + ": a FIELDMAP element requires FMAPFN.");
    }

    setCommonFields(element, values);

    FieldmapElementRep* fm = dynamic_cast<FieldmapElementRep*>(&element);

    // The geometry length is deliberately NOT set here: it comes from the map, which is only
    // loaded in FieldmapElement::initialise(). That runs on the tracked clone before the
    // element list is sorted and before PlacementResolver, so every consumer still sees it.
    fm->setFieldMapFN(values.getString("FMAPFN"));
    fm->setBScale(values.getReal("BSCALE"));
    fm->setEScale(values.getReal("ESCALE"));
    fm->setIsZReversed(values.getBool("ZREVERSE"));
}
