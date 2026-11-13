//
// Class OpalBox
//   The BOX element.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Elements/OpalBox.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/BoxRep.h"
#include "Utilities/OpalException.h"

#include <algorithm>

OpalBox::OpalBox()
    : OpalElement(
              SIZE, "BOX",
              "The \"BOX\" element is a solid block that deletes every particle inside it.") {
    itsAttr[WIDTH]  = Attributes::makeReal("WIDTH", "Full width of the block along local x [m]");
    itsAttr[HEIGHT] = Attributes::makeReal("HEIGHT", "Full height of the block along local y [m]");

    registerOwnership();

    setElement(new BoxRep("BOX"));
}

OpalBox::OpalBox(const std::string& name, OpalBox* parent) : OpalElement(name, parent) {
    setElement(new BoxRep(name));
}

OpalBox::~OpalBox() {}

OpalBox* OpalBox::clone(const std::string& name) { return new OpalBox(name, this); }

void OpalBox::update() {
    // Applies the common attributes, including DELETEONTRANSVERSEEXIT, to the
    // embedded element.
    OpalElement::update();

    const double length = Attributes::getReal(itsAttr[LENGTH]);
    const double width  = Attributes::getReal(itsAttr[WIDTH]);
    const double height = Attributes::getReal(itsAttr[HEIGHT]);

    // The block sets the element's size. With an APERTURE the shared rule would also
    // delete particles between the APERTURE and the reach, and REACH is set below.
    if (itsAttr[APERT] && !Attributes::getString(itsAttr[APERT]).empty()) {
        throw OpalException(
                "OpalBox::update()",
                "BOX \"" + getOpalName() + "\" takes WIDTH and HEIGHT, not APERTURE");
    }
    if (itsAttr[REACH]) {
        throw OpalException(
                "OpalBox::update()",
                "BOX \"" + getOpalName() + "\" takes WIDTH and HEIGHT, not REACH");
    }
    if (width <= 0.0 || height <= 0.0 || length <= 0.0) {
        throw OpalException(
                "OpalBox::update()",
                "BOX \"" + getOpalName() + "\" requires WIDTH > 0, HEIGHT > 0 and L > 0");
    }

    BoxRep* box = dynamic_cast<BoxRep*>(getElement());
    box->getGeometry().setElementLength(length);
    box->setHalfSize(0.5 * width, 0.5 * height);
    // The reach is one half-width for x and y: the larger of the two covers the block.
    box->setReach(0.5 * std::max(width, height));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(box);
}
