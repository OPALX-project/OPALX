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

#include <memory>
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
}

OpalBox::OpalBox(const std::string& name, OpalBox* parent) : OpalElement(name, parent) {}

OpalBox::~OpalBox() {}

OpalBox* OpalBox::clone(const std::string& name) { return new OpalBox(name, this); }

std::shared_ptr<ElementBase> OpalBox::build(const ElementValues& values) {
    auto element = std::make_shared<BoxRep>(values.getName());
    fill(*element, values);
    return element;
}

void OpalBox::fill(ElementBase& element, const ElementValues& values) {
    // Applies the common attributes, including DELETEONTRANSVERSEEXIT, to the
    // embedded element.
    setCommonFields(element, values);

    const double length = values.getReal("L");
    const double width  = values.getReal("WIDTH");
    const double height = values.getReal("HEIGHT");

    // The block sets the element's size. With an APERTURE the shared rule would also
    // delete particles between the APERTURE and the reach, and REACH is set below.
    if (values.hasValue("APERTURE") && !values.getString("APERTURE").empty()) {
        throw OpalException(
                "OpalBox::build",
                "BOX \"" + values.getName() + "\" takes WIDTH and HEIGHT, not APERTURE");
    }
    if (values.hasValue("REACH")) {
        throw OpalException(
                "OpalBox::build",
                "BOX \"" + values.getName() + "\" takes WIDTH and HEIGHT, not REACH");
    }
    if (width <= 0.0 || height <= 0.0 || length <= 0.0) {
        throw OpalException(
                "OpalBox::build",
                "BOX \"" + values.getName() + "\" requires WIDTH > 0, HEIGHT > 0 and L > 0");
    }

    BoxRep* box = dynamic_cast<BoxRep*>(&element);
    box->getGeometry().setElementLength(length);
    box->setHalfSize(0.5 * width, 0.5 * height);
    // The reach is one half-width for x and y: the larger of the two covers the block.
    box->setReach(0.5 * std::max(width, height));
}
