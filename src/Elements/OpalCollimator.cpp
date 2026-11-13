//
// Class OpalCollimator
//   The COLLIMATOR element.
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
#include "Elements/OpalCollimator.h"

#include <memory>
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/CollimatorRep.h"
#include "Utilities/OpalException.h"

OpalCollimator::OpalCollimator()
    : OpalElement(
              COMMON, "COLLIMATOR",
              "The \"COLLIMATOR\" element deletes particles that leave its transverse aperture.") {
    registerOwnership();

    setElement(new CollimatorRep("COLLIMATOR"));
}

OpalCollimator::OpalCollimator(const std::string& name, OpalCollimator* parent)
    : OpalElement(name, parent) {
    setElement(new CollimatorRep(name));
}

OpalCollimator::~OpalCollimator() {}

OpalCollimator* OpalCollimator::clone(const std::string& name) {
    return new OpalCollimator(name, this);
}

std::shared_ptr<ElementBase> OpalCollimator::build(const ElementValues& values) {
    auto element = std::make_shared<CollimatorRep>(values.getName());
    fill(*element, values);
    return element;
}

void OpalCollimator::update() { fill(*getElement(), getValues()); }

void OpalCollimator::fill(ElementBase& element, const ElementValues& values) {
    // Applies the common attributes, including APERTURE and
    // DELETEONTRANSVERSEEXIT, to the embedded element.
    setCommonFields(element, values);

    // Unlike every other element, a collimator exists only to scrape, so
    // an absent aperture is an error rather than the wide-open default.
    if (!values.hasValue("APERTURE") || values.getString("APERTURE").empty()) {
        throw OpalException(
                "OpalCollimator::update()",
                "COLLIMATOR \"" + values.getName() + "\" requires an APERTURE");
    }

    // The aperture is only enforced for z in [0, L), so L = 0 would never
    // mark a particle.
    if (values.getReal("L") <= 0.0) {
        throw OpalException(
                "OpalCollimator::update()",
                "COLLIMATOR \"" + values.getName() + "\" requires L > 0");
    }

    const double length = values.getReal("L");

    CollimatorRep* coll = dynamic_cast<CollimatorRep*>(&element);
    coll->getGeometry().setElementLength(length);
}
