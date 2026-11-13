//
// Class OpalMonitor
//   The MONITOR element.
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
#include "Elements/OpalMonitor.h"

#include <memory>
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MonitorRep.h"
#include "Utilities/Util.h"

OpalMonitor::OpalMonitor()
    : OpalElement(COMMON, "MONITOR", "The \"MONITOR\" element defines a monitor for both planes.") {
    registerOwnership();
}

OpalMonitor::OpalMonitor(const std::string& name, OpalMonitor* parent)
    : OpalElement(name, parent) {}

OpalMonitor::~OpalMonitor() {}

OpalMonitor* OpalMonitor::clone(const std::string& name) { return new OpalMonitor(name, this); }

std::shared_ptr<ElementBase> OpalMonitor::build(const ElementValues& values) {
    auto element = std::make_shared<MonitorRep>(values.getName());
    fill(*element, values);
    return element;
}

void OpalMonitor::fill(ElementBase& element, const ElementValues& values) {
    setCommonFields(element, values);

    MonitorRep* mon = dynamic_cast<MonitorRep*>(&element);

    double length = std::max(0.01, values.getReal("L"));
    mon->getGeometry().setElementLength(length);
    mon->setOutputFN(values.getString("OUTFN"));

    if (values.getString("TYPE") == "TEMPORAL") {
        mon->setCollectionType(CollectionType::TEMPORAL);
    } else {
        mon->setCollectionType(CollectionType::SPATIAL);
    }
}
