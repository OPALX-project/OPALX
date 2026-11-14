// ------------------------------------------------------------------------
// $RCSfile: OpalMarker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.4 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalMarker
//   The class of OPAL markers.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/01/22 15:16:02 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Elements/OpalMarker.h"

#include <memory>
#include "AbstractObjects/OpalData.h"
#include "BeamlineCore/MarkerRep.h"
#include "Utilities/Options.h"

#include <iostream>

// Class OpalMarker
// ------------------------------------------------------------------------

OpalMarker::OpalMarker()
    : OpalElement(COMMON, "MARKER", "The \"MARKER\" element defines a marker.") {}

void OpalMarker::defineBuiltins() {
    for (const char* name : {"#S", "#E"}) {
        OpalMarker* marker = new OpalMarker();
        marker->setOpalName(name);
        marker->builtin = true;
        OpalData::getInstance()->create(marker);
    }
}

OpalMarker::OpalMarker(const std::string& name, OpalMarker* parent) : OpalElement(name, parent) {}

OpalMarker::~OpalMarker() {}

OpalMarker* OpalMarker::clone(const std::string& name) { return new OpalMarker(name, this); }

void OpalMarker::print(std::ostream& os) const { OpalElement::print(os); }

std::shared_ptr<ElementBase> OpalMarker::build(const ElementValues& values) {
    return std::make_shared<MarkerRep>(values.getName());
}
