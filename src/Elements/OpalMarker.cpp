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
    : OpalElement(COMMON, "MARKER", "The \"MARKER\" element defines a marker.") {
    // Construct the begin marker for beam lines.
    OpalMarker* S = new OpalMarker("#S", this);
    S->builtin    = true;
    OpalData::getInstance()->create(S);

    // Construct the end marker for beam lines.
    OpalMarker* E = new OpalMarker("#E", this);
    E->builtin    = true;
    OpalData::getInstance()->create(E);
}

OpalMarker::OpalMarker(const std::string& name, OpalMarker* parent) : OpalElement(name, parent) {}

OpalMarker::~OpalMarker() {}

OpalMarker* OpalMarker::clone(const std::string& name) { return new OpalMarker(name, this); }

void OpalMarker::print(std::ostream& os) const { OpalElement::print(os); }

std::shared_ptr<ElementBase> OpalMarker::build(const ElementValues& values) {
    return std::make_shared<MarkerRep>(values.getName());
}
