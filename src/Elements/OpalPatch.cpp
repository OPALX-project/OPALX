// ------------------------------------------------------------------------
// $RCSfile: OpalPatch.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalPatch
//   The class of OPAL patch spaces.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalPatch.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/PatchRep.h"


// Class OpalPatch
// ------------------------------------------------------------------------

OpalPatch::OpalPatch():
    OpalElement(SIZE, "PATCH",
                "The \"PATCH\" element implements a local reference change.") {
    itsAttr[LENGTH].setReadOnly(true);
    itsAttr[HORD] = Attributes::makeReal("HORD", "Horizontal orbit shift");
    itsAttr[VERTD] = Attributes::makeReal("VERTD", "Vertical orbit shift");
    itsAttr[LONGD] = Attributes::makeReal("LONGD", "Longitudinal orbit shift");
    itsAttr[VX] = Attributes::makeReal("VX", "Rotation around x-axis.");
    itsAttr[VY] = Attributes::makeReal("VY", "Rotation around y-axis.");
    itsAttr[VS] = Attributes::makeReal("VS", "Rotation around s-axis.");

    registerOwnership();

    setElement(new PatchRep("PATCH"));
}


OpalPatch::OpalPatch(const std::string &name, OpalPatch *parent):
    OpalElement(name, parent) {
    setElement(new PatchRep(name));
}


OpalPatch::~OpalPatch()
{}


OpalPatch *OpalPatch::clone(const std::string &name) {
    return new OpalPatch(name, this);
}


bool OpalPatch::isPatch() const {
    return true;
}


void OpalPatch::update() {
    OpalElement::update();

    PatchRep *patch = static_cast<PatchRep *>(getElement());
    double dx = Attributes::getReal(itsAttr[HORD]);
    double dy = Attributes::getReal(itsAttr[VERTD]);
    double ds = Attributes::getReal(itsAttr[LONGD]);
    double vx = Attributes::getReal(itsAttr[VX]);
    double vy = Attributes::getReal(itsAttr[VY]);
    double vs = Attributes::getReal(itsAttr[VS]);
    Euclid3D shift(dx, dy, ds, vx, vy, vs);
    patch->setPatch(shift);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(patch);
}