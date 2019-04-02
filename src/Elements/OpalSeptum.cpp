// ------------------------------------------------------------------------
// $RCSfile: OpalSeptum.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalSeptum
//   The class of OPAL Septums.
//
// ------------------------------------------------------------------------
//
// $Date: 2009/09/21 10:06:06 $
// $Author: bi $
//
// ------------------------------------------------------------------------

#include "Elements/OpalSeptum.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/SeptumRep.h"
#include "Structure/OpalWake.h"
#include "Physics/Physics.h"


// Class OpalSeptum
// ------------------------------------------------------------------------

OpalSeptum::OpalSeptum():
    OpalElement(SIZE, "SEPTUM",
                "The \"SEPTUM\" element defines a Septum."),
    owk_m(NULL) {

    itsAttr[XSTART] = Attributes::makeReal
                      ("XSTART", " Start of x coordinate");
    itsAttr[XEND] = Attributes::makeReal
                    ("XEND", " End of x coordinate");
    itsAttr[YSTART] = Attributes::makeReal
                      ("YSTART", "Start of y coordinate");
    itsAttr[YEND] = Attributes::makeReal
                    ("YEND", "End of y coordinate");
    itsAttr[WIDTH] = Attributes::makeReal
                     ("WIDTH", "Width of the septum");

    registerRealAttribute("XSTART");
    registerRealAttribute("XEND");
    registerRealAttribute("YSTART");
    registerRealAttribute("YEND");
    registerRealAttribute("WIDTH");

    registerOwnership();

    setElement((new SeptumRep("SEPTUM"))->makeAlignWrapper());
}


OpalSeptum::OpalSeptum(const std::string &name, OpalSeptum *parent):
    OpalElement(name, parent),
    owk_m(NULL) {
    setElement((new SeptumRep(name))->makeAlignWrapper());
}


OpalSeptum::~OpalSeptum() {
    if(owk_m)
        delete owk_m;
}


OpalSeptum *OpalSeptum::clone(const std::string &name) {
    return new OpalSeptum(name, this);
}


void OpalSeptum::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

}


void OpalSeptum::update() {
    OpalElement::update();

    SeptumRep *sept =
        dynamic_cast<SeptumRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double xstart = Attributes::getReal(itsAttr[XSTART]);
    double xend = Attributes::getReal(itsAttr[XEND]);
    double ystart = Attributes::getReal(itsAttr[YSTART]);
    double yend = Attributes::getReal(itsAttr[YEND]);
    double width = Attributes::getReal(itsAttr[WIDTH]);

    if(itsAttr[WAKEF] && owk_m == NULL) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + std::string("_wake"));
        owk_m->initWakefunction(*sept);
        sept->setWake(owk_m->wf_m);
    }
    sept->setElementLength(length);
    sept->setDimensions(xstart, xend, ystart, yend, width);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(sept);
}