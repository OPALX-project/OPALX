// ------------------------------------------------------------------------
// $RCSfile: OpalProbe.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalProbe
//   The class of OPAL Probes.
//
// ------------------------------------------------------------------------
//
// $Date: 2009/10/07 10:06:06 $
// $Author: bi $
//
// ------------------------------------------------------------------------

#include "Elements/OpalProbe.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/ProbeRep.h"
#include "Structure/OpalWake.h"
#include "Physics/Physics.h"


// Class OpalProbe
// ------------------------------------------------------------------------

OpalProbe::OpalProbe():
    OpalElement(SIZE, "PROBE",
                "The \"PROBE\" element defines a Probe."),
    owk_m(NULL) {

    itsAttr[XSTART] = Attributes::makeReal
                      ("XSTART", " Start of x coordinate ");
    itsAttr[XEND] = Attributes::makeReal
                    ("XEND", " End of x coordinate");
    itsAttr[YSTART] = Attributes::makeReal
                      ("YSTART", "Start of y coordinate");
    itsAttr[YEND1] = Attributes::makeReal
                     ("YEND1", "Not used now");
    itsAttr[YEND] = Attributes::makeReal
                    ("YEND", "End of y coordinate");
    itsAttr[WIDTH] = Attributes::makeReal
                     ("WIDTH", "Width of the probe");

    registerRealAttribute("XSTART");
    registerRealAttribute("XEND");
    registerRealAttribute("YSTART");
    registerRealAttribute("YEND1");
    registerRealAttribute("YEND");
    registerRealAttribute("WIDTH");

    setElement((new ProbeRep("PROBE"))->makeAlignWrapper());
}


OpalProbe::OpalProbe(const string &name, OpalProbe *parent):
    OpalElement(name, parent),
    owk_m(NULL) {
    setElement((new ProbeRep(name))->makeAlignWrapper());
}


OpalProbe::~OpalProbe() {
    if(owk_m)
        delete owk_m;
}


OpalProbe *OpalProbe::clone(const string &name) {
    return new OpalProbe(name, this);
}


void OpalProbe::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

}


void OpalProbe::update() {
    ProbeRep *prob =
        dynamic_cast<ProbeRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double xstart = Attributes::getReal(itsAttr[XSTART]);
    double xend = Attributes::getReal(itsAttr[XEND]);
    double ystart = Attributes::getReal(itsAttr[YSTART]);
    double yend1 = Attributes::getReal(itsAttr[YEND1]);
    double yend = Attributes::getReal(itsAttr[YEND]);
    double width = Attributes::getReal(itsAttr[WIDTH]);

    if(itsAttr[WAKEF] && owk_m == NULL) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + string("_wake"));
        owk_m->initWakefunction(*prob);
        prob->setWake(owk_m->wf_m);
    }
    prob->setElementLength(length);
    prob->setXstart(xstart);
    prob->setXend(xend);
    prob->setYstart(ystart);
    prob->setYend(yend1);
    prob->setYend(yend);
    prob->setWidth(width);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(prob);
}
