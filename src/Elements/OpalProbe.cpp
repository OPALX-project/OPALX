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
    owk_m(nullptr) {

    itsAttr[XSTART] = Attributes::makeReal
                      ("XSTART", " Start of x coordinate [mm]");
    itsAttr[XEND] = Attributes::makeReal
                    ("XEND", " End of x coordinate [mm]");
    itsAttr[YSTART] = Attributes::makeReal
                      ("YSTART", "Start of y coordinate [mm]");
    itsAttr[YEND] = Attributes::makeReal
                    ("YEND", "End of y coordinate [mm]");
    itsAttr[WIDTH] = Attributes::makeReal
                     ("WIDTH", "Width of the probe, not used.");
    itsAttr[STEP] = Attributes::makeReal
                     ("STEP", "Step size of the probe [mm]", 1.0);

    registerRealAttribute("XSTART");
    registerRealAttribute("XEND");
    registerRealAttribute("YSTART");
    registerRealAttribute("YEND");
    registerRealAttribute("WIDTH");
    registerRealAttribute("STEP");

    registerOwnership();

    setElement((new ProbeRep("PROBE"))->makeAlignWrapper());
}


OpalProbe::OpalProbe(const std::string &name, OpalProbe *parent):
    OpalElement(name, parent),
    owk_m(nullptr) {
    setElement((new ProbeRep(name))->makeAlignWrapper());
}


OpalProbe::~OpalProbe() {
    delete owk_m;
}


OpalProbe *OpalProbe::clone(const std::string &name) {
    return new OpalProbe(name, this);
}


void OpalProbe::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
}


void OpalProbe::update() {
    OpalElement::update();

    ProbeRep *prob =
        dynamic_cast<ProbeRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double xstart = Attributes::getReal(itsAttr[XSTART]);
    double xend   = Attributes::getReal(itsAttr[XEND]);
    double ystart = Attributes::getReal(itsAttr[YSTART]);
    double yend   = Attributes::getReal(itsAttr[YEND]);
    double width  = Attributes::getReal(itsAttr[WIDTH]);
    double step   = Attributes::getReal(itsAttr[STEP]);

    if(itsAttr[WAKEF] && owk_m == nullptr) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + std::string("_wake"));
        owk_m->initWakefunction(*prob);
        prob->setWake(owk_m->wf_m);
    }
    prob->setElementLength(length); // is this needed here?
    prob->setDimensions(xstart, xend, ystart, yend, width);
    prob->setStep(step);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(prob);
}