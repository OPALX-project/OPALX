// ------------------------------------------------------------------------
// $RCSfile: OpalStripper.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalStripper
//   The class of OPAL Strippers.
//
// ------------------------------------------------------------------------
//
// $Date: 2011/07/08 08:16:01 $
// $Author: Jianjun Yang $
//
// ------------------------------------------------------------------------

#include "Elements/OpalStripper.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/StripperRep.h"
#include "Physics/Physics.h"


// Class OpalStripper
// ------------------------------------------------------------------------

OpalStripper::OpalStripper():
    OpalElement(SIZE, "STRIPPER",
                "The \"STRIPPER\" element defines a Stripper.") {

    itsAttr[XSTART] = Attributes::makeReal
                      ("XSTART", " Start of x coordinate [mm]");
    itsAttr[XEND] = Attributes::makeReal
                    ("XEND", " End of x coordinate, [mm]");
    itsAttr[YSTART] = Attributes::makeReal
                      ("YSTART", "Start of y coordinate, [mm]");
    itsAttr[YEND] = Attributes::makeReal
                    ("YEND", "End of y coordinate, [mm]");
    itsAttr[WIDTH] = Attributes::makeReal
                     ("WIDTH", "Width of the stripper [mm], NOT used yet");
    itsAttr[OPCHARGE] = Attributes::makeReal
                     ("OPCHARGE", "Charge number of the outcome particle");
    itsAttr[OPMASS] = Attributes::makeReal
                     ("OPMASS", "Mass of the outcome particle [GeV/c^2]");
    itsAttr[OPYIELD] = Attributes::makeReal
                     ("OPYIELD", "Yield (Particle number of the outcome particle) per income particle");
    itsAttr[STOP] = Attributes::makeBool
      ("STOP", "Option Whether stop tracking at the stripper. Default: true", true);

    registerRealAttribute("XSTART");
    registerRealAttribute("XEND");
    registerRealAttribute("YSTART");
    registerRealAttribute("YEND");
    registerRealAttribute("WIDTH");
    registerRealAttribute("OPCHARGE");
    registerRealAttribute("OPMASS");
    registerRealAttribute("OPYIELD");

    setElement((new StripperRep("STRIPPER"))->makeAlignWrapper());
}


OpalStripper::OpalStripper(const string &name, OpalStripper *parent):
    OpalElement(name, parent) {
    setElement((new StripperRep(name))->makeAlignWrapper());
}


OpalStripper::~OpalStripper() 
{}


OpalStripper *OpalStripper::clone(const string &name) {
    return new OpalStripper(name, this);
}


void OpalStripper::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

}


void OpalStripper::update() {
    StripperRep *strp =dynamic_cast<StripperRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double xstart = Attributes::getReal(itsAttr[XSTART]);
    double xend = Attributes::getReal(itsAttr[XEND]);
    double ystart = Attributes::getReal(itsAttr[YSTART]);
    double yend = Attributes::getReal(itsAttr[YEND]);
    double width = Attributes::getReal(itsAttr[WIDTH]);
    double opcharge = Attributes::getReal(itsAttr[OPCHARGE]);
    double opmass = Attributes::getReal(itsAttr[OPMASS]);
    double opyield = Attributes::getReal(itsAttr[OPYIELD]);
    bool   stop = Attributes::getBool(itsAttr[STOP]);

    strp->setElementLength(length);
    strp->setXstart(xstart);
    strp->setXend(xend);
    strp->setYstart(ystart);
    strp->setYend(yend);
    strp->setWidth(width);
    strp->setOPCharge(opcharge);
    strp->setOPMass(opmass);
    strp->setOPYield(opyield);
    strp->setStop(stop);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(strp);
}
