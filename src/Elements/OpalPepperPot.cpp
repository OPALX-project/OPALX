// ------------------------------------------------------------------------
// $RCSfile: OpalPepperPot.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalPepperPot
//   The class of OPAL elliptic collimators.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalPepperPot.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/FlexibleCollimatorRep.h"
#include "Structure/ParticleMatterInteraction.h"


// Class OpalPepperPot
// ------------------------------------------------------------------------

OpalPepperPot::OpalPepperPot():
    OpalElement(SIZE, "PEPPERPOT",
                "The \"PEPPERPOT\" element defines an elliptic collimator."),
    parmatint_m(NULL) {
    itsAttr[XSIZE] = Attributes::makeReal
                     ("XSIZE", "Size in x of the pepperpot in m");
    itsAttr[YSIZE] = Attributes::makeReal
                     ("YSIZE", "Size in y of the pepperpot in m");
    itsAttr[OUTFN] = Attributes::makeString
                     ("OUTFN", "Pepperpot output filename");
    itsAttr[NHOLX] = Attributes::makeReal
                     ("NHOLX", "Number of holes in x");
    itsAttr[NHOLY] = Attributes::makeReal
                     ("NHOLY", "Number of holes in y");
    itsAttr[R] = Attributes::makeReal
                 ("R", "Radios of a holes in m");

    registerStringAttribute("OUTFN");
    registerRealAttribute("XSIZE");
    registerRealAttribute("YSIZE");
    registerRealAttribute("R");
    registerRealAttribute("NHOLX");
    registerRealAttribute("NHOLY");

    registerOwnership();

    setElement((new FlexibleCollimatorRep("PEPPERPOT"))->makeAlignWrapper());
}


OpalPepperPot::OpalPepperPot(const std::string &name, OpalPepperPot *parent):
    OpalElement(name, parent),
    parmatint_m(NULL) {
    setElement((new FlexibleCollimatorRep(name))->makeAlignWrapper());
}


OpalPepperPot::~OpalPepperPot() {
    if(parmatint_m)
        delete parmatint_m;
}


OpalPepperPot *OpalPepperPot::clone(const std::string &name) {
    return new OpalPepperPot(name, this);
}


void OpalPepperPot::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);


    // const FlexibleCollimatorRep *ppo =
    //     dynamic_cast<const FlexibleCollimatorRep *>(base.removeWrappers());
    // attributeRegistry["XSIZE"]->setReal(ppo->getXsize());
    // attributeRegistry["YSIZE"]->setReal(ppo->getYsize());

}

void OpalPepperPot::update() {
    OpalElement::update();

    FlexibleCollimatorRep *ppo =
        dynamic_cast<FlexibleCollimatorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    ppo->setElementLength(length);
    ppo->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

    if (getOpalName() != "PEPPERPOT") {
        double xsize = Attributes::getReal(itsAttr[XSIZE]);
        double ysize = Attributes::getReal(itsAttr[YSIZE]);
        double diameter = 2 * Attributes::getReal(itsAttr[R]);
        int repX = Attributes::getReal(itsAttr[NHOLX]) - 1;
        int repY = Attributes::getReal(itsAttr[NHOLY]) - 1;

        double shiftx = (xsize - diameter) / repX;
        double shifty = (ysize - diameter) / repY;

        std::stringstream description;
        description << "repeat(repeat(translate(ellipse("
                    << diameter << "," << diameter << "),"
                    << -shiftx * 0.5 * repX << "," << -shifty * 0.5 * repY << "),"
                    << repX << "," << shiftx << ",0.0),"
                    << repY << ",0.0," << shifty << ")";

        std::cout << "OpalPepperPot.cpp: " << __LINE__ << "\t"
                  << description.str() << std::endl;
        ppo->setDescription(description.str());
        exit(1);
    }

    if(itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == NULL) {
        parmatint_m = (ParticleMatterInteraction::find(Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION])))->clone(getOpalName() + std::string("_parmatint"));
        parmatint_m->initParticleMatterInteractionHandler(*ppo);
        ppo->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(ppo);
}