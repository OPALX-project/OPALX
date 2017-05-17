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
#include "BeamlineCore/CollimatorRep.h"
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
    itsAttr[PITCH] = Attributes::makeReal
                     ("PITCH", "Pitch of the pepperpot in m");
    itsAttr[NHOLX] = Attributes::makeReal
                     ("NHOLX", "Number of holes in x");
    itsAttr[NHOLY] = Attributes::makeReal
                     ("NHOLY", "Number of holes in y");
    itsAttr[R] = Attributes::makeReal
                 ("R", "Radios of a holes in m");

    registerStringAttribute("OUTFN");
    registerRealAttribute("XSIZE");
    registerRealAttribute("YSIZE");
    registerRealAttribute("PITCH");
    registerRealAttribute("R");
    registerRealAttribute("NHOLX");
    registerRealAttribute("NHOLY");

    registerOwnership();

    setElement((new CollimatorRep("PEPPERPOT"))->makeAlignWrapper());
}


OpalPepperPot::OpalPepperPot(const std::string &name, OpalPepperPot *parent):
    OpalElement(name, parent),
    parmatint_m(NULL) {
    setElement((new CollimatorRep(name))->makeAlignWrapper());
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


    const CollimatorRep *ppo =
        dynamic_cast<const CollimatorRep *>(base.removeWrappers());
    attributeRegistry["XSIZE"]->setReal(ppo->getXsize());
    attributeRegistry["YSIZE"]->setReal(ppo->getYsize());

}

void OpalPepperPot::update() {
    OpalElement::update();

    CollimatorRep *ppo =
        dynamic_cast<CollimatorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    ppo->setElementLength(length);
    ppo->setOutputFN(Attributes::getString(itsAttr[OUTFN]));
    ppo->setXsize(Attributes::getReal(itsAttr[XSIZE]));
    ppo->setYsize(Attributes::getReal(itsAttr[YSIZE]));

    ppo->setRHole(Attributes::getReal(itsAttr[R]));
    ppo->setPitch(Attributes::getReal(itsAttr[PITCH]));
    ppo->setNHoles(Attributes::getReal(itsAttr[NHOLX]), Attributes::getReal(itsAttr[NHOLY]));

    ppo->setPepperPot();

    if(itsAttr[PARTICLEMATERINTERACTION] && parmatint_m == NULL) {
        parmatint_m = (ParticleMatterInteraction::find(Attributes::getString(itsAttr[PARTICLEMATERINTERACTION])))->clone(getOpalName() + std::string("_parmatint"));
        parmatint_m->initParticleMatterInteractionHandler(*ppo);
        ppo->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(ppo);
}