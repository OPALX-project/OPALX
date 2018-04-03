// ------------------------------------------------------------------------
// $RCSfile: OpalCollimator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalFlexibleCollimator
//   The class of OPAL elliptic collimators.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalFlexibleCollimator.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/FlexibleCollimatorRep.h"
#include "Structure/ParticleMatterInteraction.h"

// Class OpalFlexibleCollimator
// ------------------------------------------------------------------------

OpalFlexibleCollimator::OpalFlexibleCollimator():
    OpalElement(SIZE, "FLEXIBLECOLLIMATOR",
                "The \"FLEXIBLECOLLIMATOR\" element defines a slit."),
    partMatInt_m(NULL) {
    itsAttr[FNAME] = Attributes::makeString
                     ("FNAME", "File name containing description of holes");
    itsAttr[OUTFN] = Attributes::makeString
                     ("OUTFN", "File name of log file for deleted particles");


    registerStringAttribute("OUTFN");
    registerStringAttribute("FNAME");

    registerOwnership();

    setElement((new FlexibleCollimatorRep("FLEXIBLECOLLIMATOR"))->makeAlignWrapper());
}


OpalFlexibleCollimator::OpalFlexibleCollimator(const std::string &name, OpalFlexibleCollimator *parent):
    OpalElement(name, parent),
    partMatInt_m(NULL) {
    setElement((new FlexibleCollimatorRep(name))->makeAlignWrapper());
}


OpalFlexibleCollimator::~OpalFlexibleCollimator() {
    if(partMatInt_m)
        delete partMatInt_m;
}


OpalFlexibleCollimator *OpalFlexibleCollimator::clone(const std::string &name) {
    return new OpalFlexibleCollimator(name, this);
}


void OpalFlexibleCollimator::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
}


void OpalFlexibleCollimator::update() {
    OpalElement::update();

    FlexibleCollimatorRep *coll =
        dynamic_cast<FlexibleCollimatorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    coll->setElementLength(length);

    std::string fname = Attributes::getString(itsAttr[FNAME]);
    if (fname != "") {
        std::ifstream it(fname);
        std::string str((std::istreambuf_iterator<char>(it)),
                        std::istreambuf_iterator<char>());
        coll->setDescription(str);
    }
    coll->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

    if(itsAttr[PARTICLEMATTERINTERACTION] && partMatInt_m == NULL) {
        partMatInt_m = (ParticleMatterInteraction::find(Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION])))->clone(getOpalName() + std::string("_parmatint"));
        partMatInt_m->initParticleMatterInteractionHandler(*coll);
        coll->setParticleMatterInteraction(partMatInt_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(coll);
}