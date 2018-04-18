// ------------------------------------------------------------------------
// $RCSfile: OpalCollimator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalECollimator
//   The class of OPAL elliptic collimators.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalECollimator.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/FlexibleCollimatorRep.h"
#include "Structure/ParticleMatterInteraction.h"

#include <boost/regex.hpp>
#include <cstdlib>

// Class OpalECollimator
// ------------------------------------------------------------------------

OpalECollimator::OpalECollimator():
    OpalElement(SIZE, "ECOLLIMATOR",
                "The \"ECOLLIMATOR\" element defines an elliptic collimator."),
    parmatint_m(NULL) {
    itsAttr[XSIZE] = Attributes::makeReal
                     ("XSIZE", "Horizontal half-aperture in m");
    itsAttr[YSIZE] = Attributes::makeReal
                     ("YSIZE", "Vertical half-aperture in m");
    itsAttr[OUTFN] = Attributes::makeString
                     ("OUTFN", "Monitor output filename");


    registerStringAttribute("OUTFN");
    registerRealAttribute("XSIZE");
    registerRealAttribute("YSIZE");

    registerOwnership();

    setElement((new FlexibleCollimatorRep("ECOLLIMATOR"))->makeAlignWrapper());
}


OpalECollimator::OpalECollimator(const std::string &name, OpalECollimator *parent):
    OpalElement(name, parent),
    parmatint_m(NULL) {
    setElement((new FlexibleCollimatorRep(name))->makeAlignWrapper());
}


OpalECollimator::~OpalECollimator() {
    if(parmatint_m)
        delete parmatint_m;
}


OpalECollimator *OpalECollimator::clone(const std::string &name) {
    return new OpalECollimator(name, this);
}


void OpalECollimator::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    const FlexibleCollimatorRep *coll =
        dynamic_cast<const FlexibleCollimatorRep *>(base.removeWrappers());
    std::string Double("(-?[0-9]+\\.?[0-9]*([Ee][+-]?[0-9]+)?)");
    std::string desc = coll->getDescription();

    boost::regex parser("ellipse\\(" + Double + "," + Double + "\\)");
    boost::smatch what;
    if (!boost::regex_match(desc, what, parser)) return;

    double width = atof(std::string(what[1]).c_str());
    double height = atof(std::string(what[3]).c_str());

    attributeRegistry["XSIZE"]->setReal(0.5 * width);
    attributeRegistry["YSIZE"]->setReal(0.5 * height);
}


void OpalECollimator::update() {
    OpalElement::update();

    FlexibleCollimatorRep *coll =
        dynamic_cast<FlexibleCollimatorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    coll->setElementLength(length);

    if (getOpalName() != "ECOLLIMATOR") {
        double width = 2 * Attributes::getReal(itsAttr[XSIZE]);
        double height = 2 * Attributes::getReal(itsAttr[YSIZE]);
        std::stringstream description;
        description << "ellipse(" << width << "," << height << ")";
        coll->setDescription(description.str());
    }

    coll->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

    if(itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == NULL) {
        const std::string matterDescriptor = Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION]);
        ParticleMatterInteraction *orig = ParticleMatterInteraction::find(matterDescriptor);
        parmatint_m = orig->clone(getOpalName() + std::string("_parmatint"));
        parmatint_m->initParticleMatterInteractionHandler(*coll);
        coll->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(coll);
}