// ------------------------------------------------------------------------
// $RCSfile: OpalRCollimator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalRCollimator
//   The class of OPAL rectangular collimators.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalRCollimator.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/FlexibleCollimatorRep.h"
#include "Structure/ParticleMatterInteraction.h"

#include <boost/regex.hpp>
#include <cstdlib>

// Class OpalRCollimator
// ------------------------------------------------------------------------

OpalRCollimator::OpalRCollimator():
    OpalElement(SIZE, "RCOLLIMATOR",
                "The \"RCOLLIMATOR\" element defines a rectangular collimator."),
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

    setElement((new FlexibleCollimatorRep("RCOLLIMATOR"))->makeAlignWrapper());
}


OpalRCollimator::OpalRCollimator(const std::string &name, OpalRCollimator *parent):
    OpalElement(name, parent),
    parmatint_m(NULL) {
    setElement((new FlexibleCollimatorRep(name))->makeAlignWrapper());
}


OpalRCollimator::~OpalRCollimator() {
    if(parmatint_m)
        delete parmatint_m;
}


OpalRCollimator *OpalRCollimator::clone(const std::string &name) {
    return new OpalRCollimator(name, this);
}


void OpalRCollimator::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    const FlexibleCollimatorRep *coll =
        dynamic_cast<const FlexibleCollimatorRep *>(base.removeWrappers());

    std::string Double("(-?[0-9]+\\.?[0-9]*([Ee][+-]?[0-9]+)?)");
    std::string desc = coll->getDescription();

    boost::regex parser("rectangle\\(" + Double + "," + Double + "\\)");
    boost::smatch what;
    if (!boost::regex_match(desc, what, parser)) return;

    double width = atof(std::string(what[1]).c_str());
    double height = atof(std::string(what[3]).c_str());

    attributeRegistry["XSIZE"]->setReal(0.5 * width);
    attributeRegistry["YSIZE"]->setReal(0.5 * height);
}


void OpalRCollimator::update() {
    OpalElement::update();

    FlexibleCollimatorRep *coll =
        dynamic_cast<FlexibleCollimatorRep *>(getElement()->removeWrappers());
    coll->setElementLength(Attributes::getReal(itsAttr[LENGTH]));

    if (getOpalName() != "RCOLLIMATOR") {
        double width = 2 * Attributes::getReal(itsAttr[XSIZE]);
        double height = 2 * Attributes::getReal(itsAttr[YSIZE]);
        std::stringstream description;
        description << "rectangle(" << width << "," << height << ")";
        coll->setDescription(description.str());

        coll->setOutputFN(Attributes::getString(itsAttr[OUTFN]));
    }

    if(itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == NULL) {
        parmatint_m = (ParticleMatterInteraction::find(Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION])))->clone(getOpalName() + std::string("_parmatint"));
        parmatint_m->initParticleMatterInteractionHandler(*coll);
        coll->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(coll);
}