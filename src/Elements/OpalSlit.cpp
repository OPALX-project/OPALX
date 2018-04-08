// ------------------------------------------------------------------------
// $RCSfile: OpalCollimator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalSlit
//   The class of OPAL elliptic collimators.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalSlit.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/FlexibleCollimatorRep.h"
#include "Structure/ParticleMatterInteraction.h"

#include <boost/regex.hpp>
#include <cstdlib>

// Class OpalSlit
// ------------------------------------------------------------------------

OpalSlit::OpalSlit():
    OpalElement(SIZE, "SLIT",
                "The \"SLIT\" element defines a slit."),
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

    setElement((new FlexibleCollimatorRep("SLIT"))->makeAlignWrapper());
}


OpalSlit::OpalSlit(const std::string &name, OpalSlit *parent):
    OpalElement(name, parent),
    parmatint_m(NULL) {
    setElement((new FlexibleCollimatorRep(name))->makeAlignWrapper());
}


OpalSlit::~OpalSlit() {
    if(parmatint_m)
        delete parmatint_m;
}


OpalSlit *OpalSlit::clone(const std::string &name) {
    return new OpalSlit(name, this);
}


void OpalSlit::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
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


void OpalSlit::update() {
    OpalElement::update();

    FlexibleCollimatorRep *coll =
        dynamic_cast<FlexibleCollimatorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    coll->setElementLength(length);

    if (getOpalName() != "SLIT") {
        double width = 2 * Attributes::getReal(itsAttr[XSIZE]);
        double height = 2 * Attributes::getReal(itsAttr[YSIZE]);
        std::stringstream description;
        description << "rectangle(" << width << "," << height << ")";
        coll->setDescription(description.str());
    }

    coll->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

    if(itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == NULL) {
        parmatint_m = (ParticleMatterInteraction::find(Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION])))->clone(getOpalName() + std::string("_parmatint"));
        parmatint_m->initParticleMatterInteractionHandler(*coll);
        coll->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    std::vector<double> apert = {Attributes::getReal(itsAttr[XSIZE]),
                                 Attributes::getReal(itsAttr[YSIZE]),
                                 1.0};
    coll->setAperture(ElementBase::RECTANGULAR, apert );

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(coll);
}