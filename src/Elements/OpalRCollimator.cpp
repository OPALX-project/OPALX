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
#include "BeamlineCore/CollimatorRep.h"
#include "Structure/SurfacePhysics.h"

// Class OpalRCollimator
// ------------------------------------------------------------------------

OpalRCollimator::OpalRCollimator():
    OpalElement(SIZE, "RCOLLIMATOR",
                "The \"RCOLLIMATOR\" element defines a rectangular collimator."),
    sphys_m(NULL) {
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

    setElement((new CollimatorRep("RCOLLIMATOR"))->makeAlignWrapper());
}


OpalRCollimator::OpalRCollimator(const std::string &name, OpalRCollimator *parent):
    OpalElement(name, parent),
    sphys_m(NULL) {
    setElement((new CollimatorRep(name))->makeAlignWrapper());
}


OpalRCollimator::~OpalRCollimator() {
    if(sphys_m)
        delete sphys_m;
}


OpalRCollimator *OpalRCollimator::clone(const std::string &name) {
    return new OpalRCollimator(name, this);
}


void OpalRCollimator::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    const CollimatorRep *coll =
        dynamic_cast<const CollimatorRep *>(base.removeWrappers());
    attributeRegistry["XSIZE"]->setReal(coll->getXsize());
    attributeRegistry["YSIZE"]->setReal(coll->getYsize());
}


void OpalRCollimator::update() {
    OpalElement::update();

    CollimatorRep *coll =
        dynamic_cast<CollimatorRep *>(getElement()->removeWrappers());
    coll->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    coll->setXsize(Attributes::getReal(itsAttr[XSIZE]));
    coll->setYsize(Attributes::getReal(itsAttr[YSIZE]));
    coll->setOutputFN(Attributes::getString(itsAttr[OUTFN]));
    coll->setRColl();

    if(itsAttr[SURFACEPHYSICS] && sphys_m == NULL) {
        sphys_m = (SurfacePhysics::find(Attributes::getString(itsAttr[SURFACEPHYSICS])))->clone(getOpalName() + std::string("_sphys"));
        sphys_m->initSurfacePhysicsHandler(*coll);
        coll->setSurfacePhysics(sphys_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(coll);
}