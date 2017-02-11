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
#include "BeamlineCore/CollimatorRep.h"
#include "Structure/SurfacePhysics.h"

// Class OpalSlit
// ------------------------------------------------------------------------

OpalSlit::OpalSlit():
    OpalElement(SIZE, "SLIT",
                "The \"SLIT\" element defines a slit."),
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

    setElement((new CollimatorRep("SLIT"))->makeAlignWrapper());
}


OpalSlit::OpalSlit(const std::string &name, OpalSlit *parent):
    OpalElement(name, parent),
    sphys_m(NULL) {
    setElement((new CollimatorRep(name))->makeAlignWrapper());
}


OpalSlit::~OpalSlit() {
    if(sphys_m)
        delete sphys_m;
}


OpalSlit *OpalSlit::clone(const std::string &name) {
    return new OpalSlit(name, this);
}


void OpalSlit::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    const CollimatorRep *coll =
        dynamic_cast<const CollimatorRep *>(base.removeWrappers());
    attributeRegistry["XSIZE"]->setReal(coll->getXsize());
    attributeRegistry["YSIZE"]->setReal(coll->getYsize());
}


void OpalSlit::update() {
    OpalElement::update();

    CollimatorRep *coll =
        dynamic_cast<CollimatorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    coll->setElementLength(length);
    coll->setXsize(Attributes::getReal(itsAttr[XSIZE]));
    coll->setYsize(Attributes::getReal(itsAttr[YSIZE]));
    coll->setOutputFN(Attributes::getString(itsAttr[OUTFN]));
    coll->setSlit();

    if(itsAttr[SURFACEPHYSICS] && sphys_m == NULL) {
        sphys_m = (SurfacePhysics::find(Attributes::getString(itsAttr[SURFACEPHYSICS])))->clone(getOpalName() + std::string("_sphys"));
        sphys_m->initSurfacePhysicsHandler(*coll);
        coll->setSurfacePhysics(sphys_m->handler_m);
    }

    std::vector<double> apert = {Attributes::getReal(itsAttr[XSIZE]),
                                 Attributes::getReal(itsAttr[YSIZE])};
    coll->setAperture(ElementBase::RECTANGULAR, apert );

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(coll);
}
