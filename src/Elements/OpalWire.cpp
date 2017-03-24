// ------------------------------------------------------------------------
// $RCSfile: OpalWire.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalWire
//   The class of OPAL wire collimators.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalWire.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/CollimatorRep.h"
#include "Structure/SurfacePhysics.h"

// Class OpalWire
// ------------------------------------------------------------------------

OpalWire::OpalWire():
    OpalElement(SIZE, "WIRE",
                "The \"WIRE\" element defines a wire."),
    sphys_m(NULL) {
    itsAttr[XSIZE] = Attributes::makeReal
                     ("XSIZE", "Horizontal half-aperture in m");
    itsAttr[YSIZE] = Attributes::makeReal
                     ("YSIZE", "Vertical half-aperture in m");
    itsAttr[XPOS] = Attributes::makeReal
                    ("XPOS", "Horizontal position in m");
    itsAttr[YPOS] = Attributes::makeReal
                    ("YPOS", "Vertical position in m");
    itsAttr[OUTFN] = Attributes::makeString
                     ("OUTFN", "Wire output filename");

    itsAttr[DX] = Attributes::makeReal
      ("DX", "Misalignment in x direction",0.0);
    itsAttr[DY] = Attributes::makeReal
      ("DY", "Misalignment in y direction",0.0);

    registerStringAttribute("OUTFN");
    registerRealAttribute("XSIZE");
    registerRealAttribute("YSIZE");
    registerRealAttribute("XPOS");
    registerRealAttribute("YPOS");
    registerRealAttribute("DX");
    registerRealAttribute("DY");

    registerOwnership();

    setElement((new CollimatorRep("WIRE"))->makeAlignWrapper());
}


OpalWire::OpalWire(const std::string &name, OpalWire *parent):
    OpalElement(name, parent),
    sphys_m(NULL) {
    setElement((new CollimatorRep(name))->makeAlignWrapper());
}


OpalWire::~OpalWire() {
    if(sphys_m)
        delete sphys_m;
}


OpalWire *OpalWire::clone(const std::string &name) {
    return new OpalWire(name, this);
}


void OpalWire::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    const CollimatorRep *coll =
        dynamic_cast<const CollimatorRep *>(base.removeWrappers());
    attributeRegistry["XSIZE"]->setReal(coll->getXsize());
    attributeRegistry["YSIZE"]->setReal(coll->getYsize());

    double dx, dy, dz;
    coll->getMisalignment(dx, dy, dz);
    attributeRegistry["DX"]->setReal(dx);
    attributeRegistry["DY"]->setReal(dy);

    //  attributeRegistry["XPOS"]->setReal(coll->getXpos());
    // attributeRegistry["YPOS"]->setReal(coll->getYpos());
}


void OpalWire::update() {
    OpalElement::update();

    double dx = Attributes::getReal(itsAttr[DX]);
    double dy = Attributes::getReal(itsAttr[DY]);

    CollimatorRep *coll =
        dynamic_cast<CollimatorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    coll->setElementLength(length);
    coll->setXsize(Attributes::getReal(itsAttr[XSIZE]));
    coll->setYsize(Attributes::getReal(itsAttr[YSIZE]));
    coll->setXpos(Attributes::getReal(itsAttr[XPOS]));
    coll->setYpos(Attributes::getReal(itsAttr[YPOS]));
    coll->setOutputFN(Attributes::getString(itsAttr[OUTFN]));
    coll->setMisalignment(dx, dy, 0.0);
    coll->setWire();

    /*
    std::vector<double> apert = getApert();
    double apert_major = -1., apert_minor = -1.;
    if(apert.size() > 0) {
        apert_major = apert[0];
        if(apert.size() > 1) {
            apert_minor = apert[1];
        } else {
            apert_minor = apert[0];
        }
    }
    */
    if(itsAttr[SURFACEPHYSICS] && sphys_m == NULL) {
        sphys_m = (SurfacePhysics::find(Attributes::getString(itsAttr[SURFACEPHYSICS])))->clone(getOpalName() + std::string("_sphys"));
        sphys_m->initSurfacePhysicsHandler(*coll);
        coll->setSurfacePhysics(sphys_m->handler_m);
    }


    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(coll);
}