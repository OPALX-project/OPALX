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
#include "BeamlineCore/CollimatorRep.h"
#include "Structure/SurfacePhysics.h"


// Class OpalECollimator
// ------------------------------------------------------------------------

OpalECollimator::OpalECollimator():
    OpalElement(SIZE, "ECOLLIMATOR",
                "The \"ECOLLIMATOR\" element defines an elliptic collimator."),
    sphys_m(NULL) {
    itsAttr[XSIZE] = Attributes::makeReal
                     ("XSIZE", "Horizontal half-aperture in m");
    itsAttr[YSIZE] = Attributes::makeReal
                     ("YSIZE", "Vertical half-aperture in m");
    itsAttr[OUTFN] = Attributes::makeString
                     ("OUTFN", "Monitor output filename");
    itsAttr[DX] = Attributes::makeReal
      ("DX", "Misalignment in x direction",0.0);
    itsAttr[DY] = Attributes::makeReal
      ("DY", "Misalignment in y direction",0.0);
    itsAttr[DZ] = Attributes::makeReal
      ("DZ", "Misalignment in z direction",0.0);

    registerStringAttribute("OUTFN");
    registerRealAttribute("XSIZE");
    registerRealAttribute("YSIZE");
    registerRealAttribute("DX");
    registerRealAttribute("DY");
    registerRealAttribute("DZ");
    setElement((new CollimatorRep("ECOLLIMATOR"))->makeAlignWrapper());
}


OpalECollimator::OpalECollimator(const string &name, OpalECollimator *parent):
    OpalElement(name, parent),
    sphys_m(NULL) {
    setElement((new CollimatorRep(name))->makeAlignWrapper());
}


OpalECollimator::~OpalECollimator() {
    if(sphys_m)
        delete sphys_m;
}


OpalECollimator *OpalECollimator::clone(const string &name) {
    return new OpalECollimator(name, this);
}


void OpalECollimator::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    const CollimatorRep *coll =
        dynamic_cast<const CollimatorRep *>(base.removeWrappers());
    attributeRegistry["XSIZE"]->setReal(coll->getXsize());
    attributeRegistry["YSIZE"]->setReal(coll->getYsize());
    double dx, dy, dz;
    coll->getMisalignment(dx, dy, dz);
    attributeRegistry["DX"]->setReal(dx);
    attributeRegistry["DY"]->setReal(dy);
    attributeRegistry["DZ"]->setReal(dz);
}


void OpalECollimator::update() {

    double dx = Attributes::getReal(itsAttr[DX]);
    double dy = Attributes::getReal(itsAttr[DY]);
    double dz = Attributes::getReal(itsAttr[DZ]);

    CollimatorRep *coll =
        dynamic_cast<CollimatorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    coll->setElementLength(length);
    coll->setXsize(Attributes::getReal(itsAttr[XSIZE]));
    coll->setYsize(Attributes::getReal(itsAttr[YSIZE]));
    coll->setOutputFN(Attributes::getString(itsAttr[OUTFN]));
    coll->setMisalignment(dx, dy, dz);
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

    if(itsAttr[SURFACEPHYSICS] && sphys_m == NULL) {
        sphys_m = (SurfacePhysics::find(Attributes::getString(itsAttr[SURFACEPHYSICS])))->clone(getOpalName() + string("_sphys"));
        sphys_m->initSurfacePhysicsHandler(*coll, apert_major, apert_minor);
        coll->setSurfacePhysics(sphys_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(coll);
}
