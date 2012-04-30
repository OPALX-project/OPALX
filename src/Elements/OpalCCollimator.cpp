// ------------------------------------------------------------------------
// $RCSfile: OpalCCollimator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalCCollimator
//   The class of OPAL Cyclotron collimators.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann, Jianjun Yang $
//
// ------------------------------------------------------------------------

#include "Elements/OpalCCollimator.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/CollimatorRep.h"
#include "Structure/SurfacePhysics.h"
#include "Physics/Physics.h"

using Physics::pi;

// Class OpalCCollimator
// ------------------------------------------------------------------------

OpalCCollimator::OpalCCollimator():
    OpalElement(SIZE, "CCOLLIMATOR",
                "The \"CCOLLIMATOR\" element defines a rectangular-shape cyclotron collimator"),
    sphys_m(NULL) {
    itsAttr[XSTART] = Attributes::makeReal
                      ("XSTART", " Start of x coordinate [mm]");
    itsAttr[XEND] = Attributes::makeReal
                    ("XEND", " End of x coordinate, [mm]");
    itsAttr[YSTART] = Attributes::makeReal
                      ("YSTART", "Start of y coordinate, [mm]");
    itsAttr[YEND] = Attributes::makeReal
                    ("YEND", "End of y coordinate, [mm]");
    itsAttr[ZSTART] = Attributes::makeReal
      ("ZSTART", "Start of vertical coordinate, [mm], default value: -100",-100.0);
    itsAttr[ZEND] = Attributes::makeReal
      ("ZEND", "End of vertical coordinate, [mm], default value: 100", 100.0);
    itsAttr[WIDTH] = Attributes::makeReal
                     ("WIDTH", "Width of the collimator [mm]");
    itsAttr[OUTFN] = Attributes::makeString
                     ("OUTFN", "CCollimator output filename");

    registerRealAttribute("XSTART");
    registerRealAttribute("XEND");
    registerRealAttribute("YSTART");
    registerRealAttribute("YEND");
    registerRealAttribute("ZSTART");
    registerRealAttribute("ZEND");
    registerRealAttribute("WIDTH");
    registerStringAttribute("OUTFN");

    setElement((new CollimatorRep("CCOLLIMATOR"))->makeAlignWrapper());
}


OpalCCollimator::OpalCCollimator(const string &name, OpalCCollimator *parent):
    OpalElement(name, parent),
    sphys_m(NULL) {
    setElement((new CollimatorRep(name))->makeAlignWrapper());
}


OpalCCollimator::~OpalCCollimator() {
    if(sphys_m)
        delete sphys_m;
}


OpalCCollimator *OpalCCollimator::clone(const string &name) {
    return new OpalCCollimator(name, this);
}


void OpalCCollimator::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
}


void OpalCCollimator::update() {
    CollimatorRep *coll =
        dynamic_cast<CollimatorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double xstart = Attributes::getReal(itsAttr[XSTART]);
    double xend = Attributes::getReal(itsAttr[XEND]);
    double ystart = Attributes::getReal(itsAttr[YSTART]);
    double yend = Attributes::getReal(itsAttr[YEND]);
    double zstart = Attributes::getReal(itsAttr[ZSTART]);
    double zend = Attributes::getReal(itsAttr[ZEND]);
    double width = Attributes::getReal(itsAttr[WIDTH]);
    coll->setElementLength(length);
    coll->setXStart(xstart);
    coll->setXEnd(xend);
    coll->setYStart(ystart);
    coll->setYEnd(yend);
    coll->setZStart(zstart);
    coll->setZEnd(zend);
    coll->setWidth(width);
    coll->setOutputFN(Attributes::getString(itsAttr[OUTFN]));
    coll->setCColl();

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
