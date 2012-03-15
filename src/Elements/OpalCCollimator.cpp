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
// $Author: Andreas Adelmann $
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
                "The \"CCOLLIMATOR\" element defines a cyclotron collimator."),
    sphys_m(NULL) {
    itsAttr[ANGSTART] = Attributes::makeReal
                        ("ANGSTART", "Start of angle in rad");
    itsAttr[ANGEND] = Attributes::makeReal
                      ("ANGEND", "End of angle in rad");
    itsAttr[RSTART] = Attributes::makeReal
                      ("RSTART", "Start of radius in mm");
    itsAttr[REND] = Attributes::makeReal
                    ("REND", "End of radius in mm");
    itsAttr[WIDTH] = Attributes::makeReal
                     ("WIDTH", "Not used now");
    itsAttr[OUTFN] = Attributes::makeString
                     ("OUTFN", "CCollimator output filename");

    registerRealAttribute("ANGSTART");
    registerRealAttribute("ANGEND");
    registerRealAttribute("RSTART");
    registerRealAttribute("REND");
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
    double angstart = Attributes::getReal(itsAttr[ANGSTART]);
    double angend = Attributes::getReal(itsAttr[ANGEND]);
    double rstart = Attributes::getReal(itsAttr[RSTART]);
    double rend = Attributes::getReal(itsAttr[REND]);

    if (angstart > angend){
      double temp = angstart;
      angstart=angend;
      angend=temp;
      
      temp = rstart;
      rstart=rend;
      rend=temp;
    }
    if (angstart < 0.0){
      angstart = angstart + 2.0*pi;
      angend   =  angend   + 2.0*pi;
    }
    rstart = fabs(rstart);
    rend = fabs(rend);

    coll->setElementLength(length);
    coll->setAngStart(angstart);
    coll->setAngEnd(angend);
    coll->setRStart(rstart);
    coll->setREnd(rend);
    coll->setWidth(Attributes::getReal(itsAttr[WIDTH]));
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
