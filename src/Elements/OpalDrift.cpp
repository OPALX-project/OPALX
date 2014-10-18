// ------------------------------------------------------------------------
// $RCSfile: OpalDrift.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalDrift
//   The class of OPAL drift spaces.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalDrift.h"
#include "Structure/BoundaryGeometry.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/DriftRep.h"
#include "Structure/OpalWake.h"
#include "Structure/SurfacePhysics.h"

// Class OpalDrift
// ------------------------------------------------------------------------

OpalDrift::OpalDrift():
    OpalElement(SIZE, "DRIFT",
                "The \"DRIFT\" element defines a drift space."),
    owk_m(NULL),
    sphys_m(NULL),
    obgeo_m(NULL) {
    // CKR: the following 3 lines are redundant: OpalElement does this already!
    //      they prevent drift from working properly
    //
    //     itsAttr[LENGTH] = Attributes::makeReal
    //         ("LENGTH", "Drift length");

    //     registerRealAttribute("LENGTH");
    itsAttr[GEOMETRY] = Attributes::makeString
                        ("GEOMETRY", "BoundaryGeometry for Drifts");
    registerStringAttribute("GEOMETRY");
    setElement(new DriftRep("DRIFT"));
}


OpalDrift::OpalDrift(const std::string &name, OpalDrift *parent):
    OpalElement(name, parent),
    owk_m(NULL),
    sphys_m(NULL),
    obgeo_m(NULL) {
    setElement(new DriftRep(name));
}


OpalDrift::~OpalDrift() {
    if(owk_m)
        delete owk_m;
    if(sphys_m)
        delete sphys_m;
    if(obgeo_m)
	delete obgeo_m;
}


OpalDrift *OpalDrift::clone(const std::string &name) {
    return new OpalDrift(name, this);
}


bool OpalDrift::isDrift() const {
    return true;
}


void OpalDrift::update() {

    DriftRep *drf = static_cast<DriftRep *>(getElement());
    drf->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    if(itsAttr[WAKEF] && owk_m == NULL) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + std::string("_wake"));
        owk_m->initWakefunction(*drf);
        drf->setWake(owk_m->wf_m);
    }
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
        sphys_m->initSurfacePhysicsHandler(*drf);
        drf->setSurfacePhysics(sphys_m->handler_m);
    }
    if(itsAttr[GEOMETRY] && obgeo_m == NULL) {
        obgeo_m = (BoundaryGeometry::find(Attributes::getString(itsAttr[GEOMETRY])))->clone(getOpalName() + std::string("_geometry"));
        if(obgeo_m) {
	    drf->setBoundaryGeometry(obgeo_m);
        }
    }


    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(drf);
}
