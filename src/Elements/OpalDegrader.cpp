// ------------------------------------------------------------------------
// $RCSfile: OpalDegrader.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalDegrader
//   The class of OPAL Degrader.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalDegrader.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/DegraderRep.h"
#include "Structure/SurfacePhysics.h"


// Class OpalDegrader
// ------------------------------------------------------------------------

OpalDegrader::OpalDegrader():
    OpalElement(SIZE, "DEGRADER",
                "The \"DEGRADER\" element defines a degrader."),
    sphys_m(NULL) {
    itsAttr[XSIZE] = Attributes::makeReal
        ("XSIZE", "not used",0.0);
    itsAttr[YSIZE] = Attributes::makeReal
        ("YSIZE", "not used",0.0);
    itsAttr[OUTFN] = Attributes::makeString
        ("OUTFN", "Degrader output filename");
    itsAttr[DX] = Attributes::makeReal
        ("DX", "not used",0.0);
    itsAttr[DY] = Attributes::makeReal
        ("DY", "not used",0.0);

    registerStringAttribute("OUTFN");
    registerRealAttribute("DX");
    registerRealAttribute("DY");

    registerOwnership();

    setElement((new DegraderRep("DEGRADER"))->makeAlignWrapper());
}


OpalDegrader::OpalDegrader(const std::string &name, OpalDegrader *parent):
    OpalElement(name, parent),
    sphys_m(NULL) {
    setElement((new DegraderRep(name))->makeAlignWrapper());
}


OpalDegrader::~OpalDegrader() {
    if(sphys_m)
        delete sphys_m;
}


OpalDegrader *OpalDegrader::clone(const std::string &name) {
    return new OpalDegrader(name, this);
}


void OpalDegrader::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    const DegraderRep *deg =
        dynamic_cast<const DegraderRep *>(base.removeWrappers());

    double dx, dy, dz;
    deg->getMisalignment(dx, dy, dz);
}


void OpalDegrader::update() {
    OpalElement::update();

    DegraderRep *deg =
        dynamic_cast<DegraderRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    deg->setElementLength(length);
    deg->setZSize(length);

    deg->setOutputFN(Attributes::getString(itsAttr[OUTFN]));
    deg->setMisalignment(0.0, 0.0, 0.0);

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
        sphys_m->initSurfacePhysicsHandler(*deg);
        deg->setSurfacePhysics(sphys_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(deg);
}