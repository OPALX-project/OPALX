// ------------------------------------------------------------------------
// $RCSfile: OpalPepperPot.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalPepperPot
//   The class of OPAL elliptic collimators.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalPepperPot.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/CollimatorRep.h"
#include "Structure/SurfacePhysics.h"


// Class OpalPepperPot
// ------------------------------------------------------------------------

OpalPepperPot::OpalPepperPot():
    OpalElement(SIZE, "PEPPERPOT",
                "The \"PEPPERPOT\" element defines an elliptic collimator."),
    sphys_m(NULL) {
    itsAttr[XSIZE] = Attributes::makeReal
                     ("XSIZE", "Size in x of the pepperpot in m");
    itsAttr[YSIZE] = Attributes::makeReal
                     ("YSIZE", "Size in y of the pepperpot in m");
    itsAttr[OUTFN] = Attributes::makeString
                     ("OUTFN", "Pepperpot output filename");
    itsAttr[PITCH] = Attributes::makeReal
                     ("PITCH", "Pitch of the pepperpot in m");
    itsAttr[NHOLX] = Attributes::makeReal
                     ("NHOLX", "Number of holes in x");
    itsAttr[NHOLY] = Attributes::makeReal
                     ("NHOLY", "Number of holes in y");
    itsAttr[R] = Attributes::makeReal
                 ("R", "Radios of a holes in m");

    itsAttr[DX] = Attributes::makeReal
      ("DX", "Misalignment in x direction",0.0);
    itsAttr[DY] = Attributes::makeReal
      ("DY", "Misalignment in y direction",0.0);

    registerStringAttribute("OUTFN");
    registerRealAttribute("XSIZE");
    registerRealAttribute("YSIZE");
    registerRealAttribute("PITCH");
    registerRealAttribute("R");
    registerRealAttribute("NHOLX");
    registerRealAttribute("NHOLY");
    registerRealAttribute("DX");
    registerRealAttribute("DY");

    setElement((new CollimatorRep("PEPPERPOT"))->makeAlignWrapper());
}


OpalPepperPot::OpalPepperPot(const string &name, OpalPepperPot *parent):
    OpalElement(name, parent),
    sphys_m(NULL) {
    setElement((new CollimatorRep(name))->makeAlignWrapper());
}


OpalPepperPot::~OpalPepperPot() {
    if(sphys_m)
        delete sphys_m;
}


OpalPepperPot *OpalPepperPot::clone(const string &name) {
    return new OpalPepperPot(name, this);
}


void OpalPepperPot::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);


    const CollimatorRep *ppo =
        dynamic_cast<const CollimatorRep *>(base.removeWrappers());
    attributeRegistry["XSIZE"]->setReal(ppo->getXsize());
    attributeRegistry["YSIZE"]->setReal(ppo->getYsize());
    double dx, dy, dz;
    ppo->getMisalignment(dx, dy, dz);
    attributeRegistry["DX"]->setReal(dx);
    attributeRegistry["DY"]->setReal(dy);
}

void OpalPepperPot::update() {
    CollimatorRep *ppo =
        dynamic_cast<CollimatorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    ppo->setElementLength(length);
    ppo->setOutputFN(Attributes::getString(itsAttr[OUTFN]));
    ppo->setXsize(Attributes::getReal(itsAttr[XSIZE]));
    ppo->setYsize(Attributes::getReal(itsAttr[YSIZE]));

    ppo->setRHole(Attributes::getReal(itsAttr[R]));
    ppo->setPitch(Attributes::getReal(itsAttr[PITCH]));
    ppo->setNHoles(Attributes::getReal(itsAttr[NHOLX]), Attributes::getReal(itsAttr[NHOLY]));

    double dx = Attributes::getReal(itsAttr[DX]);
    double dy = Attributes::getReal(itsAttr[DY]);

    ppo->setMisalignment(dx, dy, 0.0);

    ppo->setPepperPot();
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
        sphys_m = (SurfacePhysics::find(Attributes::getString(itsAttr[SURFACEPHYSICS])))->clone(getOpalName() + string("_sphys"));
        sphys_m->initSurfacePhysicsHandler(*ppo);
        ppo->setSurfacePhysics(sphys_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(ppo);
}
