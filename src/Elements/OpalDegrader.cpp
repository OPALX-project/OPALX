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
    setElement((new DegraderRep("DEGRADER"))->makeAlignWrapper());
}


OpalDegrader::OpalDegrader(const string &name, OpalDegrader *parent):
    OpalElement(name, parent),
    sphys_m(NULL) {
    setElement((new DegraderRep(name))->makeAlignWrapper());
}


OpalDegrader::~OpalDegrader() {
    if(sphys_m)
        delete sphys_m;
}


OpalDegrader *OpalDegrader::clone(const string &name) {
    return new OpalDegrader(name, this);
}


void OpalDegrader::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    const DegraderRep *deg =
        dynamic_cast<const DegraderRep *>(base.removeWrappers());

    double dx, dy, dz;
    deg->getMisalignment(dx, dy, dz);
    attributeRegistry["DZ"]->setReal(dz);
}


void OpalDegrader::update() {

    double dx = Attributes::getReal(itsAttr[DX]);
    double dy = Attributes::getReal(itsAttr[DY]);
    double dz = Attributes::getReal(itsAttr[DZ]);

    DegraderRep *deg =
        dynamic_cast<DegraderRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    deg->setElementLength(length);
    deg->setOutputFN(Attributes::getString(itsAttr[OUTFN]));
    deg->setMisalignment(dx, dy, dz);
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
        sphys_m->initSurfacePhysicsHandler(*deg, apert_major, apert_minor);
        deg->setSurfacePhysics(sphys_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(deg);
}
