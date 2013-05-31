// ------------------------------------------------------------------------
// $RCSfile: OpalSBend.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalSBend
//   The class of OPAL rectangular bend magnets.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalSBend.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/SBendRep.h"
#include "Fields/BMultipoleField.h"
#include "ComponentWrappers/SBendWrapper.h"
#include "Physics/Physics.h"
#include "Structure/OpalWake.h"
#include "Structure/SurfacePhysics.h"
#include <cmath>

extern Inform *gmsg;

// Class OpalSBend
// ------------------------------------------------------------------------

OpalSBend::OpalSBend():
    OpalBend("SBEND",
             "The \"SBEND\" element defines a sector bending magnet."),
    owk_m(NULL),
    sphys_m(NULL) {
    setElement((new SBendRep("SBEND"))->makeWrappers());
}


OpalSBend::OpalSBend(const std::string &name, OpalSBend *parent):
    OpalBend(name, parent),
    owk_m(NULL),
    sphys_m(NULL) {
    setElement((new SBendRep(name))->makeWrappers());
}


OpalSBend::~OpalSBend() {
    if(owk_m)
        delete owk_m;
    if(sphys_m)
        delete sphys_m;
}


OpalSBend *OpalSBend::clone(const std::string &name) {
    return new OpalSBend(name, this);
}


void OpalSBend::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    // Get the desired field.
    const SBendWrapper *bend =
        dynamic_cast<const SBendWrapper *>(base.removeAlignWrapper());
    BMultipoleField field;

    // Get the desired field.
    if(flag == ERROR_FLAG) {
        field = bend->errorField();
    } else if(flag == ACTUAL_FLAG) {
        field = bend->getField();
    } else if(flag == IDEAL_FLAG) {
        field = bend->getDesign().getField();
    }

    double length = getLength();
    double scale = Physics::c / OpalData::getInstance()->getP0();
    if(length != 0.0) scale *= length;

    for(int i = 1; i <= field.order(); ++i) {
        std::string normName("K0L");
        normName[1] += (i - 1);
        attributeRegistry[normName]->setReal(scale * field.normal(i));

        std::string skewName("K0SL");
        skewName[1] += (i - 1);
        attributeRegistry[skewName]->setReal(scale * field.skew(i));
        scale *= double(i);
    }

    // Store pole face information.
    attributeRegistry["E1"]->setReal(bend->getEntryFaceRotation());
    attributeRegistry["E2"]->setReal(bend->getExitFaceRotation());
    attributeRegistry["H1"]->setReal(bend->getEntryFaceCurvature());
    attributeRegistry["H2"]->setReal(bend->getExitFaceCurvature());

    // Store integration parameters.
    attributeRegistry["SLICES"]->setReal(bend->getSlices());
    attributeRegistry["STEPSIZE"]->setReal(bend->getStepsize());
    //attributeRegistry["FMAPFN"]->setString(bend->getFieldMapFN());
}


void OpalSBend::update() {

    // Define geometry.
    SBendRep *bend = dynamic_cast<SBendRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double angle  = Attributes::getReal(itsAttr[ANGLE]);
    PlanarArcGeometry &geometry = bend->getGeometry();

    if(length) {
        geometry = PlanarArcGeometry(length, angle / length);
    } else {
        geometry = PlanarArcGeometry(angle);
    }

    // Define pole face angles.
    bend->setEntryFaceRotation(Attributes::getReal(itsAttr[E1]));
    bend->setExitFaceRotation(Attributes::getReal(itsAttr[E2]));
    bend->setEntryFaceCurvature(Attributes::getReal(itsAttr[H1]));
    bend->setExitFaceCurvature(Attributes::getReal(itsAttr[H2]));

    // Define integration parameters.
    bend->setSlices(Attributes::getReal(itsAttr[SLICES]));
    bend->setStepsize(Attributes::getReal(itsAttr[STEPSIZE]));

    // Define field.
    double factor = OpalData::getInstance()->getP0() / Physics::c;
    BMultipoleField field;
    double k0 = itsAttr[K0] ? Attributes::getReal(itsAttr[K0]) :
                length ? 2 * sin(angle / 2) / length : angle;
    double k0s = itsAttr[K0S] ? Attributes::getReal(itsAttr[K0S]) : 0.0;
    //JMJ 4/10/2000: above line replaced
    //    length ? angle / length : angle;
    // to avoid closed orbit created by SBEND with defalt K0.
    field.setNormalComponent(1, factor * k0);
    field.setSkewComponent(1, factor * Attributes::getReal(itsAttr[K0S]));
    field.setNormalComponent(2, factor * Attributes::getReal(itsAttr[K1]));
    field.setSkewComponent(2, factor * Attributes::getReal(itsAttr[K1S]));
    field.setNormalComponent(3, factor * Attributes::getReal(itsAttr[K2]) / 2.0);
    field.setSkewComponent(3, factor * Attributes::getReal(itsAttr[K2S]) / 2.0);
    field.setNormalComponent(4, factor * Attributes::getReal(itsAttr[K3]) / 6.0);
    field.setSkewComponent(4, factor * Attributes::getReal(itsAttr[K3S]) / 6.0);
    bend->setField(field);

    // Set field amplitude or bend angle and the magnet rotation about the z axis.
    if(itsAttr[ANGLE])
        bend->setBendAngle(Attributes::getReal(itsAttr[ANGLE]));
    else
        bend->setAmplitudem(sqrt(k0 * k0 + k0s * k0s));

    if(itsAttr[ROTATION])
        bend->setLongitudinalRotation(Attributes::getReal(itsAttr[ROTATION]));
    else if(itsAttr[K0] || itsAttr[K0S])
        bend->setLongitudinalRotation(k0, k0s);
    else
        bend->setLongitudinalRotation(0.0);

    if(itsAttr[FMAPFN])
        bend->setFieldMapFN(Attributes::getString(itsAttr[FMAPFN]));
    else if(bend->getName() != "SBEND")
        ERRORMSG(bend->getName() << ": No filename for a field map given" << endl);

    if(itsAttr[E1])
        bend->setAlpha(Attributes::getReal(itsAttr[E1]));
    else if(itsAttr[ALPHA])
        bend->setAlpha(Attributes::getReal(itsAttr[ALPHA]));
    else
        bend->setAlpha(0.0);

    if(itsAttr[BETA])
        bend->setBeta(Attributes::getReal(itsAttr[BETA]));
    else
        bend->setBeta(0.0);

    if(itsAttr[DESIGNENERGY]) {
        bend->setDesignEnergy(Attributes::getReal(itsAttr[DESIGNENERGY]) * 1e6); // convert from MeV to eV
    }

    if(itsAttr[GAP])
        bend->setFullGap(Attributes::getReal(itsAttr[GAP]));
    else
        bend->setFullGap(0.0);

    if(itsAttr[LENGTH])
        bend->setLength(Attributes::getReal(itsAttr[LENGTH]));
    else
        bend->setLength(0.0);

    if(itsAttr[E2])
        bend->setExitAngle(Attributes::getReal(itsAttr[E2]));
    else if(itsAttr[EXITANGLE])
        bend->setExitAngle(Attributes::getReal(itsAttr[EXITANGLE]));
    //        bend->setExitFaceSlope(tan(Physics::pi * Attributes::getReal(itsAttr[EXITANGLE]) / 180.0));
    else
        bend->setExitAngle(0.0);

    if(itsAttr[WAKEF] && itsAttr[DESIGNENERGY] && owk_m == NULL) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + std::string("_wake"));
        owk_m->initWakefunction(*bend);
        bend->setWake(owk_m->wf_m);
    }

    if(itsAttr[K1])
        bend->setK1(Attributes::getReal(itsAttr[K1]));
    else
        bend->setK1(0.0);
    //  if (itsAttr[L])
    //   bend->setL(Attributes::getReal(itsAttr[L]));
    //  else
    //    bend->setL(0.0);
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
        sphys_m = (SurfacePhysics::find(Attributes::getString(itsAttr[SURFACEPHYSICS])))->clone(getOpalName() + std::string("_sphys"));
        sphys_m->initSurfacePhysicsHandler(*bend, apert_major, apert_minor);
        bend->setSurfacePhysics(sphys_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(bend);
}
