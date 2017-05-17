// ------------------------------------------------------------------------
// $RCSfile: OpalRBend.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalRBend
//   The class of OPAL rectangular bend magnets.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalRBend.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/RBendRep.h"
#include "Fields/BMultipoleField.h"
#include "ComponentWrappers/RBendWrapper.h"
#include "Physics/Physics.h"
#include "Structure/OpalWake.h"
#include "Structure/ParticleMatterInteraction.h"
#include "Utilities/OpalException.h"
#include <cmath>

// Class OpalRBend
// ------------------------------------------------------------------------

OpalRBend::OpalRBend():
    OpalBend("RBEND",
             "The \"RBEND\" element defines a rectangular bending magnet."),
    owk_m(0),
    parmatint_m(NULL) {

    registerOwnership();

    setElement((new RBendRep("RBEND"))->makeWrappers());
}


OpalRBend::OpalRBend(const std::string &name, OpalRBend *parent):
    OpalBend(name, parent),
    owk_m(0),
    parmatint_m(NULL) {
    setElement((new RBendRep(name))->makeWrappers());
}


OpalRBend::~OpalRBend() {
    if(owk_m)
        delete owk_m;
    if(parmatint_m)
        delete parmatint_m;
}


OpalRBend *OpalRBend::clone(const std::string &name) {
    return new OpalRBend(name, this);
}


void OpalRBend::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    // Get the desired field.
    const RBendWrapper *bend =
        dynamic_cast<const RBendWrapper *>(base.removeAlignWrapper());
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
}


void OpalRBend::update() {
    OpalElement::update();

    // Define geometry.
    RBendRep *bend =
        dynamic_cast<RBendRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double angle  = Attributes::getReal(itsAttr[ANGLE]);
    double e1     = Attributes::getReal(itsAttr[E1]);
    RBendGeometry &geometry = bend->getGeometry();
    geometry.setElementLength(length);
    if (angle < 0) {

    }
    geometry.setBendAngle(angle);

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
    double k0 =
        itsAttr[K0] ? Attributes::getReal(itsAttr[K0]) :
        length ? 2 * sin(angle / 2) / length : angle;
    double k0s = itsAttr[K0S] ? Attributes::getReal(itsAttr[K0S]) : 0.0;
    //JMJ 4/10/2000: above line replaced
    //    length ? angle / length : angle;
    // to avoid closed orbit created by RBEND with defalt K0.
    field.setNormalComponent(1, factor * k0);
    field.setSkewComponent(1, factor * Attributes::getReal(itsAttr[K0S]));
    field.setNormalComponent(2, factor * Attributes::getReal(itsAttr[K1]));
    field.setSkewComponent(2, factor * Attributes::getReal(itsAttr[K1S]));
    field.setNormalComponent(3, factor * Attributes::getReal(itsAttr[K2]) / 2.0);
    field.setSkewComponent(3, factor * Attributes::getReal(itsAttr[K2S]) / 2.0);
    field.setNormalComponent(4, factor * Attributes::getReal(itsAttr[K3]) / 6.0);
    field.setSkewComponent(4, factor * Attributes::getReal(itsAttr[K3S]) / 6.0);
    bend->setField(field);

    // Set field amplitude or bend angle.
    if(itsAttr[ANGLE]) {
        if (bend->isPositioned() && angle < 0.0) {
            e1 = -e1;
            angle = -angle;

            Quaternion rotAboutZ(0, 0, 0, 1);
            CoordinateSystemTrafo g2l = bend->getCSTrafoGlobal2Local();
            bend->releasePosition();
            bend->setCSTrafoGlobal2Local(CoordinateSystemTrafo(g2l.getOrigin(),
                                                               rotAboutZ * g2l.getRotation()));
            bend->fixPosition();
        }
        bend->setBendAngle(angle);
    } else {
        bend->setFieldAmplitude(k0, k0s);
    }
        bend->setEntranceAngle(e1);

    if(itsAttr[ROTATION])
        throw OpalException("OpalRBend::update",
                            "ROTATION not supported any more; use PSI instead");


    if(itsAttr[FMAPFN])
        bend->setFieldMapFN(Attributes::getString(itsAttr[FMAPFN]));
    else if(bend->getName() != "RBEND") {
        ERRORMSG(bend->getName() << ": No filename for a field map given. "
                 "Will assume the default map "
                 "\"1DPROFILE1-DEFAULT\"."
                 << endl);
        bend->setFieldMapFN("1DPROFILE1-DEFAULT");
    }

    // Energy in eV.
    if(itsAttr[DESIGNENERGY]) {
        bend->setDesignEnergy(Attributes::getReal(itsAttr[DESIGNENERGY]), false);
    }

    bend->setFullGap(Attributes::getReal(itsAttr[GAP]));

    // if(itsAttr[APERT])
    //     throw OpalException("OpalRBend::fillRegisteredAttributes",
    //                         "APERTURE in RBEND not supported; use GAP and HAPERT instead");

    if(itsAttr[HAPERT]) {
        double hapert = Attributes::getReal(itsAttr[HAPERT]);
        bend->setAperture(ElementBase::RECTANGULAR, std::vector<double>({hapert, hapert, 1.0}));
    }

    if(itsAttr[LENGTH])
        bend->setLength(Attributes::getReal(itsAttr[LENGTH]));
    else
        bend->setLength(0.0);

    if(itsAttr[WAKEF] && itsAttr[DESIGNENERGY] && owk_m == NULL) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + std::string("_wake"));
        owk_m->initWakefunction(*bend);
        bend->setWake(owk_m->wf_m);
    }

    if(itsAttr[K1])
        bend->setK1(Attributes::getReal(itsAttr[K1]));
    else
        bend->setK1(0.0);

    if(itsAttr[PARTICLEMATERINTERACTION] && parmatint_m == NULL) {
        parmatint_m = (ParticleMatterInteraction::find(Attributes::getString(itsAttr[PARTICLEMATERINTERACTION])))->clone(getOpalName() + std::string("_parmatint"));
        parmatint_m->initParticleMatterInteractionHandler(*bend);
        bend->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(bend);
}