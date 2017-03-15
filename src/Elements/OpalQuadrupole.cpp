// ------------------------------------------------------------------------
// $RCSfile: OpalQuadrupole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalQuadrupole
//   The class of OPAL Quadrupoles.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:07 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Elements/OpalQuadrupole.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MultipoleRep.h"
#include "ComponentWrappers/MultipoleWrapper.h"
#include "Fields/BMultipoleField.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include "Structure/SurfacePhysics.h"

#include <cmath>
#include <iostream>
#if defined(__GNUC__) && __GNUC__ < 3
#include <strstream>
#else
#include <sstream>
#endif


// Class OpalQuadrupole
// ------------------------------------------------------------------------

OpalQuadrupole::OpalQuadrupole():
    OpalElement(SIZE, "QUADRUPOLE",
                "The \"QUADRUPOLE\" element defines a Quadrupole."),
    sphys_m(NULL) {
    itsAttr[K1] = Attributes::makeReal
                  ("K1", "Normalised upright quadrupole coefficient in m^(-2)");
    itsAttr[DK1] = Attributes::makeReal
                  ("DK1", "Normalised upright quadrupole coefficient error in m^(-2)");
    itsAttr[K1S] = Attributes::makeReal
                   ("K1S", "Normalised skew quadrupole coefficient in m^(-2)");

    itsAttr[DK1S] = Attributes::makeReal
                   ("DK1S", "Normalised skew quadrupole coefficient error in m^(-2)");

    registerRealAttribute("K1");
    registerRealAttribute("DK1");
    registerRealAttribute("K1S");
    registerRealAttribute("DK1S");

    registerOwnership();

    setElement((new MultipoleRep("QUADRUPOLE"))->makeWrappers());
}


OpalQuadrupole::OpalQuadrupole(const std::string &name, OpalQuadrupole *parent):
    OpalElement(name, parent),
    sphys_m(NULL) {
    setElement((new MultipoleRep(name))->makeWrappers());
}


OpalQuadrupole::~OpalQuadrupole() {
    if(sphys_m)
        delete sphys_m;
}


OpalQuadrupole *OpalQuadrupole::clone(const std::string &name) {
    return new OpalQuadrupole(name, this);
}


void OpalQuadrupole::print(std::ostream &os) const {
    OpalElement::print(os);
}


void OpalQuadrupole::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    // Get the desired field.
    const MultipoleWrapper *mult =
        dynamic_cast<const MultipoleWrapper *>(base.removeAlignWrapper());
    BMultipoleField field;

    // Get the desired field.
    if(flag == ERROR_FLAG) {
        field = mult->errorField();
    } else if(flag == ACTUAL_FLAG) {
        field = mult->getField();
    } else if(flag == IDEAL_FLAG) {
        field = mult->getDesign().getField();
    }

    double length = getLength();
    double scale = Physics::c / OpalData::getInstance()->getP0();
    if(length != 0.0) scale *= length;

    for(int order = 1; order <= field.order(); ++order) {
#if defined(__GNUC__) && __GNUC__ < 3
        char buffer[10];
        std::ostrstream ss(buffer, 10);
#else
        std::ostringstream ss;
#endif
        ss << (order - 1) << std::ends;
#if defined(__GNUC__) && __GNUC__ < 3
        std::string orderString(buffer);
#else
        std::string orderString = ss.str();
#endif

        std::string normName = "K" + orderString + "L";
        registerRealAttribute(normName)->setReal(scale * field.normal(order));

        std::string skewName = "K" + orderString + "SL";
        registerRealAttribute(skewName)->setReal(scale * field.skew(order));

        scale *= double(order);
    }
}


void OpalQuadrupole::update() {
    OpalElement::update();

    MultipoleRep *quad =
        dynamic_cast<MultipoleRep *>(getElement()->removeWrappers());
    quad->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    double factor = OpalData::getInstance()->getP0() / Physics::c;

    BMultipoleField field;
    field.setNormalComponent(2, factor * Attributes::getReal(itsAttr[K1]));   // this is for the maps
    field.setSkewComponent(2, factor * Attributes::getReal(itsAttr[K1S]));    // this is for the maps
    quad->setField(field);
    quad->setNormalComponent(2, Attributes::getReal(itsAttr[K1]), Attributes::getReal(itsAttr[DK1]));
    quad->setSkewComponent(2, Attributes::getReal(itsAttr[K1S]), Attributes::getReal(itsAttr[DK1S]));

    if(itsAttr[SURFACEPHYSICS] && sphys_m == NULL) {
        sphys_m = (SurfacePhysics::find(Attributes::getString(itsAttr[SURFACEPHYSICS])))->clone(getOpalName() + std::string("_sphys"));
        sphys_m->initSurfacePhysicsHandler(*quad);
        quad->setSurfacePhysics(sphys_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(quad);
}