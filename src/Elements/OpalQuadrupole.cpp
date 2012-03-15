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
#include "AbstractObjects/DoomDB.h"
#include "AbstractObjects/DoomWriter.h"
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
    itsAttr[K1S] = Attributes::makeReal
                   ("K1S", "Normalised skew quadrupole coefficient in m^(-2)");

    setElement((new MultipoleRep("QUADRUPOLE"))->makeWrappers());
}


OpalQuadrupole::OpalQuadrupole(const string &name, OpalQuadrupole *parent):
    OpalElement(name, parent),
    sphys_m(NULL) {
    setElement((new MultipoleRep(name))->makeWrappers());
}


OpalQuadrupole::~OpalQuadrupole() {
    if(sphys_m)
        delete sphys_m;
}


OpalQuadrupole *OpalQuadrupole::clone(const string &name) {
    return new OpalQuadrupole(name, this);
}


void OpalQuadrupole::doomPut(DoomWriter &writer) const {
    // Save the OPAL-9 data.
    OpalElement::doomPut(writer);

    // Save the tilt angle.
    double k1   = Attributes::getReal(itsAttr[K1]);
    double k1s  = Attributes::getReal(itsAttr[K1S]);

    static Attribute tilt(Attributes::makeReal("TILT", ""));
    Attributes::setReal(tilt, (k1s == 0.0) ? 0.0 : atan2(k1s, k1) / 2.0);
    static int index = DoomDB::getAttributeIndex("TILT");
    tilt.doomPut(writer, index);
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
    double scale = Physics::c / OPAL.getP0();
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
        string orderString(buffer);
#else
        std::string orderString = ss.str();
#endif

        string normName = "K" + orderString + "L";
        registerRealAttribute(normName)->setReal(scale * field.normal(order));

        string skewName = "K" + orderString + "SL";
        registerRealAttribute(skewName)->setReal(scale * field.skew(order));

        scale *= double(order);
    }
}


void OpalQuadrupole::update() {
    MultipoleRep *quad =
        dynamic_cast<MultipoleRep *>(getElement()->removeWrappers());
    quad->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    double factor = OPAL.getP0() / Physics::c;

    BMultipoleField field;
    field.setNormalComponent(2, factor * Attributes::getReal(itsAttr[K1]));   // this is for the maps
    field.setSkewComponent(2, factor * Attributes::getReal(itsAttr[K1S]));    // this is for the maps
    quad->setField(field);
    quad->setNormalComponent(2, Attributes::getReal(itsAttr[K1]));
    quad->setSkewComponent(2, Attributes::getReal(itsAttr[K1S]));

    vector<double> apert = getApert();
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
        sphys_m->initSurfacePhysicsHandler(*quad, apert_major, apert_minor);
        quad->setSurfacePhysics(sphys_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(quad);
}
