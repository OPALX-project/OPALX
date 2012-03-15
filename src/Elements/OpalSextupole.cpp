// ------------------------------------------------------------------------
// $RCSfile: OpalSextupole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalSextupole
//   The class of OPAL Sextupoles.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Elements/OpalSextupole.h"
#include "AbstractObjects/DoomDB.h"
#include "AbstractObjects/DoomWriter.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MultipoleRep.h"
#include "ComponentWrappers/MultipoleWrapper.h"
#include "Fields/BMultipoleField.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include <cmath>
#include <iostream>
#if defined(__GNUC__) && __GNUC__ < 3
#include <strstream>
#else
#include <sstream>
#endif


// Class OpalSextupole
// ------------------------------------------------------------------------

OpalSextupole::OpalSextupole():
    OpalElement(SIZE, "SEXTUPOLE",
                "The \"SEXTUPOLE\" element defines a Sextupole.") {
    itsAttr[K2] = Attributes::makeReal
                  ("K2", "Normalised upright sextupole coefficient in m^(-3)");
    itsAttr[K2S] = Attributes::makeReal
                   ("K2S", "Normalised skew sextupole coefficient in m^(-3)");

    setElement((new MultipoleRep("SEXTUPOLE"))->makeWrappers());
}


OpalSextupole::OpalSextupole(const string &name, OpalSextupole *parent):
    OpalElement(name, parent) {
    setElement((new MultipoleRep(name))->makeWrappers());
}


OpalSextupole::~OpalSextupole()
{}


OpalSextupole *OpalSextupole::clone(const string &name) {
    return new OpalSextupole(name, this);
}


void OpalSextupole::doomPut(DoomWriter &writer) const {
    // Save the OPAL-9 data.
    OpalElement::doomPut(writer);

    // Save the tilt angle.
    double k2   = Attributes::getReal(itsAttr[K2]);
    double k2s  = Attributes::getReal(itsAttr[K2S]);

    Attribute tilt(Attributes::makeReal
                   ("TILT", ""));
    Attributes::setReal(tilt, (k2s == 0.0) ? 0.0 : atan2(k2s, k2) / 3.0);
    static int index = DoomDB::getAttributeIndex("TILT");
    tilt.doomPut(writer, index);
}


void OpalSextupole::print(std::ostream &os) const {
        OpalElement::print(os);
}


void OpalSextupole::
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


void OpalSextupole::update() {
    MultipoleRep *sext =
        dynamic_cast<MultipoleRep *>(getElement()->removeWrappers());
    sext->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    double factor = OPAL.getP0() / (Physics::c * 2.0);
    BMultipoleField field;
    field.setNormalComponent(3, factor * Attributes::getReal(itsAttr[K2]));
    field.setSkewComponent(3, factor * Attributes::getReal(itsAttr[K2S]));
    sext->setField(field);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(sext);
}
