// ------------------------------------------------------------------------
// $RCSfile: OpalOctupole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalOctupole
//   The class of OPAL Octupoles.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:07 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Elements/OpalOctupole.h"
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

// Class OpalOctupole
// ------------------------------------------------------------------------

OpalOctupole::OpalOctupole():
    OpalElement(SIZE, "OCTUPOLE",
                "The \"OCTUPOLE\" element defines a Octupole.") {
    itsAttr[K3] = Attributes::makeReal
                  ("K3", "Normalised upright octupole coefficient in m^(-4)");
    itsAttr[K3S] = Attributes::makeReal
                   ("K3S", "Normalised skew octupole coefficient in m^(-4)");
    itsAttr[DX] = Attributes::makeReal
      ("DX", "Misalignment in x direction",0.0);
    itsAttr[DY] = Attributes::makeReal
      ("DY", "Misalignment in y direction",0.0);

    registerRealAttribute("DX");
    registerRealAttribute("DY");

    setElement((new MultipoleRep("OCTUPOLE"))->makeWrappers());
}


OpalOctupole::OpalOctupole(const std::string &name, OpalOctupole *parent):
    OpalElement(name, parent) {
    setElement((new MultipoleRep(name))->makeWrappers());
}


OpalOctupole::~OpalOctupole()
{}


OpalOctupole *OpalOctupole::clone(const std::string &name) {
    return new OpalOctupole(name, this);
}


void OpalOctupole::print(std::ostream &os) const {
    OpalElement::print(os);
}


void OpalOctupole::
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


void OpalOctupole::update() {
    MultipoleRep *oct =
        dynamic_cast<MultipoleRep *>(getElement()->removeWrappers());
    oct->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    double factor = OpalData::getInstance()->getP0() / (Physics::c * 6.0);
    BMultipoleField field;
    field.setNormalComponent(4, factor * Attributes::getReal(itsAttr[K3]));
    field.setSkewComponent(4, factor * Attributes::getReal(itsAttr[K3S]));
    oct->setField(field);

    oct->setNormalComponent(2, Attributes::getReal(itsAttr[K3]));
    oct->setSkewComponent(2, Attributes::getReal(itsAttr[K3S]));

    double dx = Attributes::getReal(itsAttr[DX]);
    double dy = Attributes::getReal(itsAttr[DY]);
    oct->setMisalignment(dx, dy, 0.0);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(oct);
}
