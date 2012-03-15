// ------------------------------------------------------------------------
// $RCSfile: OpalKicker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalKicker
//   The class of OPAL horizontal orbit correctors.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:32:23 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Elements/OpalKicker.h"
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/DoomDB.h"
#include "AbstractObjects/DoomReader.h"
#include "AbstractObjects/DoomWriter.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
// JMJ 18/12/2000 no longer need this, see code commented out below.
//#include "Utilities/Options.h"
#include "BeamlineCore/CorrectorRep.h"
#include "ComponentWrappers/CorrectorWrapper.h"
#include "Physics/Physics.h"


// Class OpalKicker
// ------------------------------------------------------------------------

OpalKicker::OpalKicker():
    OpalElement(SIZE, "KICKER",
                "The \"KICKER\" element defines a closed orbit corrector "
                "acting on both planes.") {
    itsAttr[HKICK] = Attributes::makeReal
                     ("HKICK", "Horizontal deflection in rad");
    itsAttr[VKICK] = Attributes::makeReal
                     ("VKICK", "Vertical deflection in rad");

    registerRealAttribute("HKICK");
    registerRealAttribute("VKICK");

    setElement((new CorrectorRep("KICKER"))->makeWrappers());
}


OpalKicker::OpalKicker(const string &name, OpalKicker *parent):
    OpalElement(name, parent) {
    setElement((new CorrectorRep(name))->makeWrappers());
}


OpalKicker::~OpalKicker()
{}


OpalKicker *OpalKicker::clone(const string &name) {
    return new OpalKicker(name, this);
}


void OpalKicker::doomGet(const DoomReader &reader) {
    // Read the kicker length.
    itsAttr[LENGTH].getHandler().doomGet(itsAttr[LENGTH], reader,
                                         DoomDB::getAttributeIndex("L"));

    // Read the corrector strengths.
    int hkick = DoomDB::getAttributeIndex("K0");
    int vkick = DoomDB::getAttributeIndex("K0S");
    itsAttr[HKICK].getHandler().doomGet(itsAttr[HKICK], reader, hkick);
    itsAttr[VKICK].getHandler().doomGet(itsAttr[VKICK], reader, vkick);
}


void OpalKicker::doomPut(DoomWriter &writer) const {
    // Write the kicker length.
    itsAttr[LENGTH].getHandler().doomPut(itsAttr[LENGTH], writer,
                                         DoomDB::getAttributeIndex("L"));

    // Write the corrector strengths.
    int hkick = DoomDB::getAttributeIndex("K0");
    int vkick = DoomDB::getAttributeIndex("K0S");
    itsAttr[HKICK].getHandler().doomPut(itsAttr[HKICK], writer, hkick);
    itsAttr[VKICK].getHandler().doomPut(itsAttr[VKICK], writer, vkick);
}


void OpalKicker::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
    const CorrectorWrapper *corr =
        dynamic_cast<const CorrectorWrapper *>(base.removeAlignWrapper());
    BDipoleField field;

    if(flag == ERROR_FLAG) {
        field = corr->errorField();
    } else if(flag == ACTUAL_FLAG) {
        field = corr->getField();
    } else if(flag == IDEAL_FLAG) {
        field = corr->getDesign().getField();
    }

    double scale = Physics::c / OPAL.getP0();
    attributeRegistry["HKICK"]->setReal(- field.getBy() * scale);
    attributeRegistry["VKICK"]->setReal(+ field.getBx() * scale);
}


void OpalKicker::update() {
    CorrectorRep *corr =
        dynamic_cast<CorrectorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double factor = OPAL.getP0() / Physics::c;
    double hKick = Attributes::getReal(itsAttr[HKICK]);
    double vKick = Attributes::getReal(itsAttr[VKICK]);
    corr->setElementLength(length);
    corr->setBy(- hKick * factor);
    corr->setBx(vKick * factor);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(corr);
}
