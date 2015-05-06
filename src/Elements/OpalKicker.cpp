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


OpalKicker::OpalKicker(const std::string &name, OpalKicker *parent):
    OpalElement(name, parent) {
    setElement((new CorrectorRep(name))->makeWrappers());
}


OpalKicker::~OpalKicker()
{}


OpalKicker *OpalKicker::clone(const std::string &name) {
    return new OpalKicker(name, this);
}


void OpalKicker::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    Inform m("fillRegisteredAttributes ");

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

    double scale = Physics::c / OpalData::getInstance()->getP0();
    attributeRegistry["HKICK"]->setReal(- field.getBy() * scale);
    attributeRegistry["VKICK"]->setReal(+ field.getBx() * scale);

    m << "P= " << OpalData::getInstance()->getP0()
      << " Bx= " << field.getBx()
      << " By= " << field.getBy() << endl;
}


void OpalKicker::update() {
    CorrectorRep *corr =
        dynamic_cast<CorrectorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double factor = OpalData::getInstance()->getP0() / Physics::c;
    double hKick = Attributes::getReal(itsAttr[HKICK]);
    double vKick = Attributes::getReal(itsAttr[VKICK]);
    corr->setElementLength(length);
    corr->setBy(- hKick * factor);
    corr->setBx(vKick * factor);

    corr->setKickX(Attributes::getReal(itsAttr[HKICK]));
    corr->setKickY(Attributes::getReal(itsAttr[VKICK]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(corr);
}