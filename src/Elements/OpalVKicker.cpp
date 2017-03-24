// ------------------------------------------------------------------------
// $RCSfile: OpalVKicker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalVKicker
//   The class of OPAL vertical orbit correctors.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:32:24 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Elements/OpalVKicker.h"
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/YCorrectorRep.h"
#include "ComponentWrappers/CorrectorWrapper.h"
#include "Physics/Physics.h"


// Class OpalVKicker
// ------------------------------------------------------------------------

OpalVKicker::OpalVKicker():
    OpalElement(SIZE, "VKICKER",
                "The \"VKICKER\" element defines a closed orbit corrector "
                "acting on the vertical plane.") {
    itsAttr[KICK] = Attributes::makeReal
                    ("KICK", "Vertical deflection in rad");

    registerRealAttribute("VKICK");

    registerOwnership();

    setElement((new YCorrectorRep("VKICKER"))->makeWrappers());
}


OpalVKicker::OpalVKicker(const std::string &name, OpalVKicker *parent):
    OpalElement(name, parent) {
    setElement((new YCorrectorRep(name))->makeWrappers());
}


OpalVKicker::~OpalVKicker()
{}


OpalVKicker *OpalVKicker::clone(const std::string &name) {
    return new OpalVKicker(name, this);
}


void OpalVKicker::
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


void OpalVKicker::update() {
    OpalElement::update();

    YCorrectorRep *corr = dynamic_cast<YCorrectorRep *>(getElement()->removeWrappers());

    double length = Attributes::getReal(itsAttr[LENGTH]);
    double factor = OpalData::getInstance()->getP0() / Physics::c;
    double kick = Attributes::getReal(itsAttr[KICK]);

    corr->setElementLength(length);
    corr->setBx(kick * factor);
    corr->setKickY(Attributes::getReal(itsAttr[KICK]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(corr);
}