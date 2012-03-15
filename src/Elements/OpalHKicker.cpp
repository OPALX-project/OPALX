// ------------------------------------------------------------------------
// $RCSfile: OpalHKicker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalHKicker
//   The class of OPAL horizontal orbit correctors.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalHKicker.h"
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/DoomDB.h"
#include "AbstractObjects/DoomWriter.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/XCorrectorRep.h"
#include "ComponentWrappers/CorrectorWrapper.h"
#include "Physics/Physics.h"


// Class OpalHKicker
// ------------------------------------------------------------------------

OpalHKicker::OpalHKicker():
    OpalElement(SIZE, "HKICKER",
                "The \"HKICKER\" element defines a closed orbit corrector "
                "acting on the horizontal plane.")
{
    itsAttr[KICK] = Attributes::makeReal
        ("KICK", "Horizontal deflection in rad");

    registerRealAttribute("HKICK");

    setElement((new XCorrectorRep("HKICKER"))->makeWrappers());
}


OpalHKicker::OpalHKicker(const string &name, OpalHKicker *parent):
    OpalElement(name, parent)
{
    setElement((new XCorrectorRep(name))->makeWrappers());
}


OpalHKicker::~OpalHKicker()
{}


OpalHKicker *OpalHKicker::clone(const string &name)
{
    return new OpalHKicker(name, this);
}


void OpalHKicker::doomGet(const DoomReader &reader)
{
    // Read the kicker length.
    itsAttr[LENGTH].getHandler().doomGet(itsAttr[LENGTH], reader,
                                         DoomDB::getAttributeIndex("L"));

    // Read the corrector strength(s).
    int hkick = DoomDB::getAttributeIndex("K0");
    itsAttr[KICK].getHandler().doomGet(itsAttr[KICK], reader, hkick);
}


void OpalHKicker::doomPut(DoomWriter &writer) const
{
    // Write the kicker length.
    itsAttr[LENGTH].getHandler().doomPut(itsAttr[LENGTH], writer,
                                         DoomDB::getAttributeIndex("L"));

    // Write the kicker strength.
    int hkick = DoomDB::getAttributeIndex("K0");
    itsAttr[KICK].getHandler().doomPut(itsAttr[KICK], writer, hkick);
}


void OpalHKicker::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag)
{
    OpalElement::fillRegisteredAttributes(base, flag);
    const CorrectorWrapper *corr =
        dynamic_cast<const CorrectorWrapper *>(base.removeAlignWrapper());
    BDipoleField field;

    if (flag == ERROR_FLAG) {
        field = corr->errorField();
    } else if (flag == ACTUAL_FLAG) {
        field = corr->getField();
    } else if (flag == IDEAL_FLAG) {
        field = corr->getDesign().getField();
    }

    double scale = Physics::c / OPAL.getP0();
    attributeRegistry["HKICK"]->setReal(- field.getBy() * scale);
    attributeRegistry["VKICK"]->setReal(+ field.getBx() * scale);
}


void OpalHKicker::update()
{
    XCorrectorRep *corr =
        dynamic_cast<XCorrectorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double factor = OPAL.getP0() / Physics::c;
    double kick = Attributes::getReal(itsAttr[KICK]);
    corr->setElementLength(length);
    corr->setBy(- kick * factor);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(corr);
}
