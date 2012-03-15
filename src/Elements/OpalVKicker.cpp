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
#include "AbstractObjects/DoomDB.h"
#include "AbstractObjects/DoomReader.h"
#include "AbstractObjects/DoomWriter.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
// JMJ 18/12/2000 no longer need this, see code commented out below.
//#include "Utilities/Options.h"
#include "BeamlineCore/YCorrectorRep.h"
#include "ComponentWrappers/CorrectorWrapper.h"
#include "Physics/Physics.h"


// Class OpalVKicker
// ------------------------------------------------------------------------

OpalVKicker::OpalVKicker():
    OpalElement(SIZE, "VKICKER",
                "The \"VKICKER\" element defines a closed orbit corrector "
                "acting on the vertical plane.")
{
    itsAttr[KICK] = Attributes::makeReal
        ("KICK", "Vertical deflection in rad");

    registerRealAttribute("VKICK");

    setElement((new YCorrectorRep("VKICKER"))->makeWrappers());
}


OpalVKicker::OpalVKicker(const string &name, OpalVKicker *parent):
    OpalElement(name, parent)
{
    setElement((new YCorrectorRep(name))->makeWrappers());
}


OpalVKicker::~OpalVKicker()
{}


OpalVKicker *OpalVKicker::clone(const string &name)
{
    return new OpalVKicker(name, this);
}


void OpalVKicker::doomGet(const DoomReader &reader)
{
    // Read the kicker length.
    itsAttr[LENGTH].getHandler().doomGet(itsAttr[LENGTH], reader,
                                         DoomDB::getAttributeIndex("L"));

    // Read the corrector strength(s).
    int vkick = DoomDB::getAttributeIndex("K0S");
    itsAttr[KICK].getHandler().doomGet(itsAttr[KICK], reader, vkick);
}


void OpalVKicker::doomPut(DoomWriter &writer) const
{
    // Write the kicker length.
    itsAttr[LENGTH].getHandler().doomPut(itsAttr[LENGTH], writer,
                                         DoomDB::getAttributeIndex("L"));

    // Write the kicker strength.
    int vkick = DoomDB::getAttributeIndex("K0S");
    itsAttr[KICK].getHandler().doomPut(itsAttr[KICK], writer, vkick);
}



// JMJ 18/12/2000 Following method not needed, commented out, delete after next CVS commit.
//BEGIN JMJ 14/10/2000 adding new print method by analogy with OpalSextupole.cc
// This is a block of new code

//void OpalVKicker::print(std::ostream &os) const
//{
//  if (Options::opal8) {
//    os << "JMJdebug OpalVKicker OPAL8";
//    Object *parent = getParent();
//    string head = getOpalName() + ':';
//    if (parent->getParent()) {
//      head += parent->getOpalName();
//    } else {
//      head += "VKICKER";
//    }
//    os << head;
//    int len = head.length();
//
//    if (getLength() != 0.0) {
//      printAttribute(os, "LRAD", itsAttr[LENGTH].getImage(), len);
//    }
//
//    itsAttr[TYPE].print(os, len);
//    os << "in OpalVKicker::print" << endl ;
//    printAttribute(os, "reallyVKICK", itsAttr[KICK].getImage(), len);
//    //JMJ 30/11/2000 added harmless semi-colon in OPAL8 output
//    os << ';' << std::endl;
//  } else {
//    OpalElement::print(os);
//  }
//}
//
//END   JMJ 14/10/2000 adding new print method by analogy with OpalSextupole.cc



void OpalVKicker::
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


void OpalVKicker::update()
{
    YCorrectorRep *corr =
        dynamic_cast<YCorrectorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double factor = OPAL.getP0() / Physics::c;
    double kick = Attributes::getReal(itsAttr[KICK]);
    corr->setElementLength(length);
    corr->setBx(kick * factor);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(corr);
}
