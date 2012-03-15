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

// Class OpalOctupole
// ------------------------------------------------------------------------

OpalOctupole::OpalOctupole():
  OpalElement(SIZE, "OCTUPOLE",
	     "The \"OCTUPOLE\" element defines a Octupole.")
{
  itsAttr[K3] = Attributes::makeReal
    ("K3", "Normalised upright octupole coefficient in m^(-4)");
  itsAttr[K3S] = Attributes::makeReal
    ("K3S", "Normalised skew octupole coefficient in m^(-4)");

  setElement((new MultipoleRep("OCTUPOLE"))->makeWrappers());
}


OpalOctupole::OpalOctupole(const string &name, OpalOctupole *parent):
  OpalElement(name, parent)
{
  setElement((new MultipoleRep(name))->makeWrappers());
}


OpalOctupole::~OpalOctupole()
{}


OpalOctupole *OpalOctupole::clone(const string &name)
{
  return new OpalOctupole(name, this);
}


void OpalOctupole::doomPut(DoomWriter &writer) const
{
  // Save the OPAL-9 data.
  OpalElement::doomPut(writer);

  // Save the tilt angle.
  double k3   = Attributes::getReal(itsAttr[K3]);
  double k3s  = Attributes::getReal(itsAttr[K3S]);

  static Attribute tilt(Attributes::makeReal("TILT", ""));
  Attributes::setReal(tilt, (k3s == 0.0) ? 0.0 : atan2(k3s, k3) / 4.0);
  static int index = DoomDB::getAttributeIndex("TILT");
  tilt.doomPut(writer, index);
}


void OpalOctupole::print(std::ostream &os) const
{
  if (Options::opal8) {
    Object *parent = getParent();
    string head = getOpalName() + ':';
    if (parent->getParent()) {
      head += parent->getOpalName();
    } else {
      head += "MULTIPOLE";
    }
    os << head;
    int len = head.length();

    if (getLength() != 0.0) {
      printAttribute(os, "LRAD", itsAttr[LENGTH].getImage(), len);
    }
    
    itsAttr[TYPE].print(len);
    printMultipoleStrength(os, 3, len, "K3L", "T3",
			   itsAttr[LENGTH], itsAttr[K3], itsAttr[K3S]);
    os << std::endl;
  } else {
    OpalElement::print(os);
  }
}


void OpalOctupole::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag)
{
  OpalElement::fillRegisteredAttributes(base, flag);

  // Get the desired field.
  const MultipoleWrapper *mult =
    dynamic_cast<const MultipoleWrapper *>(base.removeAlignWrapper());
  BMultipoleField field;

  // Get the desired field.
  if (flag == ERROR_FLAG) {
    field = mult->errorField();
  } else if (flag == ACTUAL_FLAG) {
    field = mult->getField();
  } else if (flag == IDEAL_FLAG) {
    field = mult->getDesign().getField();
  }

  double length = getLength();
  double scale = Physics::c / OPAL.getP0();
  if (length != 0.0) scale *= length;

  for (int order = 1; order <= field.order(); ++order) {
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
    std::string orderString=ss.str();
#endif

    string normName = "K" + orderString + "L";
    registerRealAttribute(normName)->setReal(scale * field.normal(order));

    string skewName = "K" + orderString + "SL";
    registerRealAttribute(skewName)->setReal(scale * field.skew(order));

    scale *= double(order);
  }
}


void OpalOctupole::update()
{
  MultipoleRep *oct =
    dynamic_cast<MultipoleRep *>(getElement()->removeWrappers());
  oct->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
  double factor = OPAL.getP0() / (Physics::c * 6.0);
  BMultipoleField field;
  field.setNormalComponent(4, factor * Attributes::getReal(itsAttr[K3]));
  field.setSkewComponent  (4, factor * Attributes::getReal(itsAttr[K3S]));
  oct->setField(field);

  // Transmit "unknown" attributes.
  OpalElement::updateUnknown(oct);
}
