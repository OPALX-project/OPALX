// ------------------------------------------------------------------------
// $RCSfile: RealVariable.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RealVariable
//   Implements a OPAL REAL_VARIABLE definition.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:26:42 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "ValueDefinitions/RealVariable.h"
#include "AbstractObjects/DoomWriter.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Utilities/Options.h"
#include <iostream>


// Class RealVariable
// ------------------------------------------------------------------------

RealVariable::RealVariable():
  ValueDefinition(1, "REAL_VARIABLE",
		  "The \"REAL VARIABLE\" statement defines a global "
		  "real variable:\n"
		  "\tREAL VARIABLE <name>=<real-expression>;\n")
{
  itsAttr[0] = Attributes::makeReal("VALUE", "The variable value", 0.0);

  // Construct the P0 variable.
  RealVariable *p0 = new RealVariable("P0", this, 1.0);
  OPAL.create(p0);
  p0->setDirty(true);
  OPAL.setP0(p0);
}


RealVariable::RealVariable(const string &name, RealVariable *parent,
			   double value):
  ValueDefinition(name, parent)
{
  Attributes::setReal(itsAttr[0], value);
}


RealVariable::RealVariable(const string &name, RealVariable *parent):
  ValueDefinition(name, parent)
{}


RealVariable::~RealVariable()
{}


bool RealVariable::canReplaceBy(Object *object)
{
  // Replace only by another variable.
  return (dynamic_cast<RealVariable*>(object) != 0);
}


RealVariable *RealVariable::clone(const string &name)
{
  return new RealVariable(name, this);
}


double RealVariable::getReal() const
{
  return Attributes::getReal(itsAttr[0]);
}


void RealVariable::doomGet(const DoomReader &reader)
{
  itsAttr[0].doomGet(reader, 0);
}


void RealVariable::doomPut(DoomWriter &writer) const
{
  itsAttr[0].doomPut(writer, 0);
}


//JMJ 30/11/2000 added harmless semi-colon in OPAL8 output
void RealVariable::print(std::ostream &os) const
{
  if (Options::opal8) {
    os << getOpalName() << ":=" << itsAttr[0] << ';';
  } else {
    os << "REAL " << getOpalName()
       << (itsAttr[0].isExpression() ? ":=" : "=") << itsAttr[0] << ';';
  }
  os << std::endl;
}
