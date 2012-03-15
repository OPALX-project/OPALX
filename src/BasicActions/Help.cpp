// ------------------------------------------------------------------------
// $RCSfile: Help.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Help
//   The class for OPAL HELP commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "BasicActions/Help.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "Attributes/Attributes.h"
#include <iostream>


// Class Help
// ------------------------------------------------------------------------

Help::Help():
  Action(1, "HELP",
	 "The \"HELP\" statement displays the purpose and attribute "
	 "types of an object.")
{
  itsAttr[0] =
    Attributes::makeString("NAME", "Name of object for which help is wanted");
}


Help::Help(const string &name, Help *parent):
  Action(name, parent)
{}


Help::~Help()
{}


Help *Help::clone(const string &name)
{
  return new Help(name, this);
}


void Help::execute()
{
  if (itsAttr[0]) {
    string name = Attributes::getString(itsAttr[0]);

    if (Object *object = OPAL.find(name)) {
      object->printHelp(std::cerr);
    } else {
      std::cerr << '\n'	<< *this << "Unknown object \"" << name << "\".\n"
		<< std::endl ;
    }
  } else {
    printHelp(std::cerr);
  }
}


void Help::parse(Statement &statement)
{
  parseShortcut(statement);
}
