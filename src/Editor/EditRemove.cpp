// ------------------------------------------------------------------------
// $RCSfile: EditRemove.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditRemove
//   The class for the OPAL sequence editor REMOVE command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Editor/EditRemove.h"
#include "AbstractObjects/PlaceRep.h"
#include "Attributes/Attributes.h"
#include "Editor/Edit.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include <iostream>


// Class EditRemove
// ------------------------------------------------------------------------

// The attributes of class EditRemove.
namespace {
  enum {
    SELECTED,  // If true, remove selected elements.
    CLASS,     // The class to be removed.
    SIZE
  };
}


EditRemove::EditRemove():
  Editor(SIZE, "REMOVE",
	 "The \"REMOVE\" sub-command removes element(s) from the sequence "
	 "being edited.")
{
  itsAttr[SELECTED] = Attributes::makeBool
    ("SELECTED", "If true, all selected elements are removed");
  itsAttr[CLASS] = Attributes::makePlace
    ("CLASS", "Name of element class to be removed");
}


EditRemove::EditRemove(const string &name, EditRemove *parent):
  Editor(name, parent)
{}


EditRemove::~EditRemove()
{}


EditRemove *EditRemove::clone(const string &name)
{
  return new EditRemove(name, this);
}


void EditRemove::execute()
{
  int count = 0;
  
  if (Attributes::getBool(itsAttr[SELECTED])) {
    count = Edit::block->removeMultiple();
  } else if (itsAttr[CLASS]) {
    const PlaceRep pos = Attributes::getPlace(itsAttr[CLASS]);
    count = Edit::block->removeSingle(pos);
  } else {
    throw OpalException("EditRemove::execute()",
		       "\"class=selected\" or \"class=<name>\" "
		       "is required for \"REMOVE\".");
  }
  
  if (Options::info) {
    if (count == 0) {
      std::cerr << "\nNo elements";
    } else if (count == 1) {
      std::cerr << "\n1 element";
    } else {
      std::cerr << '\n' << count << " elements";
    }

    std::cerr << " removed.\n" << std::endl;
  }
}
