// ------------------------------------------------------------------------
// $RCSfile: EditCycle.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditCycle
//   The class for the OPAL sequence editor CYCLE command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Editor/EditCycle.h"
#include "AbstractObjects/PlaceRep.h"
#include "Attributes/Attributes.h"
#include "Editor/Edit.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include <iostream>


// Class EditCycle
// ------------------------------------------------------------------------


EditCycle::EditCycle():
  Editor(1, "CYCLE",
	 "The \"CYCLE\" sub-command changes the starting point "
	 "by cyclic interchange.")
{
  itsAttr[0] = Attributes::makePlace
    ("START", "The new start position for the sequence");
}


EditCycle::EditCycle(const string &name, EditCycle *parent):
  Editor(name, parent)
{}


EditCycle::~EditCycle()
{}


EditCycle *EditCycle::clone(const string &name)
{
  return new EditCycle(name, this);
}


void EditCycle::execute()
{
  const PlaceRep start = Attributes::getPlace(itsAttr[0]);

  if (Edit::block->cycle(start)) {
    if (Options::info) {
      std::cerr << "\nStart position changed to \"" << itsAttr[0] << "\".\n"
		<< std::endl;
    }
  } else {
    throw OpalException("EditCycle::execute()", "Could not find place \"" +
		       itsAttr[0].getImage() + "\" in top-level sequence.");
  }
}
