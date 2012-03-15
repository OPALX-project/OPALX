// ------------------------------------------------------------------------
// $RCSfile: MPWriter.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: FieldReader
//   Ancillary class for writing align errors to DOOM data base.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/MPWriter.h"
#include "AbstractObjects/DoomDB.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "Fields/BMultipoleField.h"


// Class MPWriter
// ------------------------------------------------------------------------


MPWriter::MPWriter(const string &name)
{
  // Set DOOM environment for field error.
  DOOM_DB.setEnvironment("USE", name);
}


MPWriter::~MPWriter()
{}


void MPWriter::fieldError(const string &name, int occur,
			  const BMultipoleField &designField,
			  BMultipoleField &errorField)
{
  // Extract the element length.
  Element *elem = Element::find(name);
  double length = elem->getLength();
  if (length == 0.0) length = 1.0;

  // Write the field error.
  DOOM_DB.writeField(name, occur, length, designField, errorField);
}
