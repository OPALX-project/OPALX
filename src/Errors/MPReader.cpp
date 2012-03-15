// ------------------------------------------------------------------------
// $RCSfile: MPReader.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MPReader
//   Ancillary class for reading field errors from DOOM data base.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/MPReader.h"
#include "AbstractObjects/DoomDB.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "Fields/BMultipoleField.h"


// Class MPReader
// ------------------------------------------------------------------------


MPReader::MPReader(const string &name) {
    // Set DOOM environment for field error.
    DOOM_DB.setEnvironment("USE", name);
}


MPReader::~MPReader()
{}


void MPReader::fieldError(const string &name, int occur,
                          const BMultipoleField &designField,
                          BMultipoleField &errorField) {
    // Extract the element length.
    Element *elem = Element::find(name);
    double length = elem->getLength();
    if(length == 0.0) length = 1.0;

    // Read the field error.
    DOOM_DB.readField(name, occur, length, designField, errorField);
}
