// ------------------------------------------------------------------------
// $RCSfile: Dump.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Dump
//   The class for the OPAL "DUMP" command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Utility/IpplInfo.h"
#include "BasicActions/Dump.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/ObjectFunction.h"
#include <iostream>


// Ancillary class DumpObject
// Dump an object to standard error.
// ------------------------------------------------------------------------

struct DumpObject: ObjectFunction {
    virtual void operator()(Object *) const;
};


void DumpObject::operator()(Object *object) const {
    ERRORMSG(*object << endl);
}


// Class Dump
// ------------------------------------------------------------------------

Dump::Dump():
    Action(0, "DUMP",
           "The \"DUMP\" statement dumps the data structure of OPAL.")
{}


Dump::Dump(const string &name, Dump *parent):
    Action(name, parent)
{}


Dump::~Dump()
{}


void Dump::execute() {
    // Dump complete directory in alphabetical order to cerr.
    OpalData::getInstance()->apply(DumpObject());
}


Dump *Dump::clone(const string &name) {
    return new Dump(name, this);
}
