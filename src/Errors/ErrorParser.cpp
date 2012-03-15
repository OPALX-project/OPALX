// ------------------------------------------------------------------------
// $RCSfile: ErrorParser.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorParser
//   The parser class for the OPAL error module.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/04/19 12:06:04 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "Errors/ErrorParser.h"
#include "Errors/ErrorAlign.h"
#include "AbstractObjects/OpalData.h"  // JMJ: added 10/4/2000
#include "Errors/ErrorComp.h"
#include "Errors/ErrorEnd.h"
#include "Errors/ErrorField.h"
#include "Errors/ErrorPrint.h"
#include "Errors/ErrorSave.h"
#include "Errors/ErrorSelect.h"
#include "BasicActions/Option.h"


// Class ErrorParser
// ------------------------------------------------------------------------


ErrorParser::ErrorParser():
    ErrorDirectory() {
    // Error defining commands.
    ErrorDirectory.insert("EALIGN",     new ErrorAlign());
    ErrorDirectory.insert("EFCOMP",     new ErrorComp());
    ErrorDirectory.insert("ENDERROR",   new ErrorEnd());
    ErrorDirectory.insert("EFIELD",     new ErrorField());
    ErrorDirectory.insert("EPRINT",     new ErrorPrint());
    ErrorDirectory.insert("ESAVE",      new ErrorSave());
    ErrorDirectory.insert("SELECT",     new ErrorSelect());
    ErrorDirectory.insert("OPTION",     new Option());
    // JMJ 10/4/2000: added following commands here, seems harmless
    ErrorDirectory.insert("VALUE", OPAL.find("VALUE"));
    ErrorDirectory.insert("SYSTEM", OPAL.find("SYSTEM"));
    //   MatchDirectory.insert("OPTION", OPAL.find("OPTION")); // does not work
}


ErrorParser::~ErrorParser()
{}


Object *ErrorParser::find(const string &name) const {
    return ErrorDirectory.find(name);
}
