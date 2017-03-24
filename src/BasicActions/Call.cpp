// ------------------------------------------------------------------------
// $RCSfile: Call.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Call
//   The class for OPAL "CALL" commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "BasicActions/Call.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "OpalParser/OpalParser.h"
#include "Parser/FileStream.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include <iostream>

using std::cerr;
using std::endl;


// Class Call
// ------------------------------------------------------------------------

Call::Call():
    Action(1, "CALL",
           "The \"CALL\" statement switches input temporarily to the "
           "named file.") {
    itsAttr[0] =
        Attributes::makeString("FILE", "Name of file to be read", "CALL");

    registerOwnership(AttributeHandler::STATEMENT);
}


Call::Call(const std::string &name, Call *parent):
    Action(name, parent)
{}


Call::~Call()
{}


Call *Call::clone(const std::string &name) {
    return new Call(name, this);
}


void Call::execute() {
    std::string file = Attributes::getString(itsAttr[0]);

    if(Options::info) {
        cerr << endl
             << "Start reading input stream \"" << file << "\"." << endl
             << endl;
    }

    OpalParser().run(new FileStream(file));

    if(Options::info) {
        cerr << endl
             << "End reading input stream \"" << file << "\"." << endl
             << endl;
    }
}


void Call::parse(Statement &statement) {
    parseShortcut(statement);
}
