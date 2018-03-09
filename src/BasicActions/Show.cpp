// ------------------------------------------------------------------------
// $RCSfile: Show.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Show
//   The class for OPAL SHOW commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "BasicActions/Show.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "Attributes/Attributes.h"
#include <iostream>


// Class Show
// ------------------------------------------------------------------------

Show::Show():
    Action(1, "SHOW",
           "The \"SHOW\" statement displays all object names matching a "
           " regular expression.") {
    itsAttr[0] = Attributes::makeString
                 ("PATTERN", "Regular expression for pattern match");

    registerOwnership(AttributeHandler::STATEMENT);
}


Show::Show(const std::string &name, Show *parent):
    Action(name, parent)
{}


Show::~Show()
{}


Show *Show::clone(const std::string &name) {
    return new Show(name, this);
}


void Show::execute() {
    if(itsAttr[0]) {
        std::string pattern = Attributes::getString(itsAttr[0]);
        OpalData::getInstance()->printNames(std::cerr, pattern);
    } else {
        printHelp(std::cerr);
    }
}


void Show::parse(Statement &statement) {
    parseShortcut(statement);
}