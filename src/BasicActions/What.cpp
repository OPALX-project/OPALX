// ------------------------------------------------------------------------
// $RCSfile: What.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: What
//   The class for the OPAL WHAT command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "BasicActions/What.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include <iostream>


// Class What
// ------------------------------------------------------------------------

What::What():
    Action(1, "WHAT",
           "The \"WHAT\" statement displays the definition and attribute"
           " values of an object.") {
    itsAttr[0] = Attributes::makeString
                 ("NAME", "Name of object to be displayed");

    registerOwnership(AttributeHandler::STATEMENT);
}


What::What(const std::string &name, What *parent):
    Action(name, parent)
{}


What::~What()
{}


What *What::clone(const std::string &name) {
    return new What(name, this);
}


void What::execute() {
    if(itsAttr[0]) {
        std::string name = Attributes::getString(itsAttr[0]);

        if(Object *object = OpalData::getInstance()->find(name)) {
            std::cerr << *object << std::endl;
        } else {
            std::cerr << '\n' << *this << "\nUnknown object \"" << name << "\".\n"
                      << std::endl;
        }
    } else {
        printHelp(std::cerr);
    }
}


void What::parse(Statement &statement) {
    parseShortcut(statement);
}