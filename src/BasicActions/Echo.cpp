// ------------------------------------------------------------------------
// $RCSfile: Echo.cpp,v $
// ------------------------------------------------------------------------
// $Release$
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Echo
//   The class for OPAL ECHO commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "BasicActions/Echo.h"
#include "Attributes/Attributes.h"
#include "Utility/IpplInfo.h"
#include <iostream>


// Class Echo
// ------------------------------------------------------------------------

Echo::Echo():
    Action(1, "ECHO",
           "The \"ECHO\" statement sends a message to the ECHO file.") {
    itsAttr[0] = Attributes::makeString("MESSAGE", "The message to be sent.");

    registerOwnership(AttributeHandler::STATEMENT);
}


Echo::Echo(const std::string &name, Echo *parent):
    Action(name, parent)
{}


Echo::~Echo()
{}


Echo *Echo::clone(const std::string &name) {
    return new Echo(name, this);
}


void Echo::execute() {
    if (Ippl::myNode() == 0)
        std::cerr << Attributes::getString(itsAttr[0]) << std::endl;
}


void Echo::parse(Statement &statement) {
    parseShortcut(statement);
}