// ------------------------------------------------------------------------
// $RCSfile: System.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: System
//   The class for the OPAL SYSTEM command.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/01/17 22:18:36 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "BasicActions/System.h"

#include "Ippl.h"
#include "Attributes/Attributes.h"

#include <cstdlib>

// Class System
// ------------------------------------------------------------------------

System::System():
    Action(1, "SYSTEM",
           "The \"SYSTEM\" statement sends a command string to the "
           "operating system.") {
    itsAttr[0] = Attributes::makeString
                 ("CMD", "A system command to be executed");

    registerOwnership(AttributeHandler::COMMAND);
}


System::System(const std::string &name, System *parent):
    Action(name, parent)
{}


System::~System()
{}


System *System::clone(const std::string &name) {
    return new System(name, this);
}

/// node 0 is executing only
void System::execute() {
  if(Ippl::myNode() == 0) {
    int res = system(Attributes::getString(itsAttr[0]).c_str());
    if (res!=0)
      ERRORMSG("System call failed" << endl);
  }
}

void System::parse(Statement &statement) {
    parseShortcut(statement);
}
