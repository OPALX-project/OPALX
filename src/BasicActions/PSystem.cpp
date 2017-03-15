// ------------------------------------------------------------------------
// $RCSfile: PSystem.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: PSystem
//   The class for the OPAL SYSTEM command.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/01/17 22:18:36 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "BasicActions/PSystem.h"

#include "Ippl.h"
#include "Attributes/Attributes.h"
#include "AbstractObjects/Object.h"
#include "AbstractObjects/OpalData.h"

#include <cstdlib>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
// Class PSystem
// ------------------------------------------------------------------------

PSystem::PSystem():
    Action(1, "PSYSTEM",
           "The \"PSYSTEM\" statement sends a command string to the "
           "operating system from all nodes.") {
    itsAttr[0] = Attributes::makeString
                 ("CMD", "A system command to be executed");

    registerOwnership(AttributeHandler::COMMAND);
}


PSystem::PSystem(const std::string &name, PSystem *parent):
    Action(name, parent)
{}


PSystem::~PSystem()
{}


PSystem *PSystem::clone(const std::string &name) {
    return new PSystem(name, this);
}


void PSystem::execute() {
    std::string command = Attributes::getString(itsAttr[0]);
    boost::regex variablesRe("\\$\\$\\{([^\\}]+?)\\}");
    boost::smatch match;
    OpalData *opal = OpalData::getInstance();
    while (boost::regex_search(command, match, variablesRe)) {
        Object *object = 0;

        std::string variableName = match[1];
        boost::to_upper(variableName);
        if ((object = opal->find(variableName)) == 0) {
            ERRORMSG("SYSTEM: Could not replace $${" << match[1] << "}, "
                     << "unknown variable" << endl);

            std::string unknownReplacement = "\\\\$\\\\${" + match[1] + "}";
            boost::regex unknownRe("\\$\\$\\{" + match[1] + "\\}");
            command = boost::regex_replace(command, unknownRe, unknownReplacement);
            continue;
        }
        std::ostringstream os;
        object->printValue(os);
        boost::regex variableRe("\\$\\$\\{" + match[1] + "\\}");
        command = boost::regex_replace(command, variableRe, os.str());
    }

    int res = system(command.c_str());
    if (res!=0)
        ERRORMSG("SYSTEM call failed" << endl);
}

void PSystem::parse(Statement &statement) {
    parseShortcut(statement);
}