// ------------------------------------------------------------------------
// $RCSfile: VaryCmd.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: VaryCmd
//   The class for the OPAL VARY command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Match/VaryCmd.h"
#include "Attributes/Attributes.h"
#include "Match/ConcreteVar.h"
#include "Match/Match.h"
#include "Match/MatchLimits.h"
#include "Utilities/OpalOptions.h"
#include <iostream>


// Class VaryCmd
// ------------------------------------------------------------------------

// The attributes of class VaryCmd.
namespace {
    enum {
        NAME,    // the variable reference
        STEP,    // the step size
        LOWER,   // the lower limit, if any
        UPPER,   // the upper limit, if any
        SIZE
    };
}


VaryCmd::VaryCmd():
    Action(SIZE, "VARY",
           "The \"VARY\" sub-command defines a variable to be adjusted.") {
    itsAttr[NAME] = Attributes::makeReference
                    ("NAME", "the variable reference");
    itsAttr[STEP] = Attributes::makeReal
                    ("STEP", "the step size", 0.2);
    itsAttr[LOWER] = Attributes::makeReal
                     ("LOWER", "the lower limit, if any");
    itsAttr[UPPER] = Attributes::makeReal
                     ("UPPER", "the upper limit, if any");
}


VaryCmd::VaryCmd(const std::string &name, VaryCmd *parent):
    Action(name, parent)
{}


VaryCmd::~VaryCmd()
{}


VaryCmd *VaryCmd::clone(const std::string &name) {
    return new VaryCmd(name, this);
}


void VaryCmd::execute() {
    // Fetch the command attributes.
    double pars[4];
    pars[0] = Attributes::getReal(itsAttr[NAME]);
    pars[1] = Attributes::getReal(itsAttr[STEP]);
    pars[2] = pars[3] = 0.0;

    // Set up the limits flag.
    int limits = NO_LIMIT;
    if(itsAttr[LOWER]) {
        pars[2] = Attributes::getReal(itsAttr[LOWER]);
        if(itsAttr[UPPER]) {
            pars[3] = Attributes::getReal(itsAttr[UPPER]);
            limits = BOTH_LIMITS;
        } else {
            limits = LOWER_LIMIT;
        }
    } else {
        if(itsAttr[UPPER]) {
            pars[3] = Attributes::getReal(itsAttr[UPPER]);
            limits = UPPER_LIMIT;
        } else {
            limits = NO_LIMIT;
        }
    }

    // Make the match variable name.
    std::string name(itsAttr[NAME].getImage());

    // Test for previously existing match variable;
    // if found, remove it.
    if(Match::block->findVariable(name)) {
        if(Options::warn) {
            std::cerr << "\n### Warning ### "
                      << "Parameter is already variable -- new attributes used.\n"
                      << std::endl;
        }

        Match::block->deleteVariable(name);
    }

    // Add new variable to the table.
    Match::block->addVariable(new ConcreteVar(name, itsAttr[NAME], limits, pars));
}
