// ------------------------------------------------------------------------
// $RCSfile: DoomCmd.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DoomCmd
//   Defines the OPAL "DOOM" command.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/14 06:45:26 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "BasicActions/DoomCmd.h"
#include "AbstractObjects/DoomDB.h"
#include "Attributes/Attributes.h"
#include "Utilities/Options.h"
#include <iostream>


// Class DoomCmd
// ------------------------------------------------------------------------

// The attributes of class DoomCmd.
namespace {
    enum {
        OPEN,
        CLOSE,
        SHUT,
        NOUPDATE,
        DOOM_DEBUG,
        SIZE
    };
}

//JMJ 3/8/2000 grammatical corrections to following message text

DoomCmd::DoomCmd():
    Action(SIZE, "DOOM",
           "The \"DOOM\" statement allows actions on the DOOM "
           "data base.\n"
           "It can also be used to set debugging options."
           "If both \"OPEN\" and\"CLOSE\" are entered, \"CLOSE\" is ignored.") {
    itsAttr[OPEN] = Attributes::makeBool
                    ("OPEN",
                     "If true, the data base is re-opened. "
                     "All objects changed in the data base are re-read.");
    itsAttr[CLOSE] = Attributes::makeBool
                     ("CLOSE",
                      "If true, the data base is closed. "
                      "All objects changed in memory are re-rewritten.");
    itsAttr[SHUT] = Attributes::makeBool
                    ("SHUT",
                     "Same as \"CLOSE\".");
    itsAttr[NOUPDATE] = Attributes::makeBool
                        ("NOUPDATE",
                         "If true, the nothing will be saved when the OPAL program is shut down.");
    itsAttr[DOOM_DEBUG] = Attributes::makeReal
                          ("DEBUG",
                           "Allows to set debugging options:\n"
                           "\t0: no print\n"
                           "\t1: terse debug print\n"
                           "\t2: full debug print",
                           0.0);
}


DoomCmd::DoomCmd(const string &name, DoomCmd *parent):
    Action(name, parent)
{}


DoomCmd::~DoomCmd()
{}


DoomCmd *DoomCmd::clone(const string &name) {
    return new DoomCmd(name, this);
}


void DoomCmd::execute() {
    bool open = Attributes::getBool(itsAttr[OPEN]);
    bool close =
        Attributes::getBool(itsAttr[CLOSE]) || Attributes::getBool(itsAttr[SHUT]);

    if(open && close) {
        if(Options::warn) {
            std::cerr << std::endl
                      << "### Warning ### "
                      << "Conflicting options \"OPEN\" and \"CLOSE\" ignored.\n"
                      << std::endl;
        }
    } else if(open) {
        // Re-open data base.
        DOOM_DB.reOpen();
    } else if(close) {
        // Close the data base.
        DOOM_DB.shut();
    }

    // Set no-update option.
    if(itsAttr[NOUPDATE]) {
        DOOM_DB.setUpdate(! Attributes::getBool(itsAttr[NOUPDATE]));
    }

    // Set debug option.
    if(itsAttr[DOOM_DEBUG]) {
        DOOM_DB.setDebug(int(Attributes::getReal(itsAttr[DOOM_DEBUG])));
    }
}
