// ------------------------------------------------------------------------
// $RCSfile: MatchParser.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MatchParser
//   The parser class for the OPAL match module.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/04/19 12:09:33 $
// $Author: opal $
//
//  JMJ: some changes
//
// ------------------------------------------------------------------------

#include "Match/MatchParser.h"
#include "AbstractObjects/OpalData.h"
#include "Match/ConstraintCmd.h"
#include "Match/MatchEnd.h"
#include "Match/MatchOption.h"
#include "Match/LMDif.h"
#include "Match/Migrad.h"
#include "Match/Simplex.h"
#include "Match/VaryCmd.h"


// Class MatchParser
// ------------------------------------------------------------------------


MatchParser::MatchParser():
    MatchDirectory() {
    // Matching commands.
    MatchDirectory.insert("CONSTRAINT", new ConstraintCmd());
    MatchDirectory.insert("ENDMATCH",   new MatchEnd());
    MatchDirectory.insert("LMDIF",      new LMDif());
    MatchDirectory.insert("MIGRAD",     new Migrad());
    MatchDirectory.insert("OPTION",     new MatchOption());
    MatchDirectory.insert("SIMPLEX",    new Simplex());
    MatchDirectory.insert("VARY",       new VaryCmd());

    // Table commands: Use exemplars from main directory.
    OpalData *OPAL = OpalData::getInstance();
    MatchDirectory.insert("SURVEY",     OPAL->find("SURVEY"));
    MatchDirectory.insert("TWISS",      OPAL->find("TWISS"));
    MatchDirectory.insert("TWISSTRACK", OPAL->find("TWISSTRACK"));
    MatchDirectory.insert("TWISSTRACK", OPAL->find("TWISSTRACK"));
    // JMJ 10/4/2000: added following commands here, seems harmless
    MatchDirectory.insert("VALUE", OPAL->find("VALUE"));
    MatchDirectory.insert("SYSTEM", OPAL->find("SYSTEM"));
    //   MatchDirectory.insert("OPTION", OPAL->find("OPTION")); // does not work
}


MatchParser::~MatchParser()
{}


Object *MatchParser::find(const string &name) const {
    return MatchDirectory.find(name);
}
