// ------------------------------------------------------------------------
// $RCSfile: MatchCmd.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MatchCmd
//   The class for the OPAL MATCH command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Match/MatchCmd.h"
#include "Match/Match.h"


// Class MatchCmd
// ------------------------------------------------------------------------

MatchCmd::MatchCmd():
    Action(0, "MATCH",
           "The \"MATCH\" sub-command defines a match of values "
           "which can be matched.")
{
    AttributeHandler::addAttributeOwner("MATCH", AttributeHandler::COMMAND, "OPTION");
    AttributeHandler::addAttributeOwner("MATCH", AttributeHandler::COMMAND, "CONSTRAINT");
    AttributeHandler::addAttributeOwner("MATCH", AttributeHandler::COMMAND, "VARY");
    AttributeHandler::addAttributeOwner("MATCH", AttributeHandler::COMMAND, "MATCH");
    AttributeHandler::addAttributeOwner("MATCH", AttributeHandler::COMMAND, "LMDIF");
    AttributeHandler::addAttributeOwner("MATCH", AttributeHandler::COMMAND, "ENDMATCH");
    AttributeHandler::addAttributeOwner("MATCH", AttributeHandler::COMMAND, "MIGRAD");
    AttributeHandler::addAttributeOwner("MATCH", AttributeHandler::COMMAND, "SIMPLEX");
}


MatchCmd::MatchCmd(const std::string &name, MatchCmd *parent):
    Action(name, parent)
{}


MatchCmd::~MatchCmd()
{}


MatchCmd *MatchCmd::clone(const std::string &name) {
    return new MatchCmd(name, this);
}


void MatchCmd::execute() {
    // Execute match block.
    Match::block = new Match;
    Match::block->parser.run();

    // Clean up.
    delete Match::block;
    Match::block = 0;
}