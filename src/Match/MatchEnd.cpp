// ------------------------------------------------------------------------
// $RCSfile: MatchEnd.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MatchEnd
//   The class for the OPAL ENDMATCH command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Match/MatchEnd.h"
#include "Match/Match.h"
#include "Match/MatchParser.h"
#include "Match/MatchState.h"


// Class MatchEnd
// ------------------------------------------------------------------------

MatchEnd::MatchEnd():
  Action(0, "ENDMATCH",
	 "The \"ENDMATCH\" sub-command stops Matching.")
{}


MatchEnd::MatchEnd(const string &name, MatchEnd *parent):
  Action(name, parent)
{}


MatchEnd::~MatchEnd()
{}


MatchEnd *MatchEnd::clone(const string &name)
{
  return new MatchEnd(name, this);
}


void MatchEnd::execute()
{
  Match::block->print("ENDMATCH", TERMINATED);
  Match::block->parser.stop();
}
