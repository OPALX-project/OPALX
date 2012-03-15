// ------------------------------------------------------------------------
// $RCSfile: ErrorEnd.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorEnd
//   The class for the OPAL ENDERROR command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/ErrorEnd.h"
#include "Errors/Error.h"
#include "OpalParser/OpalParser.h"


// Class ErrorEnd
// ------------------------------------------------------------------------

ErrorEnd::ErrorEnd():
  Action(0, "ENDERROR",
	 "The \"ENDERROR\" sub-command ends error definition mode.")
{}


ErrorEnd::ErrorEnd(const string &name, ErrorEnd *parent):
  Action(name, parent)
{}


ErrorEnd::~ErrorEnd()
{}


ErrorEnd *ErrorEnd::clone(const string &name)
{
  return new ErrorEnd(name, this);
}


void ErrorEnd::execute()
{
  Error::block->parser.stop();
}
