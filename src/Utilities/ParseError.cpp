// ------------------------------------------------------------------------
// $RCSfile: ParseError.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParseError
//   Parse error in the OPALX parser.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/ParseError.h"

// Class ParseError
// ------------------------------------------------------------------------

ParseError::ParseError(const std::string& meth, const std::string& msg)
    : OpalException(meth, msg) {}

ParseError::ParseError(const ParseError& rhs) : OpalException(rhs) {}

ParseError::~ParseError() {}
