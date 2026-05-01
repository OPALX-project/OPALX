// ------------------------------------------------------------------------
// $RCSfile: FormatError.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: FormatError
//   Format error in Tps/Vps read.
//   Format error in Tps/Vps read.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/FormatError.h"


// Class FormatError
// ------------------------------------------------------------------------


FormatError::FormatError(const std::string &meth, const std::string &msg):
    OpalException(meth, msg)
{}


FormatError::FormatError(const FormatError &rhs):
    OpalException(rhs)
{}


FormatError::~FormatError()
{}
