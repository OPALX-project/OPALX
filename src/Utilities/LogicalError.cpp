// ------------------------------------------------------------------------
// $RCSfile: LogicalError.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: LogicalError
//   The class for all OPAL exceptions related to object Logicals.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/LogicalError.h"


// Class LogicalError
// ------------------------------------------------------------------------


LogicalError::LogicalError(const std::string &meth, const std::string &msg):
    OpalException(meth, msg)
{}


LogicalError::LogicalError(const LogicalError &rhs):
    OpalException(rhs)
{}


LogicalError::~LogicalError()
{}
