// ------------------------------------------------------------------------
// $RCSfile: OpalException.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalException
//   The base class for all OPAL exceptions.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:48 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Utilities/OpalException.h"

// Class OpalException
// ------------------------------------------------------------------------

OpalException::OpalException(const std::string& meth, const std::string& msg)
    : message(msg), method(meth) {}

OpalException::OpalException(const OpalException& rhs) : message(rhs.message), method(rhs.method) {}

OpalException::~OpalException() {}

const std::string& OpalException::what() const { return message; }

const std::string& OpalException::where() const { return method; }
