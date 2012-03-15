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

OpalException::OpalException(const string &meth, const string &msg):
    ClassicException(meth, msg)
{}


OpalException::OpalException(const OpalException &rhs):
    ClassicException(rhs)
{}


OpalException::~OpalException()
{}
