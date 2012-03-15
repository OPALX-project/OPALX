// ------------------------------------------------------------------------
// $RCSfile: RangeError.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RangeError
//   Index out of range in matrix or vector operation.
//   
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/RangeError.h"


// Class RangeError
// ------------------------------------------------------------------------


RangeError::RangeError(const string &meth, const string &msg):
  ArithmeticError(meth, msg)
{}


RangeError::RangeError(const RangeError &rhs):
  ArithmeticError(rhs)
{}


RangeError::~RangeError()
{}
