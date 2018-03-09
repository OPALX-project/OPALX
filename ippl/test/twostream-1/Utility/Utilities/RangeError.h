#ifndef CLASSIC_RangeError_HH
#define CLASSIC_RangeError_HH

// ------------------------------------------------------------------------
// $RCSfile: RangeError.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RangeError
//   
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/ArithmeticError.h"


// Class RangeError
// ------------------------------------------------------------------------
//: Range error.
//  This exception is thrown, when a CLASSIC routine or method detects an
//  index out of range.

class RangeError: public ArithmeticError {

public:

  //: The usual constructor.
  // Arguments:
  // [DL]
  // [DT][b]meth[/b]
  // [DD]the name of the method or function detecting the exception
  // [DT][b]msg [/b]
  // [DD]the message string identifying the exception
  // [/DL]
  // Construction/destruction.
  RangeError(const string &meth, const string &msg);

  RangeError(const RangeError &);
  virtual ~RangeError();

private:

  // Not implemented.
  RangeError();
};

#endif // CLASSIC_RangeError_HH
