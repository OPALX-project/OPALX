#ifndef OPAL_ArithmeticError_HH
#define OPAL_ArithmeticError_HH

// ------------------------------------------------------------------------
// $RCSfile: ArithmeticError.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ArithmeticError
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/OpalException.h"

#include <string>

// Class ArithmeticError
// ------------------------------------------------------------------------
/// The abstract base class for all OPALX arithmetic exceptions.
//  An object derived from this class is thrown when an arithmetic error
//  occurs.

class ArithmeticError : public OpalException {
protected:
    /// The usual constructor.
    // Arguments:
    // [DL]
    // [DT][b]meth[/b]
    // [DD]the name of the method or function detecting the exception
    // [DT][b]msg [/b]
    // [DD]the message string identifying the exception
    // [/DL]
    ArithmeticError(const std::string& meth, const std::string& msg);

    ArithmeticError(const ArithmeticError&);
    virtual ~ArithmeticError();

private:
    // Not implemented.
    ArithmeticError();
};

#endif  // OPAL_ArithmeticError_HH
