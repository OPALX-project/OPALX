#ifndef OPAL_LogicalError_HH
#define OPAL_LogicalError_HH

// ------------------------------------------------------------------------
// $RCSfile: LogicalError.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: LogicalError
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/OpalException.h"

#include <string>

// Class LogicalError
// ------------------------------------------------------------------------
/// Logical error exception.
//  This exception is thrown when OPALX detects an inconsistent call to
//  a routine or method.

class LogicalError : public OpalException {
public:
    /// The usual constructor.
    // Arguments:
    // [DL]
    // [DT][b]meth[/b]
    // [DD]the name of the method or function detecting the exception
    // [DT][b]msg [/b]
    // [DD]the message string identifying the exception
    // [/DL]
    LogicalError(const std::string& meth, const std::string& msg);

    LogicalError(const LogicalError&);
    virtual ~LogicalError();

private:
    // Not implemented.
    LogicalError();
};

#endif  // OPAL_LogicalError_HH
