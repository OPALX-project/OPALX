#ifndef OPAL_Error_HH
#define OPAL_Error_HH 1

// ------------------------------------------------------------------------
// $RCSfile: Error.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class Error
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/ErrorParser.h"

class Beamline;


// Class Error
// ------------------------------------------------------------------------
/// Error block.
//  This class encapsulates all data for Error definitions.
//  It acts as a communication area between the error definition commands.

class Error {

public:

    Error();
    ~Error();

    /// If true, errors are additive.
    bool addError;

    /// The line for which errors are to be generated.
    Beamline *itsLine;

    /// The parser active during error definition.
    ErrorParser parser;

    /// The static block for the error data.
    static Error *block;
};

#endif // OPAL_Error_HH
