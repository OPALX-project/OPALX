#ifndef CLASSIC_EigenvalueError_HH
#define CLASSIC_EigenvalueError_HH

// ------------------------------------------------------------------------
// $RCSfile: EigenvalueError.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EigenvalueError
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


// Class EigenvalueError
// ------------------------------------------------------------------------
/// Eigenvalue error exception.
//  This exception is thrown, when an eigenvalue routine fails to compute
//  all eigenvalues.

class EigenvalueError: public ArithmeticError {

public:

    /// The usual constructor.
    // Arguments:
    // [DL]
    // [DT][b]meth[/b]
    // [DD]the name of the method or function detecting the exception
    // [DT][b]msg [/b]
    // [DD]the message string identifying the exception
    // [/DL]
    EigenvalueError(const string &meth, const string &msg);

    EigenvalueError(const EigenvalueError &);
    virtual ~EigenvalueError();

private:

    // Not implemented.
    EigenvalueError();
};

#endif // CLASSIC_EigenvalueError_HH
