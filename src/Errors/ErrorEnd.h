#ifndef OPAL_ErrorEnd_HH
#define OPAL_ErrorEnd_HH

// ------------------------------------------------------------------------
// $RCSfile: ErrorEnd.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorEnd
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"

class ErrorParser;


// Class ErrorEnd
// ------------------------------------------------------------------------
/// The ENDERROR command.

class ErrorEnd: public Action {

public:

    /// Exemplar constructor.
    ErrorEnd();

    virtual ~ErrorEnd();

    /// Make clone.
    virtual ErrorEnd *clone(const string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    ErrorEnd(const ErrorEnd &);
    void operator=(const ErrorEnd &);

    // Clone constructor.
    ErrorEnd(const string &name, ErrorEnd *parent);
};

#endif // OPAL_ErrorEnd_HH
