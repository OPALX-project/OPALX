#ifndef OPAL_ErrorSelect_HH
#define OPAL_ErrorSelect_HH

// ------------------------------------------------------------------------
// $RCSfile: ErrorSelect.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorSelect
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class ErrorSelect
// ------------------------------------------------------------------------
/// The error SELECT command.

class ErrorSelect: public Action {

public:

    /// Exemplar constructor.
    ErrorSelect();

    virtual ~ErrorSelect();

    /// Make clone.
    virtual ErrorSelect *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    ErrorSelect(const ErrorSelect &);
    void operator=(const ErrorSelect &);

    // Clone constructor.
    ErrorSelect(const std::string &name, ErrorSelect *parent);
};

#endif // OPAL_ErrorSelect_HH
