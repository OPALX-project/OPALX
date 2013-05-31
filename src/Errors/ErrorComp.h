#ifndef OPAL_ErrorComp_HH
#define OPAL_ErrorComp_HH

// ------------------------------------------------------------------------
// $RCSfile: ErrorComp.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorComp
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"
#include "Errors/MPBase.h"

class BMultipoleField;


// Class ErrorComp
// ------------------------------------------------------------------------
/// The EFCOMP command.

class ErrorComp: public Action, private MPBase {

public:

    /// Exemplar constructor.
    ErrorComp();

    virtual ~ErrorComp();

    /// Make clone.
    virtual ErrorComp *clone(const string &name);

    /// Execute the command.
    virtual void execute();

    /// Generate error field.
    //  The first two arguments are not used, relative errors are related
    //  to [b]mainField[/b]. The result is stored in [b]errorField[/b].
    //  This class's execute() method uses a MPHandler to call this function
    //  for each element which may have a multipole field error.
    virtual void fieldError(const string &, int,
                            const BMultipoleField &mainField,
                            BMultipoleField &errorField);

private:

    // Not implemented.
    ErrorComp(const ErrorComp &);
    void operator=(const ErrorComp &);

    // Clone constructor.
    ErrorComp(const string &name, ErrorComp *parent);

    // Counters for new and existing errors.
    int oldError, newError;
};

#endif // OPAL_ErrorComp_HH
