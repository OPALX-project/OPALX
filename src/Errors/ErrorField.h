#ifndef OPAL_ErrorField_HH
#define OPAL_ErrorField_HH

// ------------------------------------------------------------------------
// $RCSfile: ErrorField.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorField
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


// Class ErrorField
// ------------------------------------------------------------------------
/// The EFIELD command.

class ErrorField: public Action, private MPBase {

public:

    /// Exemplar constructor.
    ErrorField();

    virtual ~ErrorField();

    /// Make clone.
    virtual ErrorField *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

    /// Generate error field.
    //  The first two arguments are not used, relative errors are related
    //  to [b]mainField[/b]. The result is stored in [b]errorField[/b].
    //  This class's execute() method uses a MPHandler to call this function
    //  for each element which may have a multipole field error.
    virtual void fieldError(const std::string &name, int occur,
                            const BMultipoleField &mainField,
                            BMultipoleField &errorField);

private:

    // Not implemented.
    ErrorField(const ErrorField &);
    void operator=(const ErrorField &);

    // Clone constructor.
    ErrorField(const std::string &name, ErrorField *parent);

    // Counters for new and existing errors.
    int oldError, newError;
};

#endif // OPAL_ErrorField_HH
