#ifndef OPAL_ErrorPrint_HH
#define OPAL_ErrorPrint_HH

// ------------------------------------------------------------------------
// $RCSfile: ErrorPrint.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorPrint
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"
#include "Errors/AlignBase.h"
#include "Errors/MPBase.h"
#include <fstream>

class AlignWrapper;


// Class ErrorPrint
// ------------------------------------------------------------------------
/// The EPRINT command.

class ErrorPrint:
    public Action, private AlignBase, private MPBase {

public:

    /// Exemplar constructor.
    ErrorPrint();

    virtual ~ErrorPrint();

    /// Make clone.
    virtual ErrorPrint *clone(const std::string &name);

    // Execute the command.
    virtual void execute();

    /// Print error field.
    //  Use the name and occurrence count for the current element.  Relative
    //  errors are related to [b]mainField[/b]. The field error is fetched
    //  from [b]errorField[/b].  This class's execute() method uses a MPHandler
    //  to call this method for all elements which have a field error.
    virtual void fieldError(const std::string &name, int occur,
                            const BMultipoleField &mainField,
                            BMultipoleField &errorField);

    /// Print misalignement.
    //  All information is taken from the given AlignWrapper.
    //  This class's execute() method uses an AlignHandler to call this
    //  method for all elements which have a misalignment.
    virtual void misalignment(const AlignWrapper &, int occur);

private:

    // Not implemented.
    ErrorPrint(const ErrorPrint &);
    void operator=(const ErrorPrint &);

    // Clone constructor.
    ErrorPrint(const std::string &name, ErrorPrint *parent);

    // Output stream; "mutable" is required to enable writing.
    mutable std::ofstream os;
};

#endif // OPAL_ErrorPrint_HH
