#ifndef OPAL_ErrorAlign_HH
#define OPAL_ErrorAlign_HH

// ------------------------------------------------------------------------
// $RCSfile: ErrorAlign.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorAlign
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"
#include "Errors/AlignBase.h"


// Class ErrorAlign
// ------------------------------------------------------------------------
/// The EALIGN command.

class ErrorAlign: public Action, private AlignBase {

public:

    /// Exemplar constructor.
    ErrorAlign();

    virtual ~ErrorAlign();

    /// Make clone.
    virtual ErrorAlign *clone(const string &name);

    /// Execute the command.
    virtual void execute();

    /// Generate a misalignment.
    //  Prepare a misalignment error according to the specification in this
    //  command and assign it to the given AlignWrapper.  This class's
    //  execute() method uses an AlignHandler to call this function for
    //  each AlignWrapper in a beam line.
    virtual void misalignment(const AlignWrapper &, int occur);

private:

    // Not implemented.
    ErrorAlign(const ErrorAlign &);
    void operator=(const ErrorAlign &);

    // Clone constructor.
    ErrorAlign(const string &name, ErrorAlign *parent);

    // Counters for new and existing errors.
    int oldError, newError;
};

#endif // OPAL_ErrorAlign_HH
