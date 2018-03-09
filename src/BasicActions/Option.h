#ifndef OPAL_Option_HH
#define OPAL_Option_HH

// ------------------------------------------------------------------------
// $RCSfile: Option.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Option
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"

// Class Option
// ------------------------------------------------------------------------
/// The OPTION command.
//  The user interface allowing setting of OPAL options.
//  The actual option flags are contained in namespace Options.

class Option: public Action {

public:

    /// Exemplar constructor.
    Option();

    virtual ~Option();

    /// Make clone.
    virtual Option *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:
    void handlePsDumpFrame(const std::string &dumpFrame);

    // Not implemented.
    Option(const Option &);
    void operator=(const Option &);

    // Clone constructor.
    Option(const std::string &name, Option *parent);
};

#endif // OPAL_Option_HH