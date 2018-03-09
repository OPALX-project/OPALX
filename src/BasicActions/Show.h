#ifndef OPAL_Show_HH
#define OPAL_Show_HH

// ------------------------------------------------------------------------
// $RCSfile: Show.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Show
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class Show
// -----------------------------------------------------------------------
/// The SHOW command.

class Show: public Action {

public:

    /// Exemplar constructor.
    Show();

    virtual ~Show();

    /// Make clone.
    virtual Show *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

    /// Parse command (special for one-attribute command).
    virtual void parse(Statement &);

private:

    // Not implemented.
    Show(const Show &);
    void operator=(const Show &);

    // Clone constructor.
    Show(const std::string &name, Show *parent);
};

#endif // OPAL_Show_HH
