#ifndef OPAL_Call_HH
#define OPAL_Call_HH

// ------------------------------------------------------------------------
// $RCSfile: Call.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Call
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class Call
// ------------------------------------------------------------------------
/// The CALL command.

class Call: public Action {

public:

    /// Exemplar constructor.
    Call();

    virtual ~Call();

    /// Make clone.
    virtual Call *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

    /// Parse command (special for one-attribute command).
    virtual void parse(Statement &);

private:

    // Not implemented.
    Call(const Call &);
    void operator=(const Call &);

    // Clone constructor.
    Call(const std::string &name, Call *parent);
};

#endif // OPAL_Call_HH
