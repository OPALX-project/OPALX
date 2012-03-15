#ifndef OPAL_Static_HH
#define OPAL_Static_HH

// ------------------------------------------------------------------------
// $RDSfile: Static.hh,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Static
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class Static
// ------------------------------------------------------------------------
/// The STATIC command.

class Static: public Action {

public:

    /// Exemplar constructor.
    Static();

    virtual ~Static();

    /// Make clone.
    virtual Static *clone(const string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    Static(const Static &);
    void operator=(const Static &);

    // Clone constructor.
    Static(const string &name, Static *parent);
};

#endif // OPAL_Static_HH
