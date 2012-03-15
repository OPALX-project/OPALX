#ifndef OPAL_Stop_HH
#define OPAL_Stop_HH

// ------------------------------------------------------------------------
// $RCSfile: Stop.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Stop
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class Stop
// ------------------------------------------------------------------------
/// The STOP command.

class Stop: public Action {

public:

    /// Exemplar constructor.
    Stop();

    virtual ~Stop();

    /// Make clone.
    virtual Stop *clone(const string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    Stop(const Stop &);
    void operator=(const Stop &);

    // Clone constructor.
    Stop(const string &name, Stop *parent);
};

#endif // OPAL_Stop_HH
