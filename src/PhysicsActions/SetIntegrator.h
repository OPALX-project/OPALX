#ifndef OPAL_SetIntegrator_HH
#define OPAL_SetIntegrator_HH

// ------------------------------------------------------------------------
// $RCSfile: SetIntegrator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SetIntegrator
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class SetIntegrator
// ------------------------------------------------------------------------
/// The SETINT command.

class SetIntegrator: public Action {

public:

    /// Exemplar constructor.
    SetIntegrator();

    virtual ~SetIntegrator();

    /// Make clone.
    virtual SetIntegrator *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    SetIntegrator(const SetIntegrator &);
    void operator=(const SetIntegrator &);

    // Clone constructor.
    SetIntegrator(const std::string &name, SetIntegrator *parent);
};

#endif // __SetIntegrator_HH
