#ifndef OPAL_PSystem_HH
#define OPAL_PSystem_HH

// ------------------------------------------------------------------------
// $RCSfile: PSystem.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: PSystem
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class PSystem
// ------------------------------------------------------------------------
/// The SYSTEM command.

class PSystem: public Action {

public:

    /// Exemplar constructor.
    PSystem();

    virtual ~PSystem();

    /// Make clone.
    virtual PSystem *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

    /// Parse command (special for one-attribute command).
    virtual void parse(Statement &);

private:

    // Not implemented.
    PSystem(const PSystem &);
    void operator=(const PSystem &);

    // Clone constructor.
    PSystem(const std::string &name, PSystem *parent);
};

#endif // OPAL_PSystem_HH