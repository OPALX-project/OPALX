#ifndef OPAL_VaryCmd_HH
#define OPAL_VaryCmd_HH 1

// ------------------------------------------------------------------------
// $RCSfile: VaryCmd.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: VaryCmd
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class VaryCmd
// ------------------------------------------------------------------------
/// The VARY command.

class VaryCmd: public Action {

public:

    /// Exemplar constructor.
    VaryCmd();

    virtual ~VaryCmd();

    /// Make clone.
    virtual VaryCmd *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    VaryCmd(const VaryCmd &);
    void operator=(const VaryCmd &);

    // Clone constructor.
    VaryCmd(const std::string &name, VaryCmd *parent);
};

#endif // OPAL_VaryCmd_HH
