#ifndef OPAL_DoomCmd_HH
#define OPAL_DoomCmd_HH

// ------------------------------------------------------------------------
// $RCSfile: DoomCmd.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DoomCmd
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class DoomCmd
// ------------------------------------------------------------------------
/// The DOOM command.

class DoomCmd: public Action {

public:

    /// Exemplar constructor.
    DoomCmd();

    virtual ~DoomCmd();

    /// Make clone.
    virtual DoomCmd *clone(const string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    DoomCmd(const DoomCmd &);
    void operator=(const DoomCmd &);

    // Clone constructor.
    DoomCmd(const string &name, DoomCmd *parent);
};

#endif // OPAL_DoomCmd_HH
