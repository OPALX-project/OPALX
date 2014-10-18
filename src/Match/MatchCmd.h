#ifndef OPAL_MatchCmd_HH
#define OPAL_MatchCmd_HH 1

// ------------------------------------------------------------------------
// $RCSfile: MatchCmd.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MatchCmd
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class MatchCmd
// ------------------------------------------------------------------------
/// The MATCH command.

class MatchCmd: public Action {

public:

    /// Exemplar constructor.
    MatchCmd();

    virtual ~MatchCmd();

    /// Make clone.
    virtual MatchCmd *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    MatchCmd(const MatchCmd &);
    void operator=(const MatchCmd &);

    // Clone constructor.
    MatchCmd(const std::string &name, MatchCmd *parent);
};

#endif // OPAL_MatchCmd_HH
