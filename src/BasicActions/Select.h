#ifndef OPAL_Select_HH
#define OPAL_Select_HH 1

// ------------------------------------------------------------------------
// $RCSfile: Select.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Select
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"

class Beamline;


// Class Select
// ------------------------------------------------------------------------
/// The SELECT command.

class Select: public Action {

public:

    /// Exemplar constructor.
    Select();

    virtual ~Select();

    /// Make clone.
    virtual Select *clone(const string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    Select(const Select &);
    void operator=(const Select &);

    // Clone constructor.
    Select(const string &name, Select *parent);

    // Do the selection.
    void select(const Beamline &);
};

#endif // OPAL_Select_H
