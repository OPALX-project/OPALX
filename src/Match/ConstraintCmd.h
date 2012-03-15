#ifndef OPAL_ConstraintCmd_HH
#define OPAL_ConstraintCmd_HH 1

// ------------------------------------------------------------------------
// $RCSfile: ConstraintCmd.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConstraintCmd
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:10:01 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"
#include <iostream>
using namespace std;

// Class ConstraintCmd
// ------------------------------------------------------------------------
/// The CONSTRAINT command.

class ConstraintCmd: public Action {

public:

    /// Exemplar constructor.
    ConstraintCmd();

    virtual ~ConstraintCmd();

    /// Make clone.
    virtual ConstraintCmd *clone(const string &name);

    /// Execute the command.
    virtual void execute();

    /// Parse the command.
    virtual void parse(Statement &);

    /// Print the command.
    virtual void print(std::ostream &) const;

    /// Print help for the command.
    virtual void printHelp(ostream &) const;

private:

    // Not implemented.
    ConstraintCmd(const ConstraintCmd &);
    void operator=(const ConstraintCmd &);

    // Clone constructor.
    ConstraintCmd(const string &name, ConstraintCmd *parent);

    // The value of the relational operator:
    // 0 = "==", 1 = ">", 2 = "<".
    int relation;
};

#endif // OPAL_ConstraintCmd_HH
