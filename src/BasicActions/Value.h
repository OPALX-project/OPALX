#ifndef OPAL_Value_HH
#define OPAL_Value_HH

// ------------------------------------------------------------------------
// $RCSfile: Value.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Value
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class Value
// ------------------------------------------------------------------------
/// The VALUE command.

class Value: public Action {

public:

    /// Exemplar constructor.
    Value();

    virtual ~Value();

    /// Make clone.
    virtual Value *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

    /// Parse command (special for one-attribute command).
    virtual void parse(Statement &);

private:

    // Not implemented.
    Value(const Value &);
    void operator=(const Value &);

    // Clone constructor.
    Value(const std::string &name, Value *parent);
};

#endif // OPAL_Value_HH
