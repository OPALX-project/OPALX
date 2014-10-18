#ifndef OPAL_RealVariable_HH
#define OPAL_RealVariable_HH
// ------------------------------------------------------------------------
// $RCSfile: RealVariable.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RealVariable
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:49 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/ValueDefinition.h"

// Class RealVariable
// ------------------------------------------------------------------------
/// The REAL VARIABLE definition.

class RealVariable: public ValueDefinition {

public:

    /// Exemplar constructor.
    RealVariable();

    /// Constructor for built-in variables.
    RealVariable(const std::string &name, RealVariable *parent, double value);

    virtual ~RealVariable();

    /// Test for allowed replacement.
    //  True, if [b]rhs[/b] is a real variable.
    virtual bool canReplaceBy(Object *rhs);

    /// Make clone.
    virtual RealVariable *clone(const std::string &name);

    /// Print the variable.
    virtual void print(std::ostream &) const;

    /// Return value.
    virtual double getReal() const;

private:

    // Not implemented.
    RealVariable(const RealVariable &);
    void operator=(const RealVariable &);

    // Clone constructor.
    RealVariable(const std::string &name, RealVariable *parent);
};

#endif // OPAL_RealVariable_HH
