#ifndef OPAL_ConcreteVar_HH
#define OPAL_ConcreteVar_HH 1

// ------------------------------------------------------------------------
// $RCSfile: ConcreteVar.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConcreteVar
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Match/AbstractVar.h"
#include "AbstractObjects/Attribute.h"
#include "Match/MatchLimits.h"


// Class ConcreteVar
// ------------------------------------------------------------------------
/// Concrete class for a matching variable.
//  Implements the setting and retrieving of the value in the system
//  to be adjusted.  Upper and lower limits are implemented by variable
//  transformations which map the infinite range to the limited range.

class ConcreteVar: public AbstractVar {

public:

    /// Constructor.
    //  Uses the following arguments:
    //  [ol]
    //  [li] The variable name.
    //  [li] An attribute containing a reference to thevalue to be adjusted.
    //  [li] A code for limit values: 0=no limit, 1=lower, 2=upper, 3=both.
    //  [li] Value, lower, upper, and step.
    //  [/ol]
    ConcreteVar(const string &name, Attribute &data,
                int limits, double pars[4]);

    virtual ~ConcreteVar();

    /// Get the current internal parameter value.
    //  The internal value is unlimited, it maps to the external value
    //  so as to keep the latter constrained.
    virtual double getInternalValue() const;

    /// Set the current internal parameter value.
    //  The internal value is unlimited, it maps to the external value
    //  so as to keep the latter constrained.
    virtual void setInternalValue(double);

    /// Get the current external parameter value.
    //  The external value should be consistent with the given limits.
    virtual double getExternalValue() const;

    /// Set the current external parameter value.
    //  The external value should be consistent with the given limits.
    virtual void setExternalValue(double);

    /// Print the variable name and value.
    virtual void print(std::ostream &) const;

private:

    // Not implemented.
    ConcreteVar();
    ConcreteVar(const ConcreteVar &);
    void operator=(const ConcreteVar &);

    // The attribute to be manipulated.
    Attribute itsAttr;

    // The current internal parameter step.
    double itsStep;

    // The fixed parameter limits.
    int itsLimits;
    double itsMin;
    double itsMax;
};

#endif // OPAL_ConcreteVar_HH
