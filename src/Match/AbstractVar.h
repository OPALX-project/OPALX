#ifndef OPAL_AbstractVar_HH
#define OPAL_AbstractVar_HH 1

// ------------------------------------------------------------------------
// $RCSfile: AbstractVar.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AbstractVar
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include <iosfwd>
#include <string>


// Class AbstractVar
// ------------------------------------------------------------------------
/// Abstract base for a matching variable.
//  The interface allows for variable transformations used to implement
//  upper and/or lower limits for the values.

class AbstractVar {

public:

    /// Constructor.
    //  Assign the variable name.
    AbstractVar(const std::string &name);

    virtual ~AbstractVar();

    /// Get the variable name.
    virtual const std::string &getName() const;

    /// Get the current internal parameter value.
    virtual double getInternalValue() const = 0;

    /// Set the current internal parameter value.
    virtual void setInternalValue(double) = 0;

    /// Get the current external parameter value.
    virtual double getExternalValue() const = 0;

    /// Set the current external parameter value.
    virtual void setExternalValue(double) = 0;

    /// Print the variable name and value.
    virtual void print(std::ostream &) const = 0;

protected:

    /// Name of the variable.
    const std::string itsName;

private:

    // Not implemented.
    AbstractVar();
    AbstractVar(const AbstractVar &);
    void operator=(const AbstractVar &);
};

#endif // OPAL_AbstractVar_HH
