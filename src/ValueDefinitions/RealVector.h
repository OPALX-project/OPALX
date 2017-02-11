#ifndef OPAL_RealVector_HH
#define OPAL_RealVector_HH
// ------------------------------------------------------------------------
// $RCSfile: RealVector.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RealVector
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:49 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/ValueDefinition.h"


// Class RealVector
// ------------------------------------------------------------------------
/// The REAL VECTOR definition.

class RealVector: public ValueDefinition {

public:

    /// Exemplar constructor.
    RealVector();

    virtual ~RealVector();

    /// Test for allowed replacement.
    //  True, if [b]rhs[/b] is a real vector.
    virtual bool canReplaceBy(Object *rhs);

    /// Make clone.
    virtual RealVector *clone(const std::string &name);

    /// Print the vector.
    virtual void print(std::ostream &) const;

    /// Print its value
    virtual void printValue(std::ostream &os) const;

    /// Return indexed value.
    virtual double getRealComponent(int) const;

private:

    // Not implemented.
    RealVector(const RealVector &);
    void operator=(const RealVector &);

    // Clone constructor.
    RealVector(const std::string &name, RealVector *parent);
};

#endif // OPAL_RealVector_HH