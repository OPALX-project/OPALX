#ifndef OPAL_ConcreteFun_HH
#define OPAL_ConcreteFun_HH 1

// ------------------------------------------------------------------------
// $RCSfile: ConcreteFun.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConcreteFun
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Match/AbstractFun.h"
#include "AbstractObjects/Attribute.h"
#include <vector>


// Class ConcreteFun
// ------------------------------------------------------------------------
/// A single matching constraints or an array of matching constraints.
//  The left and right hand sides and the weights of the constraint are
//  all given by array expressions.

class ConcreteFun: public AbstractFun {

public:

    /// Constructor.
    //  Uses the following arguments:
    //  [ol]
    //  [li] The left-hand side(s) for the constraint.
    //  [li] A code for the type of constraint.
    //  [li] The right-hand side(s) for the constraint.
    //  [li] The weight(s) for the constraint.
    //  [/ol]
    ConcreteFun(Attribute &lhs, int rel, Attribute &rhs, Attribute &wgt);

    virtual ~ConcreteFun();

    /// Get the number of constrained values.
    virtual int countConstraints() const;

    /// Evaluate the matching function(s).
    //  Increment [b]n[/b] for each constrained value and store the value
    //  in vector [b]f[/b].
    virtual void evaluate(Vector<double> &f, int &) const;

    /// Print the function name and value(s).
    virtual void print(std::ostream &) const;

private:

    // Not implemented.
    ConcreteFun();
    ConcreteFun(const ConcreteFun &);
    void operator=(const ConcreteFun &);

    // The constraint expression to be adjusted and its weight.
    Attribute itsLhs;
    Attribute itsRhs;
    int relation;
    std::vector<double> itsWeight;

    // The value cache.
    mutable std::vector<double> lValue;
    mutable std::vector<double> rValue;
    mutable std::vector<double> value;
};

#endif // OPAL_ConcreteFun_HH
