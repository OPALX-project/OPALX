#ifndef OPAL_BoolConstant_HH
#define OPAL_BoolConstant_HH
// ------------------------------------------------------------------------
// $RCSfile: BoolConstant.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BoolConstant
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:49 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/ValueDefinition.h"


// Class BoolConstant
// ------------------------------------------------------------------------
/// The BOOL CONSTANT definition.

class BoolConstant: public ValueDefinition {

public:

    /// Exemplar constructor.
    BoolConstant();

    virtual ~BoolConstant();

    /// Test if object can be replaced.
    //  Always false for constants.
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual BoolConstant *clone(const std::string &name);

    /// Print the constant.
    virtual void print(std::ostream &) const;

    /// Return value.
    virtual bool getBool() const;

private:

    // Not implemented.
    BoolConstant(const BoolConstant &);
    void operator=(const BoolConstant &);

    // Clone constructor.
    BoolConstant(const std::string &name, BoolConstant *parent);
};

#endif // OPAL_BoolConstant_HH
