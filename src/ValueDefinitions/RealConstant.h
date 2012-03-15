#ifndef OPAL_RealConstant_HH
#define OPAL_RealConstant_HH
// ------------------------------------------------------------------------
// $RCSfile: RealConstant.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RealConstant
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:49 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/ValueDefinition.h"


// Class RealConstant
// ------------------------------------------------------------------------
/// The REAL CONSTANT definition.

class RealConstant: public ValueDefinition {

public:

    /// Exemplar constructor.
    RealConstant();

    /// Constructor for built-in constants.
    RealConstant(const string &name, RealConstant *parent, double value);

    virtual ~RealConstant();

    /// Test if object can be replaced.
    //  Always false for constants.
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual RealConstant *clone(const string &name);

    /// Read constant from DOOM data base.
    virtual void doomGet(const DoomReader &);

    /// Write constant to DOOM data base.
    virtual void doomPut(DoomWriter &) const;

    /// Print the constant.
    virtual void print(std::ostream &) const;

    /// Return value.
    virtual double getReal() const;

private:

    // Not implemented.
    RealConstant(const RealConstant &);
    void operator=(const RealConstant &);

    // Clone constructor.
    RealConstant(const string &name, RealConstant *parent);
};

#endif // OPAL_RealConstant_HH
