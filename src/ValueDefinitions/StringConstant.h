#ifndef OPAL_StringConstant_HH
#define OPAL_StringConstant_HH
// ------------------------------------------------------------------------
// $RCSfile: StringConstant.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: StringConstant
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:49 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/ValueDefinition.h"


// Class StringConstant
// ------------------------------------------------------------------------
/// The STRING CONSTANT definition.

class StringConstant: public ValueDefinition {

public:

    /// Exemplar constructor.
    StringConstant();

    virtual ~StringConstant();

    /// Test if object can be replaced.
    //  True, if [b]rhs[/b] is a string constant.
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual StringConstant *clone(const string &name);

    /// Read the constant from the DOOM data base.
    virtual void doomGet(const DoomReader &);

    /// Write the constant to the DOOM data base.
    virtual void doomPut(DoomWriter &) const;

    /// Print the constant.
    virtual void print(std::ostream &) const;

    /// Return value.
    virtual string getString() const;

private:

    // Not implemented.
    StringConstant(const StringConstant &);
    void operator=(const StringConstant &);

    // Clone constructor.
    StringConstant(const string &name, StringConstant *parent);
};

#endif // OPAL_StringConstant_HH
