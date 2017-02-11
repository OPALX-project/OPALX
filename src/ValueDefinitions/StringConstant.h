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
    virtual StringConstant *clone(const std::string &name);

    /// Print the constant.
    virtual void print(std::ostream &) const;

    /// Print its value
    virtual void printValue(std::ostream &os) const;

    /// Return value.
    virtual std::string getString() const;

private:

    // Not implemented.
    StringConstant(const StringConstant &);
    void operator=(const StringConstant &);

    // Clone constructor.
    StringConstant(const std::string &name, StringConstant *parent);
    StringConstant(const std::string &name, StringConstant *parent, const std::string &value);
};

#endif // OPAL_StringConstant_HH