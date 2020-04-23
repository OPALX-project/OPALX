// ------------------------------------------------------------------------
// $RCSfile: UString.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: UString
//   A class used to parse string attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/UString.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/SValue.h"
#include "Utilities/OpalException.h"
#include "Utilities/ParseError.h"
#include "Utilities/Util.h"

using namespace Expressions;


// Class UString
// ------------------------------------------------------------------------

namespace Attributes {

    UString::UString(const std::string &name, const std::string &help):
        AttributeHandler(name, help, 0)
    {}


    UString::~UString()
    {}


    const std::string &UString::getType() const {
        static const std::string type("string");
        return type;
    }


    void UString::parse(Attribute &attr, Statement &stat, bool) const {
        attr.set(new SValue<std::string>(Util::toUpper(parseString(stat, "UString value expected."))));
    }

};