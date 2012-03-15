// ------------------------------------------------------------------------
// $RCSfile: String.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: String
//   A class used to parse string attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/String.h"
#include "AbstractObjects/DoomReader.h"
#include "AbstractObjects/DoomWriter.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/SValue.h"
#include "Utilities/OpalException.h"
#include "Utilities/ParseError.h"

using namespace Expressions;


// Class String
// ------------------------------------------------------------------------

namespace Attributes {

    String::String(const string &name, const string &help):
        AttributeHandler(name, help, 0)
    {}


    String::~String()
    {}


    void String::doomGet
    (Attribute &attr, const DoomReader &reader, int index) const {
        attr.set(new SValue<string>(reader.getString(index)));
    }


    void String::doomPut
    (const Attribute &attr, DoomWriter &writer, int index) const {
        if(attr) {
            SValue<string> &value = dynamic_cast<SValue<string> &>(attr.getBase());
            writer.putString(index, value.evaluate());
            writer.putInt(index, 1);
        } else {
            writer.putInt(index, 0);
        }
    }


    const string &String::getType() const {
        static const string type("string");
        return type;
    }


    void String::parse(Attribute &attr, Statement &stat, bool) const {
        attr.set(new SValue<string>(parseString(stat, "String value expected.")));
    }

};
