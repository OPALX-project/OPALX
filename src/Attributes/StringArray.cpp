// ------------------------------------------------------------------------
// $RCSfile: StringArray.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class StringArray:
//   A class used to parse array attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/StringArray.h"
#include "Attributes/Attributes.h"
#include "AbstractObjects/DoomReader.h"
#include "AbstractObjects/DoomWriter.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/AValue.h"
#include "Utilities/OpalException.h"
#include "Utilities/ParseError.h"
#include <vector>

using namespace Expressions;


// Class StringArray
// ------------------------------------------------------------------------

namespace Attributes {

    StringArray::StringArray(const string &name, const string &help):
        AttributeHandler(name, help, 0)
    {}


    StringArray::~StringArray()
    {}


    void StringArray::doomGet
    (Attribute &attr, const DoomReader &reader, int index) const {
        int size = reader.getInt(index);
        std::vector<string> value(size);

        for(int i = 0; i < size; i++) {
            value[i] = reader.getString(index + i + 1);
        }

        Attributes::setStringArray(attr, value);
    }


    void StringArray::doomPut
    (const Attribute &attr, DoomWriter &writer, int index) const {
        std::vector<string> array = Attributes::getStringArray(attr);
        int size = array.size();
        writer.putInt(index, size);

        for(int i = 0; i < size; i++) {
            writer.putString(index + i + 1, array[i]);
        }
    }


    const string &StringArray::getType() const {
        static string type = "string array";
        return type;
    }


    void StringArray::parse(Attribute &attr, Statement &stat, bool) const {
        Attributes::setStringArray(attr, Expressions::parseStringArray(stat));
    }


    void StringArray::parseComponent
    (Attribute &attr, Statement &statement, bool, int index) const {
        std::vector<string> array = Attributes::getStringArray(attr);

        if(AttributeBase *base = &attr.getBase()) {
            array = dynamic_cast<AValue<string>*>(base)->evaluate();
        }

        while(int(array.size()) < index) {
            array.push_back(string());
        }

        array[index-1] =
            Expressions::parseString(statement, "String value expected.");
        Attributes::setStringArray(attr, array);
    }

};
