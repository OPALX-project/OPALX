//
// Class UString
//   This class is used to parse attributes of type string that should
//   the upper case, i.e. it returns upper case strings.
//
// Copyright (c) 2020, Christof Metzger-Kraus, Open Sourcerer
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Attributes/UString.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/SValue.h"
#include "Utilities/OpalException.h"
#include "Utilities/ParseError.h"
#include "Utilities/Util.h"

using namespace Expressions;

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