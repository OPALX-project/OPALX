//
// Class ElementValues
//   The attribute values of one element, read from its description.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#include "Elements/ElementValues.h"

#include <utility>
#include <variant>

#include "AbstractObjects/Object.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/LogicalError.h"
#include "Utilities/OpalException.h"

ElementValues::ElementValues(opalx::run::ObjectDescription object) : object_m(std::move(object)) {}

const opalx::run::AttributeEntry* ElementValues::entry(const std::string& attribute) const {
    for (const opalx::run::AttributeEntry& entry : object_m.attributes) {
        if (entry.name == attribute) return &entry;
    }
    // Not described: the attribute has no value. Make sure the keyword has it at all.
    const Object* keyword = OpalData::getInstance()->findKeyword(object_m.keyword);
    if (keyword == nullptr || keyword->findAttribute(attribute) == nullptr) {
        throw LogicalError(
                "ElementValues::entry",
                "The " + object_m.keyword + " keyword has no attribute \"" + attribute + "\".");
    }
    return nullptr;
}

bool ElementValues::hasValue(const std::string& attribute) const {
    const opalx::run::AttributeEntry* found = entry(attribute);
    return found != nullptr && found->hasValue;
}

bool ElementValues::defaultUsed(const std::string& attribute) const {
    const opalx::run::AttributeEntry* found = entry(attribute);
    return found == nullptr || found->defaultUsed;
}

template <class T>
T ElementValues::get(const std::string& attribute, const char* type, T none) const {
    const opalx::run::AttributeEntry* found = entry(attribute);
    if (found == nullptr || !found->hasValue) return none;
    if (!found->error.empty()) throw OpalException("ElementValues::get", found->error);
    if (const T* value = std::get_if<T>(&found->value)) return *value;
    throw OpalException(
            "ElementValues::get", "Attribute \"" + attribute + "\" is not " + type + ".");
}

double ElementValues::getReal(const std::string& attribute) const {
    return get<double>(attribute, "real", 0.0);
}

bool ElementValues::getBool(const std::string& attribute) const {
    return get<bool>(attribute, "logical", false);
}

std::string ElementValues::getString(const std::string& attribute) const {
    return get<std::string>(attribute, "string", std::string());
}

std::vector<double> ElementValues::getRealArray(const std::string& attribute) const {
    return get<std::vector<double>>(attribute, "a real array", {});
}

std::vector<std::string> ElementValues::getStringArray(const std::string& attribute) const {
    return get<std::vector<std::string>>(attribute, "a string array", {});
}
