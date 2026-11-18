//
// Class ObjectValues
//   The attribute values of one input object, read from its description.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#include "Track/ObjectValues.h"

#include <utility>
#include <variant>

#include "Utilities/LogicalError.h"
#include "Utilities/OpalException.h"

ObjectValues::ObjectValues(opalx::run::ObjectDescription object) : object_m(std::move(object)) {}

const opalx::run::AttributeEntry& ObjectValues::entry(const std::string& attribute) const {
    for (const opalx::run::AttributeEntry& entry : object_m.attributes) {
        if (entry.name == attribute) return entry;
    }
    throw LogicalError(
            "ObjectValues::entry",
            "The " + object_m.keyword + " keyword has no attribute \"" + attribute + "\".");
}

bool ObjectValues::hasValue(const std::string& attribute) const {
    return entry(attribute).hasValue;
}

bool ObjectValues::defaultUsed(const std::string& attribute) const {
    return entry(attribute).defaultUsed;
}

template <class T>
T ObjectValues::get(const std::string& attribute, const char* type, T none) const {
    const opalx::run::AttributeEntry& found = entry(attribute);
    if (!found.hasValue) return none;
    if (!found.error.empty()) throw OpalException("ObjectValues::get", found.error);
    if (const T* value = std::get_if<T>(&found.value)) return *value;
    throw OpalException(
            "ObjectValues::get", "Attribute \"" + attribute + "\" is not " + type + ".");
}

double ObjectValues::getReal(const std::string& attribute) const {
    return get<double>(attribute, "real", 0.0);
}

bool ObjectValues::getBool(const std::string& attribute) const {
    return get<bool>(attribute, "logical", false);
}

std::string ObjectValues::getString(const std::string& attribute) const {
    return get<std::string>(attribute, "string", std::string());
}

std::vector<double> ObjectValues::getRealArray(const std::string& attribute) const {
    return get<std::vector<double>>(attribute, "a real array", {});
}

std::vector<std::string> ObjectValues::getStringArray(const std::string& attribute) const {
    return get<std::vector<std::string>>(attribute, "a string array", {});
}
