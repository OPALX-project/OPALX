//
// Class ElementValues
//   The attribute values of one element, read from its description.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#ifndef OPALX_ELEMENT_VALUES_H
#define OPALX_ELEMENT_VALUES_H

#include <string>
#include <vector>

#include "Track/RunDescription.h"

/**
 * @brief The attribute values of one element, read from its ObjectDescription.
 *
 * The functions follow the rules of the Attributes::get* functions on the parsed object: an
 * attribute without a value reads as 0, false, "" or an empty list. An attribute whose
 * evaluation failed throws its error when it is read.
 */
class ElementValues {
public:
    explicit ElementValues(opalx::run::ObjectDescription object);

    /// The name the input gave the element.
    const std::string& getName() const { return object_m.name; }

    /// Same as bool(attribute): the attribute has a value, from the input or a default.
    bool hasValue(const std::string& attribute) const;

    /// Same as Attribute::defaultUsed(): the input did not set the attribute.
    bool defaultUsed(const std::string& attribute) const;

    double getReal(const std::string& attribute) const;
    bool getBool(const std::string& attribute) const;
    std::string getString(const std::string& attribute) const;
    std::vector<double> getRealArray(const std::string& attribute) const;
    std::vector<std::string> getStringArray(const std::string& attribute) const;

private:
    /**
     * @brief The entry of @p attribute, or nullptr if it has no value (describeObject() leaves
     *        those out).
     * @throws LogicalError if the keyword has no attribute of that name.
     */
    const opalx::run::AttributeEntry* entry(const std::string& attribute) const;

    /// The value of @p attribute as a @p T, or @p none if it has no value.
    template <class T>
    T get(const std::string& attribute, const char* type, T none) const;

    opalx::run::ObjectDescription object_m;
};

#endif  // OPALX_ELEMENT_VALUES_H
