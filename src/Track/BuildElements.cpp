//
// Build the runtime elements of a RUN from its run description.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#include "Track/BuildElements.h"

#include <memory>
#include <string>
#include <variant>

#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Attributes/Bool.h"
#include "Attributes/BoolArray.h"
#include "Attributes/OpalString.h"
#include "Attributes/PredefinedString.h"
#include "Attributes/Real.h"
#include "Attributes/RealArray.h"
#include "Attributes/StringArray.h"
#include "Attributes/UpperCaseString.h"
#include "Attributes/UpperCaseStringArray.h"
#include "Utilities/OpalException.h"

namespace opalx::run {

    namespace {

        [[noreturn]] void cannotCopy(
                const ElementDescription& element, const AttributeEntry& entry,
                const std::string& why) {
            throw OpalException(
                    "opalx::run::buildElement",
                    "Element \"" + element.object.name + "\": attribute " + entry.name + " " + why);
        }

        /// Set @p attribute to the value the input gave it, with the setter of its handler.
        void copyValue(
                Attribute& attribute, const AttributeEntry& entry,
                const ElementDescription& element) {
            AttributeHandler* handler   = &attribute.getHandler();
            const AttributeValue& value = entry.value;
            if (dynamic_cast<Attributes::Real*>(handler) && std::holds_alternative<double>(value)) {
                Attributes::setReal(attribute, std::get<double>(value));
            } else if (
                    dynamic_cast<Attributes::Bool*>(handler)
                    && std::holds_alternative<bool>(value)) {
                Attributes::setBool(attribute, std::get<bool>(value));
            } else if (
                    dynamic_cast<Attributes::PredefinedString*>(handler)
                    && std::holds_alternative<std::string>(value)) {
                Attributes::setPredefinedString(attribute, std::get<std::string>(value));
            } else if (
                    dynamic_cast<Attributes::UpperCaseString*>(handler)
                    && std::holds_alternative<std::string>(value)) {
                Attributes::setUpperCaseString(attribute, std::get<std::string>(value));
            } else if (
                    dynamic_cast<Attributes::String*>(handler)
                    && std::holds_alternative<std::string>(value)) {
                Attributes::setString(attribute, std::get<std::string>(value));
            } else if (
                    dynamic_cast<Attributes::RealArray*>(handler)
                    && std::holds_alternative<std::vector<double>>(value)) {
                Attributes::setRealArray(attribute, std::get<std::vector<double>>(value));
            } else if (
                    dynamic_cast<Attributes::BoolArray*>(handler)
                    && std::holds_alternative<std::vector<bool>>(value)) {
                Attributes::setBoolArray(attribute, std::get<std::vector<bool>>(value));
            } else if (
                    dynamic_cast<Attributes::UpperCaseStringArray*>(handler)
                    && std::holds_alternative<std::vector<std::string>>(value)) {
                Attributes::setUpperCaseStringArray(
                        attribute, std::get<std::vector<std::string>>(value));
            } else if (
                    dynamic_cast<Attributes::StringArray*>(handler)
                    && std::holds_alternative<std::vector<std::string>>(value)) {
                Attributes::setStringArray(attribute, std::get<std::vector<std::string>>(value));
            } else {
                cannotCopy(element, entry, "has type " + entry.type + ", which cannot be copied.");
            }
        }

    }  // namespace

    std::shared_ptr<ElementBase> buildElement(const ElementDescription& element) {
        const ObjectDescription& input = element.object;
        std::unique_ptr<Object> object(
                OpalData::getInstance()->makeObject(input.keyword, input.name));
        auto* made = dynamic_cast<Element*>(object.get());
        if (made == nullptr) {
            throw OpalException(
                    "opalx::run::buildElement", "Element \"" + input.name
                                                        + "\" has no element keyword (\""
                                                        + input.keyword + "\").");
        }
        // Only what the input set is copied; every other attribute keeps the keyword's default,
        // so update() sees the same set and unset attributes as for the parsed object.
        for (const AttributeEntry& entry : input.attributes) {
            if (entry.defaultUsed || !entry.hasValue) continue;
            if (!entry.error.empty()) {
                cannotCopy(element, entry, "could not be evaluated: " + entry.error);
            }
            Attribute* attribute = made->findAttribute(entry.name);
            if (attribute == nullptr) cannotCopy(element, entry, "is not an attribute.");
            copyValue(*attribute, entry, element);
        }
        made->update();
        return made->getElementPtr();
    }

    std::vector<std::shared_ptr<ElementBase>> buildElements(const LineDescription& line) {
        std::vector<std::shared_ptr<ElementBase>> elements;
        elements.reserve(line.elements.size());
        for (const ElementDescription& element : line.elements) {
            std::shared_ptr<ElementBase> built = buildElement(element);
            const std::string& name            = built->getName();
            if (line.topology == "RING" && !(name.size() > 0 && name.front() == '#')) {
                built->setBeamlineMembership(BeamlineTopology::RING, line.name);
            }
            elements.push_back(std::move(built));
        }
        return elements;
    }

    double ringLength(const std::vector<std::shared_ptr<ElementBase>>& elements) {
        double length = 0.0;
        for (const auto& element : elements) {
            if (element->getType() == ElementType::MARKER) continue;
            length += element->getGeometry().getArcLength();
        }
        return length;
    }

}  // namespace opalx::run
