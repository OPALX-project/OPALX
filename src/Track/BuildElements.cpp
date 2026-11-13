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

#include <map>

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
#include "Elements/ElementValues.h"
#include "Elements/OpalDrift.h"
#include "Elements/OpalMarker.h"
#include "Elements/OpalMonitor.h"
#include "Elements/OpalProbe.h"
#include "Elements/OpalSolenoid.h"
#include "Utilities/LogicalError.h"
#include "Utilities/OpalException.h"

namespace opalx::run {

    namespace {

        /// What an element keyword builds: the runtime element type and the build function.
        struct ElementKind {
            ElementType type;
            BuildFunction build;  ///< nullptr: update() of a keyword object fills the element.
        };

        /// Every element keyword. One keyword, one runtime element type.
        const std::map<std::string, ElementKind>& elementKinds() {
            static const std::map<std::string, ElementKind> kinds{
                    {"BOX", {ElementType::BOX, nullptr}},
                    {"COLLIMATOR", {ElementType::COLLIMATOR, nullptr}},
                    {"CONSTANTEFIELDCAVITY", {ElementType::CONSTANTEFIELDCAVITY, nullptr}},
                    {"CONSTANTFOCUSING", {ElementType::CONSTANTFOCUSING, nullptr}},
                    {"CYCLOTRONSECTOR", {ElementType::CYCLOTRONSECTOR, nullptr}},
                    {"DRIFT", {ElementType::DRIFT, &OpalDrift::build}},
                    {"FIELDMAP", {ElementType::FIELDMAP, nullptr}},
                    {"LASER", {ElementType::LASER, nullptr}},
                    {"MARKER", {ElementType::MARKER, &OpalMarker::build}},
                    {"MONITOR", {ElementType::MONITOR, &OpalMonitor::build}},
                    {"MULTIPOLE", {ElementType::MULTIPOLE, nullptr}},
                    {"MULTIPOLET", {ElementType::MULTIPOLET, nullptr}},
                    {"PROBE", {ElementType::PROBE, &OpalProbe::build}},
                    {"QUADRUPOLE", {ElementType::MULTIPOLE, nullptr}},
                    {"RBEND", {ElementType::RBEND, nullptr}},
                    {"RFCAVITY", {ElementType::RFCAVITY, nullptr}},
                    {"SBEND", {ElementType::SBEND, nullptr}},
                    {"SCALINGFFAMAGNET", {ElementType::SCALINGFFAMAGNET, nullptr}},
                    {"SOLENOID", {ElementType::SOLENOID, &OpalSolenoid::build}},
                    {"TRAVELINGWAVE", {ElementType::TRAVELINGWAVE, nullptr}},
                    {"VARIABLE_RF_CAVITY", {ElementType::ANY, nullptr}},
                    {"VERTICALFFAMAGNET", {ElementType::ANY, nullptr}},
            };
            return kinds;
        }

        const ElementKind& kindOf(const std::string& keyword, const std::string& name) {
            const auto found = elementKinds().find(keyword);
            if (found == elementKinds().end()) {
                throw OpalException(
                        "opalx::run::buildElement",
                        "Element \"" + name + "\" has no element keyword (\"" + keyword + "\").");
            }
            return found->second;
        }

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

    std::string elementTypeOf(const std::string& keyword) {
        return ElementBase::getTypeString(kindOf(keyword, keyword).type);
    }

    std::shared_ptr<ElementBase> buildElement(const ElementDescription& element) {
        const ObjectDescription& input = element.object;
        const ElementKind& kind        = kindOf(input.keyword, input.name);
        for (const AttributeEntry& entry : input.attributes) {
            if (entry.defaultUsed || !entry.hasValue || entry.error.empty()) continue;
            cannotCopy(element, entry, "could not be evaluated: " + entry.error);
        }

        std::shared_ptr<ElementBase> built;
        if (kind.build != nullptr) {
            built = kind.build(ElementValues(input));
        } else {
            std::unique_ptr<Object> object(
                    OpalData::getInstance()->makeObject(input.keyword, input.name));
            auto* made = dynamic_cast<Element*>(object.get());
            // Only what the input set is copied; every other attribute keeps the keyword's
            // default, so update() sees the same set and unset attributes as for the parsed
            // object.
            for (const AttributeEntry& entry : input.attributes) {
                if (entry.defaultUsed || !entry.hasValue) continue;
                Attribute* attribute = made->findAttribute(entry.name);
                if (attribute == nullptr) cannotCopy(element, entry, "is not an attribute.");
                copyValue(*attribute, entry, element);
            }
            made->update();
            built = made->getElementPtr();
        }
        if (built->getType() != kind.type) {
            throw LogicalError(
                    "opalx::run::buildElement",
                    "The " + input.keyword + " keyword built a "
                            + ElementBase::getTypeString(built->getType()) + " element.");
        }
        return built;
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
