//
// Build the runtime elements of a RUN from its run description.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#include "Track/BuildElements.h"

#include <map>
#include <memory>
#include <string>

#include "AbsBeamline/ElementBase.h"
#include "Elements/ElementValues.h"
#include "Elements/OpalAsymmetricEnge.h"
#include "Elements/OpalBox.h"
#include "Elements/OpalCavity.h"
#include "Elements/OpalCollimator.h"
#include "Elements/OpalConstantEFieldCavity.h"
#include "Elements/OpalConstantFocusing.h"
#include "Elements/OpalCyclotronSector.h"
#include "Elements/OpalDrift.h"
#include "Elements/OpalEnge.h"
#include "Elements/OpalFieldmapElement.h"
#include "Elements/OpalLaser.h"
#include "Elements/OpalMarker.h"
#include "Elements/OpalMonitor.h"
#include "Elements/OpalMultipole.h"
#include "Elements/OpalMultipoleT.h"
#include "Elements/OpalPolynomialTimeDependence.h"
#include "Elements/OpalProbe.h"
#include "Elements/OpalQuadrupole.h"
#include "Elements/OpalRBend.h"
#include "Elements/OpalSBend.h"
#include "Elements/OpalScalingFFAMagnet.h"
#include "Elements/OpalSinusoidalTimeDependence.h"
#include "Elements/OpalSolenoid.h"
#include "Elements/OpalSplineTimeDependence.h"
#include "Elements/OpalTanh.h"
#include "Elements/OpalTravelingWave.h"
#include "Elements/OpalVariableRFCavity.h"
#include "Elements/OpalVerticalFFAMagnet.h"
#include "Utilities/LogicalError.h"
#include "Utilities/OpalException.h"

namespace opalx::run {

    namespace {

        /// What an element keyword builds: the runtime element type and the build function.
        struct ElementKind {
            ElementType type;
            BuildFunction build;
        };

        /// Every element keyword. One keyword, one runtime element type.
        const std::map<std::string, ElementKind>& elementKinds() {
            static const std::map<std::string, ElementKind> kinds{
                    {"BOX", {ElementType::BOX, &OpalBox::build}},
                    {"COLLIMATOR", {ElementType::COLLIMATOR, &OpalCollimator::build}},
                    {"CONSTANTEFIELDCAVITY",
                     {ElementType::CONSTANTEFIELDCAVITY, &OpalConstantEFieldCavity::build}},
                    {"CONSTANTFOCUSING",
                     {ElementType::CONSTANTFOCUSING, &OpalConstantFocusing::build}},
                    {"CYCLOTRONSECTOR",
                     {ElementType::CYCLOTRONSECTOR, &OpalCyclotronSector::build}},
                    {"DRIFT", {ElementType::DRIFT, &OpalDrift::build}},
                    {"FIELDMAP", {ElementType::FIELDMAP, &OpalFieldmapElement::build}},
                    {"LASER", {ElementType::LASER, &OpalLaser::build}},
                    {"MARKER", {ElementType::MARKER, &OpalMarker::build}},
                    {"MONITOR", {ElementType::MONITOR, &OpalMonitor::build}},
                    {"MULTIPOLE", {ElementType::MULTIPOLE, &OpalMultipole::build}},
                    {"MULTIPOLET", {ElementType::MULTIPOLET, &OpalMultipoleT::build}},
                    {"PROBE", {ElementType::PROBE, &OpalProbe::build}},
                    {"QUADRUPOLE", {ElementType::MULTIPOLE, &OpalQuadrupole::build}},
                    {"RBEND", {ElementType::RBEND, &OpalRBend::build}},
                    {"RFCAVITY", {ElementType::RFCAVITY, &OpalCavity::build}},
                    {"SBEND", {ElementType::SBEND, &OpalSBend::build}},
                    {"SCALINGFFAMAGNET",
                     {ElementType::SCALINGFFAMAGNET, &OpalScalingFFAMagnet::build}},
                    {"SOLENOID", {ElementType::SOLENOID, &OpalSolenoid::build}},
                    {"TRAVELINGWAVE", {ElementType::TRAVELINGWAVE, &OpalTravelingWave::build}},
                    {"VARIABLE_RF_CAVITY", {ElementType::ANY, &OpalVariableRFCavity::build}},
                    {"VERTICALFFAMAGNET", {ElementType::ANY, &OpalVerticalFFAMagnet::build}},
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

    }  // namespace

    std::string elementTypeOf(const std::string& keyword) {
        return ElementBase::getTypeString(kindOf(keyword, keyword).type);
    }

    std::shared_ptr<ElementBase> buildElement(const ElementDescription& element) {
        const ObjectDescription& input = element.object;
        const ElementKind& kind        = kindOf(input.keyword, input.name);
        for (const AttributeEntry& entry : input.attributes) {
            if (entry.defaultUsed || !entry.hasValue || entry.error.empty()) continue;
            throw OpalException(
                    "opalx::run::buildElement",
                    "Element \"" + input.name + "\": attribute " + entry.name
                            + " could not be evaluated: " + entry.error);
        }

        std::shared_ptr<ElementBase> built = kind.build(ElementValues(input));
        if (built->getType() != kind.type) {
            throw LogicalError(
                    "opalx::run::buildElement",
                    "The " + input.keyword + " keyword built a "
                            + ElementBase::getTypeString(built->getType()) + " element.");
        }
        return built;
    }

    namespace {

        using RegisterModel = void (*)(const ElementValues&);

        /// The keywords of time dependences and end-field models, with their registerModel().
        const std::map<std::string, RegisterModel>& modelKeywords() {
            static const std::map<std::string, RegisterModel> table = {
                    {"POLYNOMIAL_TIME_DEPENDENCE", &OpalPolynomialTimeDependence::registerModel},
                    {"SINUSOIDAL_TIME_DEPENDENCE", &OpalSinusoidalTimeDependence::registerModel},
                    {"SPLINE_TIME_DEPENDENCE", &OpalSplineTimeDependence::registerModel},
                    {"ENGE", &OpalEnge::registerModel},
                    {"TANH", &OpalTanh::registerModel},
                    {"ASYMMETRIC_ENGE", &OpalAsymmetricEnge::registerModel},
            };
            return table;
        }

    }  // namespace

    bool isModelKeyword(const std::string& keyword) { return modelKeywords().count(keyword) != 0; }

    void registerModels(const std::vector<ObjectDescription>& models) {
        for (const ObjectDescription& model : models) {
            const auto found = modelKeywords().find(model.keyword);
            if (found == modelKeywords().end()) {
                throw LogicalError(
                        "opalx::run::registerModels",
                        "\"" + model.keyword + "\" is not a time dependence or end-field model.");
            }
            found->second(ElementValues(model));
        }
    }

    std::vector<std::shared_ptr<ElementBase>> buildElements(const LineDescription& line) {
        registerModels(line.models);
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
