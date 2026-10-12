/**
 * @file ElementInteractionManager.cpp
 * @brief Implementation of the generic element-interaction dispatcher.
 */

#include "Algorithms/ElementInteractionManager.h"

#include "AbsBeamline/ElementBase.h"
#include "Utilities/OpalException.h"

void ElementInteractionManager::initialize(const std::set<std::shared_ptr<ElementBase>>& elements) {
    clear();
    std::vector<Entry> prepared;
    for (const auto& element : elements) {
        if (!element) {
            continue;
        }
        auto interaction = element->createInteraction();
        if (interaction) {
            if (!prepared.empty()) {
                throw OpalException(
                        "ElementInteractionManager::initialize",
                        "Only one stateful collective interaction (BEAMBEAM) is supported per "
                        "tracking run; found '" + prepared.front().element->getName() + "' and '"
                                + element->getName() + "'. Use a single BeamBeam element in a "
                                                       "single-pass LINE.");
            }
            prepared.push_back(Entry{element, std::move(interaction)});
        }
    }
    entries_m = std::move(prepared);
}

void ElementInteractionManager::clear() { entries_m.clear(); }

ElementInteractionResult ElementInteractionManager::execute(
        ElementInteractionPhase phase, ElementInteractionContext& context) {
    ElementInteractionResult accumulated;
    for (auto& entry : entries_m) {
        const ElementInteractionResult result = entry.interaction->execute(phase, context);
        accumulated.selfFieldHandled = accumulated.selfFieldHandled || result.selfFieldHandled;
        if (phase == ElementInteractionPhase::SelfField && result.selfFieldHandled) {
            break;
        }
    }
    return accumulated;
}

bool ElementInteractionManager::freezesFieldMesh() const noexcept {
    for (const auto& entry : entries_m) {
        if (entry.interaction->freezesFieldMesh()) {
            return true;
        }
    }
    return false;
}

bool ElementInteractionManager::suppressesDefaultSelfField() const noexcept {
    for (const auto& entry : entries_m) {
        if (entry.interaction->suppressesDefaultSelfField()) {
            return true;
        }
    }
    return false;
}
