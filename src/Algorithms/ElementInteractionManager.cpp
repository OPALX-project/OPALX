/**
 * @file ElementInteractionManager.cpp
 * @brief Implementation of the generic element-interaction dispatcher.
 */

#include "Algorithms/ElementInteractionManager.h"

#include "AbsBeamline/ElementBase.h"

void ElementInteractionManager::initialize(const std::set<std::shared_ptr<ElementBase>>& elements) {
    clear();
    entries_m.reserve(elements.size());
    for (const auto& element : elements) {
        if (!element) {
            continue;
        }
        auto interaction = element->createInteraction();
        if (interaction) {
            entries_m.push_back(Entry{element, std::move(interaction)});
        }
    }
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
