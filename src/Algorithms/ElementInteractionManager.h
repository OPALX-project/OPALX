/**
 * @file ElementInteractionManager.h
 * @brief Owns per-run interactions created by placed lattice elements.
 */

#ifndef OPALX_ElementInteractionManager_HH
#define OPALX_ElementInteractionManager_HH

#include "Algorithms/ElementInteraction.h"

#include <cstddef>
#include <memory>
#include <set>
#include <vector>

class ElementBase;

/**
 * @brief Generic owner and dispatcher for element-scoped runtime interactions.
 *
 * The manager retains each placed element alongside the interaction it created,
 * guaranteeing that the interaction's configuration reference remains valid.
 */
class ElementInteractionManager {
public:
    /** @brief Rebuild interactions for one prepared runtime beamline. */
    void initialize(const std::set<std::shared_ptr<ElementBase>>& elements);

    /** @brief Remove all per-run interactions and their state. */
    void clear();

    /** @brief Dispatch a tracker phase to all registered interactions. */
    ElementInteractionResult execute(
            ElementInteractionPhase phase, ElementInteractionContext& context);

    /** @brief True when any active interaction requires a fixed field mesh. */
    bool freezesFieldMesh() const noexcept;

    /** @brief True when any interaction suppresses the ordinary source self-field solve. */
    bool suppressesDefaultSelfField() const noexcept;

    /** @brief Number of placed elements that supplied a runtime interaction. */
    std::size_t size() const noexcept { return entries_m.size(); }

private:
    struct Entry {
        std::shared_ptr<ElementBase> element;
        std::unique_ptr<ElementInteraction> interaction;
    };

    std::vector<Entry> entries_m;
};

#endif  // OPALX_ElementInteractionManager_HH
