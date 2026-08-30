/**
 * @file ElementInteraction.h
 * @brief Generic runtime contract for element-scoped interactions.
 */

#ifndef OPALX_ElementInteraction_HH
#define OPALX_ElementInteraction_HH

#include "PartBunch/PartBunchFwd.h"

class CoordinateSystemTrafo;
class Inform;
class OrbitThreader;

/**
 * @brief Stable integration points at which a runtime element interaction may act.
 *
 * The tracker dispatches these phases without knowing which concrete element
 * interactions are present. An interaction ignores phases that it does not use.
 */
enum class ElementInteractionPhase {
    SelfField,
    AfterEmission,
    Diagnostics,
};

/**
 * @brief Non-owning view of the tracking state and services for one phase dispatch.
 *
 * Currently only `BeamBeamInteraction` uses `ElementInteractionContext`.
 * Every other element inherits the null `ElementBase::createInteraction()`
 * implementation and continues through the normal visitor/`apply()` mechanism.
 * The context provides an architectural boundary through which
 * `ParallelTracker` supplies runtime services to optional, stateful collective
 * interactions without knowing about `BeamBeamInteraction` specifically.
 *
 * `ParallelTracker` constructs a context immediately before calling
 * `ElementInteractionManager::execute()`. The manager passes that context to
 * the registered interactions in sequence. The context owns none of the
 * referenced objects: an interaction may use them only while its `execute()`
 * call is active and must not retain their addresses or references afterward.
 * `bunch` is always available; the pointer members are optional services whose
 * validity depends on the dispatched phase.
 *
 * During `SelfField`, the primary-particle positions have already been rotated
 * from the reference-path frame into the instantaneous beam frame. The two
 * coordinate transformations describe this rotation and its inverse; the
 * inverse is used to return positions and computed fields to the reference
 * frame after the solve. The orbit threader supports lookup of placed elements
 * at the current longitudinal location, and `endOfLine` lets an interaction
 * report an out-of-range lookup to the tracker.
 *
 * During `AfterEmission` and `Diagnostics`, only `bunch` and, when supplied,
 * `message` are guaranteed. Interactions must check every optional pointer
 * before use and should reject a phase invocation that lacks a service required
 * for correct execution.
 */
struct ElementInteractionContext {
    PartBunch_t& bunch;  ///< Active multi-container bunch; valid in every phase.

    /// Source reference-orbit lookup; supplied during `SelfField`.
    OrbitThreader* sourceOrbitThreader = nullptr;

    /// Reference-path to beam-frame rotation already applied before `SelfField`.
    const CoordinateSystemTrafo* referenceToBeamCSTrafo = nullptr;

    /// Inverse beam-frame to reference-path transformation for solve results.
    const CoordinateSystemTrafo* beamToReferenceCSTrafo = nullptr;

    /// Optional tracker log stream; interactions may create a local fallback.
    Inform* message = nullptr;

    /// Optional mutable tracker flag for an orbit-threader out-of-range result.
    bool* endOfLine = nullptr;
};

/** @brief Result flags accumulated across one generic interaction dispatch. */
struct ElementInteractionResult {
    bool selfFieldHandled = false;
};

/**
 * @brief Per-run behavior associated with a placed lattice element.
 *
 * Element definitions contain immutable configuration. Implementations of this
 * interface contain mutable state whose lifetime is one tracking run.
 */
class ElementInteraction {
public:
    virtual ~ElementInteraction() = default;

    /** @brief Execute one generic tracker phase. */
    virtual ElementInteractionResult execute(
            ElementInteractionPhase phase, ElementInteractionContext& context) = 0;

    /** @brief Whether the interaction currently requires its field mesh to remain fixed. */
    virtual bool freezesFieldMesh() const noexcept { return false; }

    /** @brief Whether the interaction suppresses the ordinary source self-field solve. */
    virtual bool suppressesDefaultSelfField() const noexcept { return false; }
};

#endif  // OPALX_ElementInteraction_HH
