#ifndef CLASSIC_BeamBeamDefinitions_HH
#define CLASSIC_BeamBeamDefinitions_HH

#include <optional>

/**
 * @brief Beam-beam-specific shared type vocabulary.
 *
 * This namespace centralizes the stable BeamBeam concepts used across the
 * element, tracker, and diagnostics code. The actual runtime ownership still
 * belongs to the tracker; only the type definitions live here.
 */
namespace BEAMBEAM {

/**
 * @brief Lifecycle of a single BeamBeam passage during tracking.
 */
enum class WindowState { Inactive, Active, Completed };

/**
 * @brief Static BeamBeam configuration derived from the element attributes.
 */
struct Config {
    bool copyModel = false;
    bool visualize = false;
};

/**
 * @brief Actual lab-frame/path-length geometry of the currently active placed
 * BeamBeam element.
 *
 * These values are not intrinsic element-local coordinates. They are derived
 * from the placed element range in the lattice and therefore represent the
 * actual BeamBeam window seen by the tracker at runtime.
 */
struct ActualGeometry {
    double interactionPointS = 0.0;
    double beginS            = 0.0;
    double endS              = 0.0;
    double length            = 0.0;
    Config config;
};

/**
 * @brief Tracker-owned runtime state for the BeamBeam model.
 *
 * The tracker keeps the lifecycle, the currently active actual geometry, and
 * the saved field-domain rollback state together in this structure.
 */
template <class SavedFieldDomainState>
struct Runtime {
    WindowState state = WindowState::Inactive;
    std::optional<ActualGeometry> geometry;
    std::optional<SavedFieldDomainState> savedFieldDomain;
};

/**
 * @brief Diagnostics-only BeamBeam state.
 *
 * These flags control visualization or one-shot debug output and are separate
 * from the physics/runtime lifecycle.
 */
struct Diagnostics {
    bool frameObserved          = false;
    bool entryRhoSnapshotDumped = false;
};

}  // namespace BEAMBEAM

#endif  // CLASSIC_BeamBeamDefinitions_HH
