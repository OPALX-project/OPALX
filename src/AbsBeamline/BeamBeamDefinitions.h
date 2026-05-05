#ifndef CLASSIC_BeamBeamDefinitions_HH
#define CLASSIC_BeamBeamDefinitions_HH

#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

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
        std::optional<double> xAperture;
        std::optional<double> yAperture;
        std::vector<std::size_t> witnessContainers;
    };

    /**
     * @brief Decode the numeric BeamBeam witness-container bit mask.
     *
     * The parser stores `WITNESS_CONTAINERS="1,2"` as bits 1 and 2 in a numeric
     * element attribute because CLASSIC element user attributes are scalar doubles.
     * A zero mask represents `WITNESS_CONTAINERS="NONE"` and preserves the legacy
     * behavior: only the source container samples the BeamBeam self-field.
     *
     * @param maskValue Scalar element attribute containing the bit mask.
     * @returns Container indices that passively gather the source BeamBeam field.
     */
    inline std::vector<std::size_t> decodeWitnessContainerMask(double maskValue) {
        std::vector<std::size_t> witnessContainers;
        if (maskValue <= 0.0) {
            return witnessContainers;
        }

        auto mask = static_cast<unsigned long long>(std::llround(maskValue));
        for (std::size_t index = 0; mask != 0; ++index) {
            if ((mask & 1ULL) != 0) {
                witnessContainers.push_back(index);
            }
            mask >>= 1U;
        }
        return witnessContainers;
    }

    /**
     * @brief Longitudinal shift from a target container frame into the source frame.
     *
     * Particle coordinates are stored relative to each container reference path. For a target
     * container at path length @f$s_t@f$ and the BeamBeam source at @f$s_s@f$, the source-frame
     * longitudinal coordinate is
     * @f[
     *   z_s = z_t + (s_t - s_s).
     * @f]
     *
     * @param sourceS Path length of the source BeamBeam container.
     * @param targetS Path length of the target/witness container.
     * @returns Additive longitudinal offset applied before the source-frame rotation.
     */
    inline double longitudinalOffsetToSourceFrame(double sourceS, double targetS) {
        return targetS - sourceS;
    }

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
