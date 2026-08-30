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
 * element interaction and diagnostics code. Runtime ownership belongs to the
 * per-run BeamBeamInteraction; only the shared type definitions live here.
 */
namespace BEAMBEAM {

    /** Fixed longitudinal extent of the BeamBeam collective-field solve [m]. */
    inline constexpr double fieldWindowLength = 20.0e-3;

    /**
     * Maximum longitudinal-to-transverse cell ratio in the primary rest frame.
     *
     * The integrated Green function was finite through approximately 1.95e5 and
     * non-finite at 3.84e5 for the deterministic BeamBeam manufactured source.
     * The lower value leaves headroom below the measured transition.
     */
    inline constexpr double maximumRestFrameCellAspectRatio = 1.5e5;

    /**
     * @brief Minimum transverse span required by the BeamBeam cell-aspect constraint.
     *
     * Lorentz transformation to the primary rest frame stretches the longitudinal
     * cell by @p gamma. For a fixed longitudinal window, this returns the transverse
     * span whose cell width limits @f$\gamma\,\Delta z/\Delta x@f$ to @p maxAspect.
     */
    inline double minimumTransverseSpan(
            double longitudinalLength, std::size_t longitudinalCells, std::size_t transverseCells,
            double gamma, double maxAspect = maximumRestFrameCellAspectRatio) {
        if (longitudinalCells < 2 || transverseCells < 2 || longitudinalLength <= 0.0
            || gamma <= 0.0 || maxAspect <= 0.0) {
            return 0.0;
        }
        const double restFrameLongitudinalCell =
                gamma * longitudinalLength / static_cast<double>(longitudinalCells - 1);
        return restFrameLongitudinalCell * static_cast<double>(transverseCells - 1) / maxAspect;
    }

    /**
     * @brief Lifecycle of a single BeamBeam passage during tracking.
     */
    enum class WindowState { Inactive, Active, Completed };

    /**
     * @brief Static BeamBeam configuration derived from the element attributes.
     */
    struct Config {
        bool visualize = false;
        /**
         * Suppress the BeamBeam collective kick on physical source container 0.
         *
         * The source still deposits charge and the resulting mesh field remains
         * available to the copied source model and passive witness containers.
         * This switch does not suppress ordinary tracking or external fields.
         */
        bool rigidSource = false;
        std::optional<double> copyTime;
        std::vector<std::size_t> witnessContainers;
    };

    /**
     * @brief Whether the physical source samples its BeamBeam collective field.
     *
     * Keeping this policy explicit ensures that rigid-source mode changes only
     * the source response; field deposition and witness gathering are unchanged.
     */
    inline bool sourceCollectiveKickEnabled(const Config& config) { return !config.rigidSource; }

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
     * These values are not intrinsic element-local coordinates. For ELEMEDGE
     * placement, beginS is the placed entrance and endS is beginS plus the
     * declared element length. The OrbitThreader range is intentionally not used
     * as the physical end because a short tracking horizon may clip that map.
     * For 6D-pose placement, the threaded entrance supplies the path-length
     * anchor because no ELEMEDGE coordinate exists.
     */
    struct ActualGeometry {
        double interactionPointS = 0.0;
        double beginS            = 0.0;
        double endS              = 0.0;
        double length            = 0.0;
        Config config;
    };

    /**
     * @brief Derive the BeamBeam interaction point from placed element geometry.
     *
     * The interaction point has no independent input coordinate: it is always
     * the midpoint between the placed element entrance and exit,
     * @f[
     *   s_\mathrm{IP} = s_\mathrm{begin}
     *                   + \frac{s_\mathrm{end}-s_\mathrm{begin}}{2}.
     * @f]
     * All arguments and the return value are lab-frame path lengths in metres.
     */
    inline double interactionPointAtElementMidpoint(double beginS, double endS) {
        return beginS + 0.5 * (endS - beginS);
    }

    /** @brief Upstream boundary of the fixed collective-field window [m]. */
    inline double fieldWindowBegin(double interactionPointS) {
        return interactionPointS - 0.5 * fieldWindowLength;
    }

    /** @brief Downstream boundary of the fixed collective-field window [m]. */
    inline double fieldWindowEnd(double interactionPointS) {
        return interactionPointS + 0.5 * fieldWindowLength;
    }

    /** @brief True once the complete physical primary has crossed the downstream field boundary. */
    inline bool sourceFullyExitedFieldWindow(double sourceTailS, double interactionPointS) {
        return sourceTailS > fieldWindowEnd(interactionPointS);
    }

    /**
     * @brief Test whether the tracked source bunch overlaps its copied counter-propagating source.
     *
     * In the current copy model the second high-energy source bunch is represented by mirroring the
     * tracked source about the interaction point @f$s_\mathrm{IP}@f$. A tracked absolute
     * longitudinal interval
     * @f[
     *   [s_\mathrm{tail}, s_\mathrm{head}]
     * @f]
     * therefore has the copied-source interval
     * @f[
     *   [2s_\mathrm{IP}-s_\mathrm{head}, 2s_\mathrm{IP}-s_\mathrm{tail}].
     * @f]
     * The two finite bunches overlap if the two closed intervals intersect,
     * @f[
     *   s_\mathrm{tail} \le 2s_\mathrm{IP}-s_\mathrm{tail}
     *   \quad\land\quad
     *   2s_\mathrm{IP}-s_\mathrm{head} \le s_\mathrm{head}.
     * @f]
     *
     * @param sourceTailS Trailing-edge path length of the tracked source bunch.
     * @param sourceHeadS Leading-edge path length of the tracked source bunch.
     * @param geometry Actual BeamBeam window geometry in path-length coordinates.
     * @return True while the tracked source and its copied source overlap.
     */
    inline bool copiedSourceBunchesOverlap(
            double sourceTailS, double sourceHeadS, const ActualGeometry& geometry) {
        const double copiedTailS = 2.0 * geometry.interactionPointS - sourceHeadS;
        const double copiedHeadS = 2.0 * geometry.interactionPointS - sourceTailS;
        return sourceTailS <= copiedHeadS && copiedTailS <= sourceHeadS;
    }

    /**
     * @brief Test whether the mirrored-source copy model should be active.
     *
     * A BeamBeam element may specify @c COPY_TIME in seconds. A missing value disables copied
     * fields; otherwise the copied-source contribution is included once
     * @f[
     *   t \ge t_\mathrm{copy}.
     * @f]
     *
     * @param currentTime Current simulation time in seconds.
     * @param copyTime Optional copied-source activation time in seconds.
     * @return True if copied fields should be included in the active BeamBeam solve.
     */
    inline bool copyTimeReached(double currentTime, const std::optional<double>& copyTime) {
        return copyTime.has_value() && currentTime >= *copyTime;
    }

    /**
     * @brief Interaction-owned runtime state for the BeamBeam model.
     *
     * BeamBeamInteraction keeps the lifecycle, currently active actual geometry,
     * and saved field-domain rollback state together in this structure.
     */
    template <class SavedFieldDomainState>
    struct Runtime {
        WindowState state = WindowState::Inactive;
        std::optional<ActualGeometry> geometry;
        std::optional<SavedFieldDomainState> savedFieldDomain;
        bool sourceBunchesOverlap = false;
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
