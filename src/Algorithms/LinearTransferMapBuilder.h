// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_LINEAR_TRANSFER_MAP_BUILDER_H
#define OPAL_LINEAR_TRANSFER_MAP_BUILDER_H

#include <optional>
#include "Algorithms/ExternalFieldRayTracker.h"
#include "Algorithms/IndexMap.h"
#include "Structure/LinearTransferMap.h"

class OpalBeamline;
class PartData;

/**
 * @brief Construct external-field linear maps from a sampled reference orbit.
 *
 * Owns Bishop-frame transport, boundary refinement, the perturbed rays,
 * the pivoted finite-difference solve and ordered composition. The caller owns
 * the reference pass and attaches the returned maps to their runtime elements.
 * The builder neither modifies the beamline nor writes diagnostics to stdout.
 *
 * Two independent selections are used: nominal body containment splits the sampled
 * orbit into constant-owner segments, while each private ray selects field supports
 * at its own integration stages. Owners are attachment destinations, not an isolated
 * field model. Empty-owner intervals are retained in the full ordered product.
 * BeamlineMembership is logical sequence metadata and does not participate in either
 * selection. The supplied beamline and particle data must outlive this builder.
 * See LinearTransferMap for coordinates, units and the finite-difference equations.
 */
class LinearTransferMapBuilder {
public:
    /**
     * @brief Explicit numerical configuration, independent of global OPTION state.
     *
     * Level zero uses twelve centered-difference rays per segment. With L Richardson
     * levels, use \f$12(L+1)\f$ rays at amplitudes \f$\epsilon_j/2^k\f$,
     * \f$k=0,\ldots,L\f$. The tableau is
     * \f[
     * \begin{aligned}
     * R_{k,0}&=M(\boldsymbol\epsilon/2^k),\\
     * R_{k,j}&=R_{k,j-1}\\
     * &\quad+\frac{R_{k,j-1}-R_{k-1,j-1}}{4^j-1}.
     * \end{aligned}
     * \f]
     * Here \f$1\le j\le k\le L\f$ in the extrapolation recurrence.
     * Return \f$R_{L,L}\f$, formally \f$O(\epsilon^{2(L+1)})\f$ for a sufficiently smooth ray map.
     * This extrapolates perturbation size, not DT: integration, field interpolation,
     * boundary and roundoff errors remain. Levels are limited to four to bound the
     * ray cost (60 per segment) and discourage unsupported high-order precision claims.
     * The caller must use integrationMethod for the supplied reference orbit as well.
     */
    struct Settings {
        static constexpr unsigned maximumRichardsonLevels = 4;
        /// Starting amplitudes in (x,x',y,y',zeta,delta) order; all default to 1e-3.
        /// Metres for x,y,zeta; dimensionless for x',y',delta. These are not DT values.
        std::array<double, 6> finiteDifferenceSteps{1.e-3, 1.e-3, 1.e-3, 1.e-3, 1.e-3, 1.e-3};
        /// Additional amplitude halvings, from 0 (default) through maximumRichardsonLevels.
        unsigned richardsonLevels{0};
        /// Fixed-nominal-step integrator for reference and shadow rays; default BORIS.
        ExternalFieldRayTracker::IntegrationMethod integrationMethod{
                ExternalFieldRayTracker::IntegrationMethod::BORIS};
        /// Require finite positive amplitudes that remain nonzero after refinement,
        /// delta amplitude < 1, an allowed level count and a supported integration method.
        /// @throws OpalException If any of these conditions fails.
        void validate() const;
    };

    /// Sampled reference state including its transported frame and compensated low parts.
    /// Does not cache an active set: support and nominal-body queries are independent.
    struct ReferenceSample {
        LinearTransferMapReference state;
    };
    /// One unique nominal-owner interval; a caller may attach a map copy to each owner.
    struct Segment {
        LinearTransferMap map;
        /// Runtime nominal body owners, possibly empty; NOT the field-contributor set.
        /// Shared pointers keep these occurrences alive while the result exists.
        IndexMap::value_t owners;
    };
    /// Value-owned calculation results, independent of any later attachment to elements.
    struct Result {
        /// Includes unowned intervals (possibly containing fringe fields); each overlap once.
        std::vector<Segment> segments;
        /// M_N ... M_1 in reference-path order; absent if no nondegenerate segments were built.
        std::optional<matrix6x6_t> combined;
        /// |det(combined)-1|; populated exactly when combined is populated.
        std::optional<double> determinantResidual;
        /// Canonical-form diagnostic defined in LinearTransferMap; same availability as combined.
        std::optional<double> symplecticResidual;
    };

    /// Construct with Settings defaults; does not read global OPTION values.
    LinearTransferMapBuilder(OpalBeamline& beamline, const PartData& reference, double dt);
    /// @param beamline Placed runtime elements used for body and external-field queries.
    /// @param reference Particle rest energy and charge; retained by reference.
    /// @param dt Finite nonzero nominal ray time step [s]; support events may subdivide it.
    /// @param settings Numerical configuration, copied and validated on construction.
    /// @throws OpalException If dt or settings are invalid.
    LinearTransferMapBuilder(
            OpalBeamline& beamline, const PartData& reference, double dt, const Settings& settings);
    /**
     * @brief Refine nominal boundaries, calculate segment maps and compose the full interval.
     *
     * @param samples Reference states in increasing path-length order, with consistent
     * transported frames, generated using the configured integration method. Passed by
     * value so clipping pre-roll cannot mutate a caller's retained sample vector.
     * @param startPathLength Requested entrance coordinate [m]. Earlier pre-roll is removed;
     * a bracketed entrance is numerically localized with the same ray tracker.
     * @return All nondegenerate segments through the final supplied sample and their product.
     * Fewer than two usable samples or no retained segments yields an empty result, not
     * an identity map. Every call restarts segment numbering at zero.
     *
     * Body transitions are detected between samples; the caller must resolve the orbit
     * finely enough not to skip complete bodies. This API does not sort a reversed orbit
     * or find a closed orbit. RF structures are rejected when encountered in selection;
     * failed ray-plane intersections and singular/ill-conditioned input solves also throw.
     * No maps or overlap flags on the elements are changed.
     */
    Result build(std::vector<ReferenceSample> samples, double startPathLength);

    /// Project the beamline x axis perpendicular to nonzero lab momentum, using its y axis
    /// as fallback. Returns axes only; other reference-state fields retain default values.
    static LinearTransferMapReference initialFrame(
            OpalBeamline& beamline, const Vector_t<double, 3>& momentum);
    /// Minimum-rotation transport to nonzero lab momentum; see LinearTransferMapReference.
    /// Copies all non-axis fields unchanged: the caller must separately update the ray state.
    static LinearTransferMapReference transportFrame(
            const LinearTransferMapReference& frame, const Vector_t<double, 3>& momentum);

private:
    using RayState = ExternalFieldRayTracker::State;
    OpalBeamline& itsOpalBeamline_m;
    ExternalFieldRayTracker tracker_m;
    double dt_m;
    Settings settings_m;

    /// Refine a nominal-body transition by tracked bisection; not a field-support boundary.
    LinearTransferMapReference refineBoundary(
            const std::shared_ptr<ElementBase>& element, const LinearTransferMapReference& before,
            const LinearTransferMapReference& after, bool entering);
    /// Track to the exit reference plane, allowing each ray its own arrival time and supports.
    RayState trackRayToExit(
            const RayState& initial, const LinearTransferMapReference& exit,
            double referenceFlightTime);
    LinearTransferMap makeLinearTransferMap(
            const LinearTransferMapReference& entrance, const LinearTransferMapReference& exit,
            std::size_t pass);
    /// Launch twelve rays and form D_out * inverse(D_in); condition is the unscaled input
    /// infinity-norm condition estimate, not an optics error estimate.
    matrix6x6_t makeCenteredMap(
            const LinearTransferMapReference& entrance, const LinearTransferMapReference& exit,
            const std::array<double, 6>& steps, double& condition);
    static std::array<double, 6> coordinates(
            const RayState& ray, const LinearTransferMapReference& reference);
    static RayState rayFromCoordinates(
            const std::array<double, 6>& coordinates, const LinearTransferMapReference& reference);
};
#endif
