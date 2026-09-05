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
        /// Metres for x,y,zeta; dimensionless for x',y',delta, in that interleaved order.
        std::array<double, 6> finiteDifferenceSteps{1.e-3, 1.e-3, 1.e-3, 1.e-3, 1.e-3, 1.e-3};
        unsigned richardsonLevels{0};
        ExternalFieldRayTracker::IntegrationMethod integrationMethod{
                ExternalFieldRayTracker::IntegrationMethod::BORIS};
        /// Reject invalid levels/amplitudes and unsupported integration methods.
        void validate() const;
    };

    struct ReferenceSample {
        LinearTransferMapReference state;
    };
    struct Segment {
        LinearTransferMap map;
        IndexMap::value_t owners;
    };
    struct Result {
        /// Includes unowned intervals (possibly containing fringe fields); each overlap once.
        std::vector<Segment> segments;
        std::optional<matrix6x6_t> combined;
        std::optional<double> determinantResidual;
        std::optional<double> symplecticResidual;
    };

    LinearTransferMapBuilder(OpalBeamline& beamline, const PartData& reference, double dt);
    LinearTransferMapBuilder(
            OpalBeamline& beamline, const PartData& reference, double dt, const Settings& settings);
    /// Samples are copied so pre-roll clipping cannot mutate the caller's reference orbit.
    Result build(std::vector<ReferenceSample> samples, double startPathLength);

    /// Initialize the transverse axes from the beamline frame, without a roll singularity.
    static LinearTransferMapReference initialFrame(
            OpalBeamline& beamline, const Vector_t<double, 3>& momentum);
    /// Minimum-rotation transport of an existing orthonormal frame to a new tangent.
    static LinearTransferMapReference transportFrame(
            const LinearTransferMapReference& frame, const Vector_t<double, 3>& momentum);

private:
    using RayState = ExternalFieldRayTracker::State;
    OpalBeamline& itsOpalBeamline_m;
    ExternalFieldRayTracker tracker_m;
    double dt_m;
    Settings settings_m;

    LinearTransferMapReference refineBoundary(
            const std::shared_ptr<ElementBase>& element, const LinearTransferMapReference& before,
            const LinearTransferMapReference& after, bool entering);
    RayState trackRayToExit(
            const RayState& initial, const LinearTransferMapReference& exit,
            double referenceFlightTime);
    LinearTransferMap makeLinearTransferMap(
            const LinearTransferMapReference& entrance, const LinearTransferMapReference& exit,
            std::size_t pass);
    matrix6x6_t makeCenteredMap(
            const LinearTransferMapReference& entrance, const LinearTransferMapReference& exit,
            const std::array<double, 6>& steps, double& condition);
    static std::array<double, 6> coordinates(
            const RayState& ray, const LinearTransferMapReference& reference);
    static RayState rayFromCoordinates(
            const std::array<double, 6>& coordinates, const LinearTransferMapReference& reference);
};
#endif
