/** @file SpaceChargeTransport.h @brief Host NS-controlled analytic RMS Strang transport. */
#ifndef OPALX_SPACE_CHARGE_TRANSPORT_H
#define OPALX_SPACE_CHARGE_TRANSPORT_H

#include "Algorithms/LinearTransferMapBuilder.h"
#include "Algorithms/SpaceChargeMap.h"

namespace opalx::maps {
    /** @brief Slice resolution independent of DT, finite differences and quadrature. */
    struct SpaceChargeTransportSettings {
        std::size_t ns =
                1;  ///< Positive kick count per retained nominal body OR ownerless interval.
        std::size_t maximumSlices = 4096;  ///< Explicit resource limit for the complete plan.
        LinearTransferMapBuilder::Settings
                external;  ///< Existing FD/Richardson/integrator settings.
    };

    /** @brief Cached external transports and physical station data for one Strang interval. */
    struct SpaceChargeTransportSlice {
        std::size_t sourceSegment, subdivision;
        double pathBegin, pathMidpoint, pathEnd;  ///< Actual localized reference path stations [m].
        double integrationLength;  ///< Compensated endpoint difference h [m], not a chord.
        SpaceChargeCoordinates midpoint;
        LinearTransferMapDraft leftExternal, rightExternal;  ///< Actual finite-difference halves.
        tps::Map<double> left, right;  ///< Common-time canonical maps, with explicit chart bridges.
    };

    /** @brief Immutable external-field plan, reusable for different input covariances/charges.
     * Owns no element/beamline pointers. Rebuild after lattice/species/reference edits; this is a
     * snapshot, not a lattice fingerprint. Requires active Kokkos/serialized host use, and must
     * be destroyed before runtime finalization. No global map-storage registration is performed.
     */
    class SpaceChargeTransportPlan {
    public:
        const std::vector<SpaceChargeTransportSlice>& slices() const { return slices_m; }
        const SpaceChargeCoordinates& entrance() const { return entrance_m; }
        const SpaceChargeCoordinates& exit() const { return exit_m; }
        const std::string& id() const { return id_m; }
        std::size_t ns() const { return ns_m; }
        double pathLength() const { return pathLength_m; }  ///< Total SC integration length [m].

    private:
        friend SpaceChargeTransportPlan buildSpaceChargeTransportPlan(
                OpalBeamline&, const PartData&, const TransferMapPass&, double,
                const SpaceChargeTransportSettings&, const std::string&);
        SpaceChargeTransportPlan(
                SpaceChargeCoordinates entrance, SpaceChargeCoordinates exit,
                std::vector<SpaceChargeTransportSlice> slices, std::string id, std::size_t ns,
                double pathLength);
        SpaceChargeCoordinates entrance_m, exit_m;
        std::vector<SpaceChargeTransportSlice> slices_m;
        std::string id_m;
        std::size_t ns_m;
        double pathLength_m;
    };

    /**
     * @brief Build actual external half-maps over a retained forward reference pass.
     *
     * Each retained source interval of length L is divided into NS equal path intervals h=L/NS,
     * with stations entrance, midpoint and exit. An isolated full element normally corresponds
     * to one source interval; clipped bodies and ownerless drift intervals each receive NS.
     * This explicit initial policy avoids assuming a physical element appears once per pass.
     * Zero-length records, empty/revoked passes and overlapping-field flags are rejected.
     *
     * Interior stations are localized by actual reference integration with nominal DT and tracked
     * path bisection; no position/time/momentum interpolation. Recorded pass endpoints and join
     * frames are retained exactly. The caller supplies the same beamline/species/integrator as the
     * sampled pass, with adequate DT resolution. Method mismatch is rejected. Existing reference
     * discretization/closure errors are not corrected by NS refinement.
     *
     * Both external halves use the existing centered-difference/Richardson ray builder with their
     * own fields/support crossings. No square root or fractional power of a whole map is used.
     * Field sampling uses the tracker's reference-field evaluator. Nominal-body or field-support
     * overlaps and RF encountered at reference/shadow field evaluations throw: overlap scheduling
     * is DEFERRED. This remains sampled geometry; unresolved holes/grazing features require DT
     * refinement as in ExternalFieldRayTracker. Default legacy field summation is unchanged.
     *
     * @param beamline Placed external lattice, unchanged during plan construction.
     * @param reference Species rest energy [eV] and charge [e], as used for the source pass.
     * @param pass Complete retained external pass, including ownerless intervals.
     * @param dt Positive nominal external-ray time step [s], independent of NS.
     * @param settings Positive NS, plan-size limit and external numerical settings.
     * @param occurrenceId Nonempty namespace for half-map/kick provenance, combined with pass
     * ID/NS.
     * @return Value-owned external plan; no maps are attached to the runtime lattice.
     * @throws OpalException Invalid settings/intervals, missing path resolution, unsupported
     * overlap/RF/material, failed ray integration or FD solve, invalid charts or resource limits.
     */
    SpaceChargeTransportPlan buildSpaceChargeTransportPlan(
            OpalBeamline& beamline, const PartData& reference, const TransferMapPass& pass,
            double dt, const SpaceChargeTransportSettings& settings,
            const std::string& occurrenceId);

    /** @brief Density and kick diagnostics from one actual covariance propagation. */
    struct SpaceChargeTransportStep {
        SpaceChargeMatrix midpointCovariance, exitCovariance;  ///< Common-time canonical charts.
        SpaceChargeKick kick;  ///< Recomputed from this pass's midpoint covariance, never cached.
    };
    struct SpaceChargeTransportResult {
        SpaceChargeMatrix covariance;  ///< Common-time canonical covariance at plan exit.
        tps::Map<double>
                polynomial;  ///< Frozen-along-this-pass product, canonical entrance to exit.
        std::vector<SpaceChargeTransportStep> steps;
    };

    /**
     * @brief Propagate full 6D covariance with midpoint-density analytic RMS Strang kicks.
     *
     * For each h, column-vector ordering is M_right*K(h,Sigma_mid)*M_left:
     * \f[\Sigma_{mid}=M_L\Sigma_{in}M_L^T,\quad
     *     \Sigma_{kick}=K\Sigma_{mid}K^T,\quad
     *     \Sigma_{out}=M_R\Sigma_{kick}M_R^T.\f]
     * Every kick recomputes the 3D density/rest tensor from the transported full covariance;
     * neither transverse nor longitudinal widths are frozen across slices. During the pure
     * momentum kick positions do not change, so its spatial density is constant. No implicit
     * inner iteration is needed for that subflow. Charge is the signed physical total [C].
     *
     * The result polynomial holds the self-field coefficients along THIS covariance trajectory;
     * it is not the derivative of the nonlinear covariance evolution and is not valid for a
     * different input distribution without rerunning this function. External half-maps are cached
     * in the plan; density/kicks are not. This is a single pass, not a periodic matching solve.
     * For smooth sufficiently accurate external maps/quadrature, NS refinement is second order
     * for the closed RMS model; FD, reference tracking and model errors set separate floors.
     * Initial covariance must be SPD in plan.entrance()'s common-time canonical chart. All outputs
     * are local values; exceptions publish no partial result. Active Kokkos, serialized host.
     * @param plan Cached actual external half-maps and midpoint coordinate systems.
     * @param initialCanonicalCovariance Full centered 6D covariance at the plan entrance.
     * @param profile Uniform spatial ellipsoid or Gaussian RMS closure.
     * @param totalCharge Signed physical bunch charge [C], constant across all slices.
     * @param quadrature Dimensionless tensor integration controls, independent of NS.
     * @return Exit covariance, frozen-along-this-pass polynomial and slice diagnostics.
     * @test TestSpaceChargeTransport checks actual half-maps, midpoint density, independent
     * moment-ODE convergence, provenance, zero current, rejection paths and MPI consistency.
     */
    SpaceChargeTransportResult transportSpaceCharge(
            const SpaceChargeTransportPlan& plan,
            const SpaceChargeMatrix& initialCanonicalCovariance,
            opalx::spacecharge::analytic3d::Profile profile, double totalCharge,
            const opalx::spacecharge::analytic3d::IntegrationSettings& quadrature = {});
}  // namespace opalx::maps
#endif
