// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_EXTERNAL_FIELD_RAY_TRACKER_H
#define OPAL_EXTERNAL_FIELD_RAY_TRACKER_H

#include <functional>
#include <limits>
#include <vector>
#include "OPALTypes.h"
#include "Steppers/BorisPusher.h"

class OpalBeamline;
class PartData;

/**
 * @brief Host-side external-field integration shared by reference and transfer-map rays.
 *
 * Positions are in metres, time and step size in seconds, and mechanical momentum is
 * \f$\mathbf u=\mathbf p/(mc)=\vec\beta\gamma\f$. Fields are in V/m and tesla;
 * PartData supplies rest energy in eV and charge in elementary-charge units.
 * A step is the same drift--Boris-kick--drift used by OrbitThreader. No collective
 * fields or analytic element transfer matrices enter this calculation.
 *
 * advance() resolves field-support transitions independently for each ray. A trial
 * Boris step of duration h is accepted only when the support sets at its start,
 * drift midpoint and end agree. Otherwise it is bisected, advancing the first
 * half before recomputing the second. Transition-containing substeps terminate at
 * \f[
 * c|h|\leq\max(10^{-12}c|\Delta t|,
 *     64\epsilon_{\rm mach}\max(1\,{\rm m},|\mathbf r|)).
 * \f]
 * Thus a residual straddle is bounded by the event tolerance or floating-point
 * spatial resolution, not by the nominal time step. Smooth-field steps retain
 * the second-order Boris discretization. This is not a symplecticity guarantee:
 * the subdivision depends on the ray and the fields may be interpolated.
 *
 * The trial step is also capped at one quarter of the shortest positive longitudinal
 * field-support extent divided by c, preventing a regular longitudinal crossing
 * from skipping a complete thin element. The element's isInside() query, including
 * its bend chart, is authoritative. Grazing crossings, unresolved transverse holes
 * or discontinuities hidden inside a field table still require adequate nominal
 * resolution; support queries are not a field-interpolation error estimator.
 */
class ExternalFieldRayTracker {
public:
    struct State {
        Vector_t<double, 3> position = Vector_t<double, 3>(0.0);
        Vector_t<double, 3> momentum = Vector_t<double, 3>(0.0);
        double time{0.0};
    };

    /// Return true for a material hit. The midpoint state is before the kick.
    using FieldEvaluator =
            std::function<bool(const State&, Vector_t<double, 3>&, Vector_t<double, 3>&)>;
    struct Step {
        State midpoint;
        State end;
        Vector_t<double, 3> electric = Vector_t<double, 3>(0.0);
        Vector_t<double, 3> magnetic = Vector_t<double, 3>(0.0);
        double duration{0.0};
        bool hitMaterial{false};
    };

    ExternalFieldRayTracker(OpalBeamline& beamline, const PartData& reference);

    /// One numerical step with an explicit field evaluator (also supplies path-log data).
    Step step(const State& initial, double dt, const FieldEvaluator& fields) const;

    /// Advance through support events; optionally record accepted substeps in time order.
    State advance(const State& initial, double dt, std::vector<Step>* accepted = nullptr) const;

private:
    OpalBeamline& beamline_m;
    const PartData& reference_m;
    BorisPusher integrator_m;
    double maximumStep_m{std::numeric_limits<double>::max()};
    State advanceRecursive(
            const State& initial, double dt, double tolerance, std::vector<Step>* accepted,
            unsigned depth) const;
};
#endif
