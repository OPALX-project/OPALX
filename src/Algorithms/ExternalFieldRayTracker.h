// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_EXTERNAL_FIELD_RAY_TRACKER_H
#define OPAL_EXTERNAL_FIELD_RAY_TRACKER_H

#include <functional>
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
        bool hitMaterial{false};
    };

    ExternalFieldRayTracker(OpalBeamline& beamline, const PartData& reference);

    /// One numerical step with an explicit field evaluator (also supplies path-log data).
    Step step(const State& initial, double dt, const FieldEvaluator& fields) const;

    /// Advance with the superposed external fields selected at the ray's position.
    State advance(const State& initial, double dt) const;

private:
    OpalBeamline& beamline_m;
    const PartData& reference_m;
    BorisPusher integrator_m;
};
#endif
