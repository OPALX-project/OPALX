// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPALX_BORIS_MIDPOINT_H
#define OPALX_BORIS_MIDPOINT_H

#include "OPALTypes.h"
#include "Steppers/BorisPusher.h"

namespace boris_midpoint {
    using Vector = Vector_t<double, 3>;

    /// A value-only trial state; particle views are never stored or modified here.
    struct State {
        Vector position, momentum;
    };

    /**
     * @brief Half drift in SI coordinates: r + c*dt*p/(2*sqrt(1+p.p)).
     *
     * Position is in metres, normalized momentum is p/(mc), and dt is in seconds.
     * Using -dt rolls back an uncommitted first half drift, to floating-point
     * accuracy, only when momentum has not changed. Any intervening field solve
     * must return the same physical coordinates in the same tracker frame; a
     * particle boundary-condition wrap cannot be undone by this operation.
     */
    KOKKOS_INLINE_FUNCTION Vector halfDrift(Vector r, const Vector& p, double dt) {
        r += (0.5 * Physics::c * dt / Kokkos::sqrt(1.0 + dot(p, p))) * p;
        return r;
    }

    /**
     * @brief Trial the existing Boris kick and second half drift from a true midpoint.
     *
     * The caller gathers space charge at the actual collective half-drift
     * positions and adds external fields in the same axes before this call.
     * E is in V/m, B in tesla, mass is rest energy in eV, and charge is in proton
     * charge units. All inputs are values/read-only; a rejected trial leaves the
     * live midpoint momentum unchanged. Accepting a trial does not require a new
     * transport path: the tracker retains its existing kick/push operations.
     */
    KOKKOS_INLINE_FUNCTION State endpoint(
            Vector midpoint, Vector p, const Vector& e, const Vector& b,
            double dt, double mass, double charge) {
        BorisPusher().kick(midpoint, p, e, b, dt, mass, charge);
        return {halfDrift(midpoint, p, dt), p};
    }
} // namespace boris_midpoint

#endif
