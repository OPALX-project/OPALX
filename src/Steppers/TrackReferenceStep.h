// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPALX_TRACK_REFERENCE_STEP_H
#define OPALX_TRACK_REFERENCE_STEP_H

#include <cmath>
#include <functional>
#include <stdexcept>
#include "OPALTypes.h"
#include "Steppers/BorisPusher.h"

/** Internal host-side TRACK transport, independent of particle storage.
 * Positions are metres, momenta P/(mc), times seconds, mass rest energy [eV],
 * and charge in proton-charge units. The advance preserves scaled-position
 * drift/kick/drift arithmetic.
 */
namespace track_reference {
    struct State {
        Vector_t<double, 3> position, momentum;
        bool hitMaterial = false;
    };

    /// Field callback of a reference step: adds E (V/m) and B (T) at position r, momentum p and
    /// time t, and returns whether r is in material. record is true in one call per step, at the
    /// state half a step after the start, where a MONITOR records the reference
    /// (Monitor::recordReference()).
    using FieldEvaluator = std::function<bool(
            const Vector_t<double, 3>& r, const Vector_t<double, 3>& p, double t, bool record,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B)>;

    /// @throws std::invalid_argument unless dt is finite and nonzero and endTime is finite.
    inline void checkStep(double dt, double endTime) {
        if (!std::isfinite(dt) || dt == 0 || !std::isfinite(endTime))
            throw std::invalid_argument("Invalid TRACK reference step duration/time");
    }

    /** Advance a copy by one Boris step. The field callback receives the midpoint (R,P,t) with
     * record = true and returns a material flag. A trial callback must have no output or state
     * side effects. No energy or section projection is applied.
     */
    template <class Field>
    State advance(
            State state, double dt, double endTime, double mass, double charge, Field&& field) {
        checkStep(dt, endTime);
        const BorisPusher pusher;
        const double scale = Physics::c * dt;
        state.position /= scale;
        pusher.push(state.position, state.momentum, dt);
        state.position *= scale;
        Vector_t<double, 3> electric(0), magnetic(0);
        state.hitMaterial =
                field(state.position, state.momentum, endTime - 0.5 * dt, true, electric, magnetic);
        pusher.kick(state.position, state.momentum, electric, magnetic, dt, mass, charge);
        state.position /= scale;
        pusher.push(state.position, state.momentum, dt);
        state.position *= scale;
        return state;
    }

}  // namespace track_reference
#endif
