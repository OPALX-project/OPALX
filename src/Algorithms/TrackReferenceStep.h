// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPALX_TRACK_REFERENCE_STEP_H
#define OPALX_TRACK_REFERENCE_STEP_H

#include <cmath>
#include <stdexcept>
#include "OPALTypes.h"
#include "Steppers/BorisPusher.h"

class OpalBeamline;
class PartData;

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

    /** Advance a copy. The field callback receives the midpoint (R,P,t,E,B)
     * and returns a material flag. A trial callback must have no output or state
     * side effects. No energy or section projection is applied.
     */
    template <class Field>
    State advance(
            State state, double dt, double endTime, double mass, double charge, Field&& field) {
        if (!std::isfinite(dt) || dt == 0 || !std::isfinite(endTime))
            throw std::invalid_argument("Invalid TRACK reference step duration/time");
        const BorisPusher pusher;
        const double scale = Physics::c * dt;
        state.position /= scale;
        pusher.push(state.position, state.momentum, dt);
        state.position *= scale;
        Vector_t<double, 3> electric(0), magnetic(0);
        state.hitMaterial =
                field(state.position, state.momentum, endTime - 0.5 * dt, electric, magnetic);
        pusher.kick(state.position, state.momentum, electric, magnetic, dt, mass, charge);
        state.position /= scale;
        pusher.push(state.position, state.momentum, dt);
        state.position *= scale;
        return state;
    }

    /** One reference step in the fields of the beamline (OpalBeamline::applyOnHost()).
     * diagnostics=true also lets the monitors record the reference (committed steps only);
     * diagnostics=false is for side-effect-free trials.
     */
    State advanceInBeamline(
            OpalBeamline&, const PartData&, State, double dt, double endTime, bool diagnostics);

}  // namespace track_reference
#endif
