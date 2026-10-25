// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPALX_TRACK_REFERENCE_STEP_H
#define OPALX_TRACK_REFERENCE_STEP_H

#include <mpi.h>
#include <Kokkos_NumericTraits.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include "OPALTypes.h"
#include "Steppers/BorisPusher.h"

class OpalBeamline;
class PartData;

/** Internal thin-gap event ordering shared by host reference and device bunch transport.
 * Positions are laboratory metres, momenta p/(mc), time seconds, and mass rest
 * energy [eV]. The policy supplies magnetic transport, gap geometry/kicks, visited
 * storage, and an accepted-kick callback. Trials never invoke kick or accepted.
 * Magnetic arithmetic and coordinate transforms remain owned by each caller.
 */
namespace rf_event {
    using Vector = Vector_t<double, 3>;
    enum class Failure { None = 0, FieldSupport = 1, GapKick = 2, InvalidStep = 3 };
    struct Diagnostic {
        Failure failure = Failure::None;
        int gap         = -1;
        Vector position = Vector(0), momentum = Vector(0);
        double time = 0;
        KOKKOS_INLINE_FUNCTION void set(
                Failure f, const Vector& r, const Vector& p, double t, int g = -1) {
            failure  = f;
            position = r;
            momentum = p;
            time     = t;
            gap      = g;
        }
    };
    struct Result {
        bool success;
        double duration;
    };

    /** Localize positive plane crossings with 40 bisections, then apply a full kick.
     * An initial point on the plane does not cross; equal roots select the later
     * policy index. Each accepted gap is visited once per call. accepted() may
     * stop after a complete kick (reference EKINSTOP); its duration is dt-remaining.
     * Caller validation and rollback are deliberately separate: this routine can
     * modify r/p before a later failure. Policies must keep trial transport free
     * of irreversible side effects. The algorithm owns no allocations or
     * collectives; host policies may allocate their visited storage.
     *
     * terminalEndpoint opts into a right-closed final RF event after the reference
     * has met EKINSTOP. Independent host/device/frame arithmetic can leave a
     * coincident particle just before that plane. Only a clearly upstream,
     * approaching particle whose final signed distance is within 64*N*epsilon*L
     * receives the full endpoint kick. Here N=max(1,frameSteps) counts accepted
     * frame updates including this final step, and L=max(1 m, |r_d|, |endR_d|)
     * over laboratory components d. Each step performs drift, rigid-frame round
     * trips and plane subtraction; the factor 64 conservatively budgets their
     * accumulated arithmetic error to first order, assuming 64*N*epsilon << 1
     * (a nonfinite factor or factor >=1 is rejected). Arrivals within this
     * allowance are numerically unresolved and receive a complete kick by the
     * terminal convention; this can include a physically later particle within
     * the allowance. Crossings beyond the allowance are excluded.
     * This allowance describes numerical endpoint resolution, not magnetic/PIC
     * discretization error or a bound on integration accuracy. It grows with the
     * number of frame updates even if dt shrinks. Ordinary steps keep strict
     * crossing signs. The caller must stop after this terminal step; no position
     * projection, extra drift, fractional impulse or persistent visit state is used.
     */
    template <class Policy>
    KOKKOS_INLINE_FUNCTION Result
    advance(Vector& r, Vector& p, double time, double dt, double mass, Policy& policy,
            Diagnostic* diagnostic = nullptr, bool terminalEndpoint = false,
            unsigned long long frameSteps = 1) {
        const double accumulatedRoundoff =
                64 * Kokkos::epsilon_v<double> * static_cast<double>(frameSteps ? frameSteps : 1);
        if (terminalEndpoint
            && (!Kokkos::isfinite(accumulatedRoundoff) || !(accumulatedRoundoff < 1))) {
            if (diagnostic) diagnostic->set(Failure::InvalidStep, r, p, time);
            return {false, 0};
        }
        double remaining = dt;
        while (remaining > 0) {
            Vector endR = r, endP = p;
            if (!policy.magnetic(endR, endP, remaining, time, mass)) {
                if (diagnostic) diagnostic->set(Failure::FieldSupport, r, p, time);
                return {false, dt - remaining};
            }
            double earliest = remaining;
            int selected    = -1;
            for (size_t j = 0; j < policy.count(); ++j) {
                if (policy.visited(j)) continue;
                const double startDistance = policy.distance(j, r);
                const double endDistance   = policy.distance(j, endR);
                double scale               = 1.;
                if (terminalEndpoint)
                    for (unsigned d = 0; d < 3; ++d)
                        scale = Kokkos::fmax(
                                scale, Kokkos::fmax(Kokkos::abs(r(d)), Kokkos::abs(endR(d))));
                const double roundoff      = terminalEndpoint ? accumulatedRoundoff * scale : 0;
                const bool roundedEndpoint = terminalEndpoint && startDistance < -roundoff
                                             && endDistance < 0 && endDistance >= -roundoff
                                             && endDistance > startDistance;
                if (!(startDistance < 0) || (endDistance < 0 && !roundedEndpoint)) continue;
                double low = 0, high = remaining;
                // A roundoff-short endpoint is accepted at the common final time;
                // a real bracket retains the existing forty-iteration localization.
                for (int iteration = 0; iteration < (roundedEndpoint ? 0 : 40); ++iteration) {
                    const double mid = 0.5 * (low + high);
                    Vector trialR = r, trialP = p;
                    if (!policy.magnetic(trialR, trialP, mid, time, mass)) {
                        if (diagnostic) diagnostic->set(Failure::FieldSupport, r, p, time, j);
                        return {false, dt - remaining};
                    }
                    if (policy.distance(j, trialR) < 0)
                        low = mid;
                    else
                        high = mid;
                }
                Vector trialR = r, trialP = p;
                if (!policy.magnetic(trialR, trialP, high, time, mass)) {
                    if (diagnostic) diagnostic->set(Failure::FieldSupport, r, p, time, j);
                    return {false, dt - remaining};
                }
                if (policy.supports(j, trialR) && high <= earliest) {
                    earliest = high;
                    selected = j;
                }
            }
            if (selected < 0) {
                r = endR;
                p = endP;
                return {true, dt};
            }
            if (!policy.magnetic(r, p, earliest, time, mass)) {
                if (diagnostic) diagnostic->set(Failure::FieldSupport, r, p, time, selected);
                return {false, dt - remaining};
            }
            time += earliest;
            remaining -= earliest;
            if (!policy.kick(selected, r, p, time, mass)) {
                if (diagnostic) diagnostic->set(Failure::GapKick, r, p, time, selected);
                return {false, dt - remaining};
            }
            policy.markVisited(selected);
            if (policy.accepted(selected, r, p, time)) return {true, dt - remaining};
        }
        return {true, dt};
    }

    /** Host adapter for element virtual calls and unrestricted visited storage.
     * The bridge is never a device policy. Guarded bodies let CUDA compile the
     * shared template while std::function and host elements stay on the host.
     */
    struct HostPolicy {
        size_t number = 0;
        std::function<bool(Vector&, Vector&, double, double, double)> magneticAdvance;
        std::function<double(size_t, const Vector&)> planeDistance;
        std::function<bool(size_t, const Vector&)> inProfile;
        std::function<bool(size_t, const Vector&, Vector&, double, double)> applyKick;
        std::function<bool(size_t)> wasVisited;
        std::function<void(size_t)> mark;
        std::function<bool(size_t, const Vector&, const Vector&, double)> onAccepted;

        KOKKOS_INLINE_FUNCTION size_t count() const { return number; }
        KOKKOS_INLINE_FUNCTION bool magnetic(
                Vector& r, Vector& p, double dt, double t, double mass) const {
            KOKKOS_IF_ON_HOST((return magneticAdvance(r, p, dt, t, mass);))
            KOKKOS_IF_ON_DEVICE((Kokkos::abort("RF host policy called on device");))
            return false;
        }
        KOKKOS_INLINE_FUNCTION double distance(size_t j, const Vector& r) const {
            KOKKOS_IF_ON_HOST((return planeDistance(j, r);))
            KOKKOS_IF_ON_DEVICE((Kokkos::abort("RF host policy called on device");))
            return 0;
        }
        KOKKOS_INLINE_FUNCTION bool supports(size_t j, const Vector& r) const {
            KOKKOS_IF_ON_HOST((return inProfile(j, r);))
            KOKKOS_IF_ON_DEVICE((Kokkos::abort("RF host policy called on device");))
            return false;
        }
        KOKKOS_INLINE_FUNCTION bool kick(
                size_t j, const Vector& r, Vector& p, double t, double mass) const {
            KOKKOS_IF_ON_HOST((return applyKick(j, r, p, t, mass);))
            KOKKOS_IF_ON_DEVICE((Kokkos::abort("RF host policy called on device");))
            return false;
        }
        KOKKOS_INLINE_FUNCTION bool visited(size_t j) const {
            KOKKOS_IF_ON_HOST((return wasVisited(j);))
            KOKKOS_IF_ON_DEVICE((Kokkos::abort("RF host policy called on device");))
            return false;
        }
        KOKKOS_INLINE_FUNCTION void markVisited(size_t j) const {
            KOKKOS_IF_ON_HOST((mark(j);))
            KOKKOS_IF_ON_DEVICE((Kokkos::abort("RF host policy called on device");))
        }
        KOKKOS_INLINE_FUNCTION bool accepted(
                size_t j, const Vector& r, const Vector& p, double t) const {
            KOKKOS_IF_ON_HOST((return onAccepted(j, r, p, t);))
            KOKKOS_IF_ON_DEVICE((Kokkos::abort("RF host policy called on device");))
            return false;
        }
    };
}  // namespace rf_event

/** Internal host-side TRACK transport, independent of particle storage.
 * Positions are metres, momenta P/(mc), times seconds, mass rest energy [eV],
 * and charge in proton-charge units. The legacy advance preserves scaled-position
 * drift/kick/drift arithmetic; the resolved variant shares COF's spatial solver.
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

    /** Shared field sampling for committed TRACK steps and terminal-event trials.
     * diagnostics=false omits passive monitors; supported magnetic elements use
     * their normal reference field and aperture predicates. diagnostics=true
     * retains the existing monitor sampling, performed only for the committed step.
     */
    State advanceInBeamline(
            OpalBeamline&, const PartData&, State, double dt, double endTime, bool diagnostics);

    /** Spatial, support-boundary-resolved Boris transport for static, external-field
     * TURNS reference tracking. Physical particles use DeviceExternalField.
     * Units match State; no closed-orbit projection or momentum renormalization.
     * This reference path runs on the host. Passive monitors are sampled only
     * on accepted reference substeps, never on particle or terminal-event trials.
     */
    State advanceResolvedInBeamline(
            OpalBeamline&, const PartData&, State, double dt, double endTime, bool diagnostics);

    /** Run a side-effect-free event search on rank zero and broadcast its duration
     * or diagnostic before any particle kernels/collectives execute. All ranks
     * must call together; failures become OpalException on every rank. If given,
     * onFailure runs on every rank after the diagnostic broadcast and before
     * throwing. This allows collective flushing of already accepted probe hits
     * when a terminal reference trial fails, without recording trial states.
     */
    double collectiveTerminalStep(MPI_Comm, const std::function<double()>& search);

    /// Internal overload for collective accepted-output cleanup after a failed trial.
    double collectiveTerminalStep(
            MPI_Comm, const std::function<double()>& search,
            const std::function<void()>& onFailure);
}  // namespace track_reference
#endif
