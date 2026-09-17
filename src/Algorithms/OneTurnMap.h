// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_ONE_TURN_MAP_H
#define OPAL_ONE_TURN_MAP_H

#include <array>
#include <functional>
#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/ExternalFieldRayTracker.h"

/** Fixed-section, fixed-momentum Poincare map for single-particle ring optics.
 * Coordinates are (x, px, y, py), in metres and mechanical p/(mc), NOT slopes.
 * Every launch has local z=0 and pz=sqrt(p0^2-px^2-py^2)>0. The section's
 * local +z is the required return direction, and is fixed for all perturbations.
 * A negative-side excursion arms a return; the launch and reverse crossing do
 * not count. The crossing is located by reintegration, not linear interpolation.
 *
 * This is a tracking service, not a closed-orbit solver or topology validator.
 * The caller must supply an initialized static magnetic RING, valid geometry,
 * mass/charge and integrator, and sufficiently fine DT to resolve crossings.
 * Field-free intervals are permitted. A loss propagates from the shared tracker.
 * A general self-intersecting orbit may need a more selective section: this
 * service returns the first armed positive crossing within the search budget.
 *
 * Host-only; no bunch, field solver, MPI collective, or output is constructed.
 * No energy projection is applied to returned momentum. Integration drift must
 * be checked separately from the nonlinear fixed-point residual.
 */
class OneTurnMap {
public:
    using Coordinates = std::array<double, 4>;
    using Matrix      = std::array<Coordinates, 4>;
    using State       = ExternalFieldRayTracker::State;
    using Step        = ExternalFieldRayTracker::Step;
    /// Must advance through support events, throwing on loss. Accepted substeps
    /// must be nonempty, contiguous, in time order, and cover the requested dt.
    using Advance = std::function<State(const State&, double, std::vector<Step>*)>;

    struct Settings {
        double momentum    = 0;       ///< Fixed positive total p/(mc).
        double dt          = 0;       ///< Nominal step [s], positive.
        double time        = 0;       ///< Initial time [s].
        unsigned maxSteps  = 200000;  ///< Hard nominal-step budget per traversal.
        double maxPath     = 0;       ///< Positive search limit [m], e.g. twice circumference.
        double armDistance = 1e-9;    ///< Negative-side hysteresis [m].
    };
    struct Result {
        Coordinates coordinates{};
        State ray;
        unsigned steps         = 0;
        double sectionResidual = 0;  ///< Signed local z at return [m].
    };

    /// The tracker and its beamline/reference must outlive this service.
    OneTurnMap(const ExternalFieldRayTracker&, const CoordinateSystemTrafo&, const Settings&);
    /// Injection boundary for analytic tests or another validated shared stepper.
    OneTurnMap(Advance, const CoordinateSystemTrafo&, const Settings&);
    Result operator()(const Coordinates&) const;
    /// Eight independent returns at fixed energy and in the identical section.
    /// Steps have units (m, p/(mc), m, p/(mc)); all must be finite and positive.
    Matrix jacobian(const Coordinates&, const Coordinates& steps) const;

private:
    Advance advance_m;
    CoordinateSystemTrafo section_m;
    Settings settings_m;
    double distance(const State&) const;
    State locate(const State&, double duration) const;
};
#endif
