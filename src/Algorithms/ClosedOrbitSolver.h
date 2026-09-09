// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_CLOSED_ORBIT_SOLVER_H
#define OPAL_CLOSED_ORBIT_SOLVER_H

#include <string>
#include "Algorithms/OneTurnMap.h"

/** Damped Newton fixed-point solver for a deterministic four-dimensional map.
 * With F(u)=T(u)-u and M=dT/du, solve (M-I) delta=-F. For diagonal coordinate
 * scales S, column-pivoted Householder QR solves S^-1(M-I)S v=-S^-1 F;
 * delta=S v. No matrix inverse, collective operation or tracking state mutation
 * is used. This host-only solver does not change the pusher or its conservation.
 * The callback must use the same fixed-energy section for every evaluation.
 * Coordinates are (x,px,y,py), metres and mechanical p/(mc), not slopes.
 *
 * Central differences have O(h^2) truncation and O(roundoff/h) noise. Damping
 * uses an Armijo decrease in ||S^-1 F||_2. Componentwise physical residual AND
 * undamped Newton correction must meet tolerances; one fresh traversal then
 * verifies residual closure. This proves numerical closure only, not timestep
 * convergence, uniqueness, stability or a physically accurate closed orbit.
 *
 * Callback OpalException losses during line search reject a trial; failures at
 * the initial state, derivative probes or final verification return diagnostic
 * failure statuses. Other exception types propagate (programming/system errors).
 */
class ClosedOrbitSolver {
public:
    using Coordinates = OneTurnMap::Coordinates;
    using Matrix      = OneTurnMap::Matrix;
    using Map         = std::function<Coordinates(const Coordinates&)>;
    enum class Status {
        Converged,
        SingularJacobian,
        LineSearchFailed,
        MaxIterations,
        EvaluationFailed,
        VerificationFailed
    };
    struct Settings {
        Coordinates scales{1, 1, 1, 1};  ///< Characteristic scales in coordinate units.
        Coordinates finiteDifferenceSteps{1e-6, 1e-6, 1e-6, 1e-6};  ///< Positive absolute h.
        double positionTolerance = 1e-10;                           ///< Metres, x and y.
        double momentumTolerance = 1e-10;                           ///< p/(mc), px and py.
        double rankTolerance     = 1e-12;  ///< Relative pivot cutoff versus first QR pivot.
        unsigned maxIterations   = 20;     ///< Maximum accepted corrections.
        unsigned maxBacktracks   = 20;     ///< Trials at lambda=1,1/2,...,2^-maxBacktracks.
        bool damping             = true;
    };
    struct Iteration {
        Coordinates coordinates{}, residual{}, correction{};  ///< Accepted state and step.
        double damping    = 0;
        double pivotRatio = 0;  ///< min/max |diag(R)| before step; NOT a condition number.
    };
    struct Result {
        Status status = Status::EvaluationFailed;
        Coordinates coordinates{}, residual{};
        Matrix matrix{};  ///< dT/du at last linearized state (final state on success).
        double pivotRatio    = 0;
        unsigned evaluations = 0;  ///< Includes derivative probes and failed trials.
        std::vector<Iteration> iterations;
        std::string message;
    };
    /// Invalid settings/initial coordinates throw OpalException before tracking.
    static Result solve(const Map&, const Coordinates&, const Settings&);
    /// Adapter to the general shared return-map service; does not modify it.
    static Result solve(const OneTurnMap&, const Coordinates&, const Settings&);
};
#endif
