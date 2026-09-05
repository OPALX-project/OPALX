// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_EXTERNAL_FIELD_RAY_TRACKER_H
#define OPAL_EXTERNAL_FIELD_RAY_TRACKER_H

#include <functional>
#include <limits>
#include <string>
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
 * The default is drift--Boris-kick--drift. Map computation can instead select
 * classical RK4 or the twelve-stage eighth-order DOP853 formula. No collective
 * fields or analytic element transfer matrices enter this calculation.
 *
 * Boris drifts are evaluated directly in metres:
 * \f$\mathbf r_{1/2}=\mathbf r_0+(ch/2)\mathbf u_0/\gamma_0\f$ and
 * \f$\mathbf r_1=\mathbf r_{1/2}+(ch/2)\mathbf u_1/\gamma_1\f$,
 * where \f$\gamma_i=\sqrt{1+|\mathbf u_i|^2}\f$. The signed path is advanced by
 * \f[
 * \Delta s=\frac{ch}{2}\left(\frac{|\mathbf u_0|}{\gamma_0}
 *                         +\frac{|\mathbf u_1|}{\gamma_1}\right).
 * \f]
 * This integrates speed, not the chord between endpoints. It is exact for constant
 * speed (up to momentum/roundoff error), and second-order quadrature with electric
 * acceleration. Position, time and path additions retain Kahan correction terms.
 * The Boris kick, its order and field-support tolerances are unchanged. Only this
 * host-side reference/map tracker is affected, not the production particle pusher.
 *
 * RK4 and DOP853 solve the same relativistic Lorentz equations, including path:
 * \f[
 * \dot{\mathbf r}=c\mathbf u/\gamma,\qquad
 * \dot{\mathbf u}=\frac{qc}{\mathcal E_0}
 *    \left(\mathbf E+\frac{c}{\gamma}\mathbf u\times\mathbf B\right),\qquad
 * \dot s=c|\mathbf u|/\gamma,
 * \f]
 * where \f$\mathcal E_0=mc^2\f$ is the numerical rest energy in eV and q is the
 * signed charge in elementary-charge units. Each stage samples its own position,
 * momentum and time. For \f$\mathbf y=(\mathbf r,\mathbf u,s)\f$,
 * \f[
 * k_i=h f(t_0+c_i h,\mathbf y_0+\sum_{j<i}a_{ij}k_j),\qquad
 * \mathbf y_1=\mathbf y_0+\sum_i b_i k_i.
 * \f]
 * Time is advanced exactly by h with compensation. Weighted increments and
 * position/path updates use compensation, but momentum and field evaluations
 * remain doubles. DOP853 coefficients are documented in RungeKuttaTableau.h.
 * Both RK methods use fixed nominal DT, NOT adaptive error control. DOP853's
 * embedded 5/3 estimators and dense output are not used. An independent h/2 solve
 * supplies a same-order diagnostic midpoint; it does not change the full-h result.
 * Including its field sample, a smooth trial costs 9 evaluations for RK4 and 25
 * for DOP853 (versus 1 for Boris). These methods are not symplectic and do not
 * preserve magnetic energy exactly; their smooth global orders are 4 and 8.
 *
 * advance() resolves field-support transitions independently for each ray. A trial
 * is accepted only when support sets at the start, every field-evaluation stage,
 * diagnostic midpoint and end agree. Otherwise it is bisected, advancing the first
 * half before recomputing the second. Transition-containing substeps terminate at
 * \f[
 * c|h|\leq\max(10^{-12}c|\Delta t|,
 *     64\epsilon_{\rm mach}\max(1\,{\rm m},|\mathbf r|)).
 * \f]
 * Thus a residual straddle is bounded by the event tolerance or floating-point
 * spatial resolution, not by the nominal time step. Smooth-field steps retain
 * the selected discretization order. This is not a symplecticity guarantee:
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
    /// Numerical ray stepping, independent of finite-difference/Richardson order.
    /// RK choices apply to maps only; fixed DT, no fallback for unknown methods.
    enum class IntegrationMethod { BORIS, RK4, DOP853 };
    static IntegrationMethod parseIntegrationMethod(const std::string& name);
    static std::string integrationMethodName(IntegrationMethod method);

    struct State {
        Vector_t<double, 3> position = Vector_t<double, 3>(0.0);
        Vector_t<double, 3> momentum = Vector_t<double, 3>(0.0);
        double time{0.0};
        double pathLength{0.0};
        /// Represented values are the displayed high parts minus these corrections.
        Vector_t<double, 3> positionCorrection = Vector_t<double, 3>(0.0);
        double timeCorrection{0.0};
        double pathLengthCorrection{0.0};
    };

    /// Return true for a material hit. Called at every RK stage or Boris drift midpoint.
    using FieldEvaluator =
            std::function<bool(const State&, Vector_t<double, 3>&, Vector_t<double, 3>&)>;
    struct Step {
        /// Boris: pre-kick drift midpoint. RK: independently integrated half-duration state.
        State midpoint;
        State end;
        Vector_t<double, 3> electric = Vector_t<double, 3>(0.0);
        Vector_t<double, 3> magnetic = Vector_t<double, 3>(0.0);
        double duration{0.0};
        bool hitMaterial{false};
    };

    ExternalFieldRayTracker(
            OpalBeamline& beamline, const PartData& reference,
            IntegrationMethod method = IntegrationMethod::BORIS);

    /// One numerical step with an explicit field evaluator (also supplies path-log data).
    Step step(const State& initial, double dt, const FieldEvaluator& fields) const;

    /// Advance through support events; optionally record accepted substeps in time order.
    State advance(const State& initial, double dt, std::vector<Step>* accepted = nullptr) const;

    /**
     * @brief Locate a signed path target within a time interval by tracked bisection.
     * Reintegrates each trial, including support subdivision; does not linearly
     * interpolate position/time under acceleration. The target must be bracketed
     * by initial and advance(initial, dt). Uses the same time-resolution floor as
     * support localization, and throws for an unbracketed target.
     */
    State advanceToPathLength(const State& initial, double dt, double target) const;

private:
    OpalBeamline& beamline_m;
    const PartData& reference_m;
    BorisPusher integrator_m;
    IntegrationMethod integrationMethod_m;
    double maximumStep_m{std::numeric_limits<double>::max()};
    Step borisStep(const State& initial, double dt, const FieldEvaluator& fields) const;
    Step rungeKuttaStep(const State& initial, double dt, const FieldEvaluator& fields) const;
    State advanceRecursive(
            const State& initial, double dt, double tolerance, std::vector<Step>* accepted,
            unsigned depth) const;
};
#endif
