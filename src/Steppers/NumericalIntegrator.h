//
// Class NumericalIntegrator
//   The particle part of one time step, as an exchangeable integration scheme. Each scheme
//   family has a file of its own next to this one, e.g. BorisIntegrators.h.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPALX_NUMERICAL_INTEGRATOR_H
#define OPALX_NUMERICAL_INTEGRATOR_H

#include "Distribution/EmittingSamplers.h"
#include "PartBunch/PartBunchFwd.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "Utility/IpplTimings.h"

#include <memory>
#include <string>

class OpalBeamline;

namespace opalx::spacecharge {
    class SpaceChargeSolver;
}  // namespace opalx::spacecharge

/**
 * @brief The order of the particle part of one time step: drifts, kicks, the self-field solve,
 *        emission and external fields.
 *
 * ParallelTracker calls advance() once per step and owns the rest of the step: the bunch time,
 * the reference particles, losses, output, stop rules and checkpoints. Between the operations all
 * particle data are in the reference frame; frame changes happen and close inside
 * solveSelfFields() (the solve frame) and addExternalFields() (the element frames).
 */
class NumericalIntegrator {
public:
    /// What the operations of a step act on. All borrowed; they must outlive the integrator.
    struct Context {
        PartBunch_t& bunch;
        OpalBeamline& beamline;
        opalx::spacecharge::SpaceChargeSolver& spaceChargeSolver;
        EmittingSamplers& emittingSamplers;
        /// Source plane of applySourcePlaneLoss().
        opalx::spacecharge::DirichletPlaneConfig dirichletPlane;
        int minBinEmitted        = 0;      ///< OPTION MINBINEMITTED.
        bool aggressiveStateSync = false;  ///< OPTION AGGRESSIVE_STATE_SYNC.
    };

    virtual ~NumericalIntegrator() = default;

    /// @brief Advance all active particles from the bunch time t to t + dt.
    virtual void advance() = 0;

protected:
    explicit NumericalIntegrator(Context context);

    /// @brief Drift every particle over half its own dt (PartBunch::pushParticles()).
    void halfDrift();

    /// @brief Momentum kick of every particle over its own dt with the E and B it holds, between
    ///        two spin half rotations (PartBunch::evolveSpinHalfStep()).
    void kick();

    /**
     * @brief Replace E and B of the primary container with its self field at the current
     *        positions and zero those of the other containers (PartBunch::resetFields(), then
     *        SpaceChargeSolver::solveStep()). Collective.
     *
     * Elements add to E and B, so this must precede addExternalFields() of the same kick.
     * Particles migrate between ranks only inside this operation.
     */
    void solveSelfFields();

    /**
     * @brief Append the particles born in this step in state @p start, each with its fractional
     *        dt (EmittingSamplers::emit()).
     *
     * A new particle carries no self field in its creation step: the scheme calls this after
     * solveSelfFields(), and new particles start with zero E and B.
     */
    void emit(EmissionStart start);

    /**
     * @brief Delete the particles moving backwards behind the source plane
     *        (OpalBeamline::markBackwardParticlesAtSourcePlane()). Collective.
     *
     * Space-charge kicks can push particles behind the cathode, and the loss cannot wait for the
     * cleanup after the step: a scheme calls it before it adds the external fields, so these
     * particles get no further kick.
     */
    void applySourcePlaneLoss();

    /**
     * @brief Add the fields of the elements to E and B (OpalBeamline::applyExternalFields()).
     *        Collective.
     * @param fieldFraction Elements sample their fields at t + fieldFraction dt.
     * @param record Whether MONITORs record. Monitor::apply() assumes the positions half a step
     *        after the step start.
     */
    void addExternalFields(double fieldFraction = 0.5, bool record = true);

private:
    Context context_m;

    /// Created in the step order of BORIS; they do not overlap.
    IpplTimings::TimerRef pushTimer_m;
    IpplTimings::TimerRef spaceChargeTimer_m;
    IpplTimings::TimerRef emissionTimer_m;
    IpplTimings::TimerRef externalFieldsTimer_m;
    IpplTimings::TimerRef spinTimer_m;
    IpplTimings::TimerRef kickTimer_m;
};

/**
 * @brief The integration scheme called @p name, acting on @p context: "BORIS" is BorisMidpoint
 *        (the default), "BORIS-PRESTEP" is BorisPrestep.
 * @throws OpalException for an unknown name.
 */
std::unique_ptr<NumericalIntegrator> makeNumericalIntegrator(
        const std::string& name, NumericalIntegrator::Context context);

#endif  // OPALX_NUMERICAL_INTEGRATOR_H
