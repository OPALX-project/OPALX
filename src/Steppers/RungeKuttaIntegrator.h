//
// Class RungeKutta4
//   Classical fourth-order Runge-Kutta for the particles of the bunch, with the self field solved
//   once per step and held through the stages.
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
#ifndef OPALX_RUNGE_KUTTA_INTEGRATOR_H
#define OPALX_RUNGE_KUTTA_INTEGRATOR_H

#include "Steppers/NumericalIntegrator.h"
#include "VectorMath.h"

#include <cstddef>
#include <vector>

/**
 * @brief Classical fourth-order Runge-Kutta for y = (R, P) of every particle (PhaseSpaceRate),
 *        with the external fields at every stage and the self field solved once per step and
 *        held. Scheme "RK4".
 *
 * One step: applySourcePlaneLoss() and solveSelfFields() at the start positions, emit() at the
 * birth state, then the four stages of runge_kutta::rk4. Stage i adds the element fields at its
 * positions and at the time t + c_i dt to the held self field and advances every particle over
 * its own dt. Like ASTRA, which treats the space-charge field like an external field taken from
 * a stored solve, the scheme is of fourth order in the external fields and of first order in
 * the self-field term (the self-field timing of BORIS-PRESTEP). MONITORs record at stage 1 only
 * (Monitor::apply()). Spin is not transported; validate() refuses a polarized beam.
 *
 * A particle born in the step goes through the four stages over its fractional dt, but the
 * elements sample its fields at the stage times of the bunch, not at its own; that costs the
 * fourth order in its creation step, an approximation of the kind the Boris schemes make.
 *
 * Stage 0 samples the element fields at the start positions. A particle a rounding behind the
 * start of a field region, e.g. the cathode plane at the entrance of an RF gun map, gets no field
 * in k_0, where the Boris schemes sample half a drift later.
 *
 * advanceReference() takes the same stages for one particle in the external fields, so a particle
 * on its reference particle stays on it.
 */
class RungeKutta4 final : public NumericalIntegrator {
public:
    explicit RungeKutta4(Context context);

    void advance() override;

    /// The stages of advanceStage() on the host, at the times endTime - dt + c_i dt; a stage in
    /// material marks the step.
    track_reference::State advanceReference(
            const track_reference::State& start, double dt, double endTime, double mass,
            double charge, const track_reference::FieldEvaluator& field) const override;

    // The kernels of advance(), public only because CUDA requires the function that encloses a
    // device lambda to be public.

    /// Copy E and B of the active containers, the self field after solveSelfFields() and
    /// emit(), into the scratch, which is sized for the current local particle counts.
    void holdSelfField();
    /// Copy the held self field back into E and B, before the element fields of a stage.
    void restoreSelfField();
    /**
     * @brief The derivative k_i from P, E and B, and the next stage state y0 + c_{i+1} dt k_i,
     *        or after the last stage y0 + dt sum_j b_j k_j, with the particle's own dt.
     * @throws OpalException if the local particle count of an active container differs from
     *         the one holdSelfField() saw.
     */
    void advanceStage(unsigned i);

private:
    /**
     * @brief Per-particle data of one container for one step: the start state, the sum of
     *        b_i k_i and the held self field, indexed like the particles.
     *
     * Plain views, not particle attributes: they are filled after the self-field solve and the
     * emission, the only operations of the step that create, delete, migrate or reorder
     * particles, and are read until the last stage. Between those, OpalBeamline::
     * applyExternalFields() changes frames in place and the element apply() overrides only add
     * fields or record (Monitor); losses are marked after the step. advanceStage() checks the
     * particle count.
     */
    struct Scratch {
        Kokkos::View<Vector_t<double, 3>*> startR, startP, sumR, sumP, heldE, heldB;
        std::size_t count = 0;  ///< Local particles held; 0 for an inactive container.
    };

    std::vector<Scratch> scratch_m;  ///< One per container; the views only grow.
    IpplTimings::TimerRef stageTimer_m;
};

#endif  // OPALX_RUNGE_KUTTA_INTEGRATOR_H
