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
#include "Steppers/RungeKuttaIntegrator.h"

#include "PartBunch/PartBunch.h"
#include "Steppers/PhaseSpaceRate.h"
#include "Steppers/RungeKuttaTableau.h"
#include "Utilities/OpalException.h"

#include <string>
#include <utility>

namespace {
    using runge_kutta::rk4;
    constexpr unsigned numStages = rk4.c.size();

    /// Whether stage i of @p table is y0 + c_i h k_{i-1}: the scheme keeps per particle only
    /// y0, the last k and the sum of b_j k_j.
    template <unsigned Stages>
    constexpr bool usesOnlyThePreviousStage(const runge_kutta::Tableau<Stages>& table) {
        if (table.c[0] != 0.0) return false;
        for (unsigned i = 1; i < Stages; ++i) {
            for (unsigned j = 0; j < i; ++j) {
                if (table.a[i][j] != (j + 1 == i ? table.c[i] : 0.0)) return false;
            }
        }
        return true;
    }
    static_assert(usesOnlyThePreviousStage(rk4));

    /// The stage at which MONITORs record; Monitor::apply() expects the state half a step after
    /// the start.
    constexpr unsigned recordingStage = 1;
    static_assert(rk4.c[recordingStage] == 0.5);
}  // namespace

RungeKutta4::RungeKutta4(Context context)
    : NumericalIntegrator(std::move(context)),
      stageTimer_m(IpplTimings::getTimer("Step: RK4 stages")) {}

void RungeKutta4::advance() {
    applySourcePlaneLoss();
    solveSelfFields();
    emit(EmissionStart::Birth);
    holdSelfField();
    for (unsigned i = 0; i < numStages; ++i) {
        if (i > 0) restoreSelfField();
        addExternalFields(rk4.c[i], i == recordingStage);
        advanceStage(i);
    }
}

void RungeKutta4::holdSelfField() {
    IpplTimings::startTimer(stageTimer_m);
    scratch_m.resize(bunch().getNumParticleContainers());
    for (std::size_t ci = 0; ci < scratch_m.size(); ++ci) {
        Scratch& scratch = scratch_m[ci];
        const auto pc    = bunch().getParticleContainer(ci);
        scratch.count    = bunch().isPcActive(ci) ? pc->getLocalNum() : 0;
        if (scratch.heldE.extent(0) < scratch.count) {
            for (auto* view :
                 {&scratch.startR, &scratch.startP, &scratch.sumR, &scratch.sumP, &scratch.heldE,
                  &scratch.heldB}) {
                Kokkos::realloc(Kokkos::WithoutInitializing, *view, scratch.count);
            }
        }
        const auto range = std::make_pair(std::size_t{0}, scratch.count);
        Kokkos::deep_copy(
                Kokkos::subview(scratch.heldE, range), Kokkos::subview(pc->E.getView(), range));
        Kokkos::deep_copy(
                Kokkos::subview(scratch.heldB, range), Kokkos::subview(pc->B.getView(), range));
    }
    IpplTimings::stopTimer(stageTimer_m);
}

void RungeKutta4::restoreSelfField() {
    IpplTimings::startTimer(stageTimer_m);
    for (std::size_t ci = 0; ci < scratch_m.size(); ++ci) {
        if (!bunch().isPcActive(ci)) continue;
        const Scratch& scratch = scratch_m[ci];
        const auto pc          = bunch().getParticleContainer(ci);
        const auto range       = std::make_pair(std::size_t{0}, scratch.count);
        Kokkos::deep_copy(
                Kokkos::subview(pc->E.getView(), range), Kokkos::subview(scratch.heldE, range));
        Kokkos::deep_copy(
                Kokkos::subview(pc->B.getView(), range), Kokkos::subview(scratch.heldB, range));
    }
    IpplTimings::stopTimer(stageTimer_m);
}

void RungeKutta4::advanceStage(unsigned i) {
    IpplTimings::startTimer(stageTimer_m);
    const bool first     = i == 0;
    const bool last      = i + 1 == numStages;
    const double weight  = rk4.b[i];
    const double next    = last ? 0.0 : rk4.c[i + 1];
    for (std::size_t ci = 0; ci < scratch_m.size(); ++ci) {
        if (!bunch().isPcActive(ci)) continue;
        const Scratch& scratch = scratch_m[ci];
        const auto pc          = bunch().getParticleContainer(ci);
        if (pc->getLocalNum() != scratch.count) {
            throw OpalException(
                    "RungeKutta4::advanceStage",
                    "The local particle count of container " + std::to_string(ci) + " changed from "
                            + std::to_string(scratch.count) + " to "
                            + std::to_string(pc->getLocalNum()) + " within the step.");
        }
        const PartData& reference = *pc->getReference();
        const double mass         = reference.getM();
        const double charge       = reference.getQ();

        auto R      = pc->R.getView();
        auto P      = pc->P.getView();
        auto E      = pc->E.getView();
        auto B      = pc->B.getView();
        auto dt     = pc->dt.getView();
        auto startR = scratch.startR;
        auto startP = scratch.startP;
        auto sumR   = scratch.sumR;
        auto sumP   = scratch.sumP;
        Kokkos::parallel_for(
                "RungeKutta4::advanceStage", scratch.count, KOKKOS_LAMBDA(const std::size_t j) {
                    if (first) {
                        startR(j) = R(j);
                        startP(j) = P(j);
                    }
                    const PhaseSpaceRate k = phaseSpaceRate(P(j), E(j), B(j), charge, mass);
                    if (first) {
                        sumR(j) = weight * k.dR;
                        sumP(j) = weight * k.dP;
                    } else {
                        sumR(j) += weight * k.dR;
                        sumP(j) += weight * k.dP;
                    }
                    const double h = dt(j);
                    if (last) {
                        R(j) = startR(j) + h * sumR(j);
                        P(j) = startP(j) + h * sumP(j);
                    } else {
                        R(j) = startR(j) + (next * h) * k.dR;
                        P(j) = startP(j) + (next * h) * k.dP;
                    }
                });
        Kokkos::fence();
    }
    IpplTimings::stopTimer(stageTimer_m);
}
