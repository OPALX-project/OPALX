//
// Class NumericalIntegrator
//   The operations of a step and the factory makeNumericalIntegrator(); the schemes are in
//   their own files, e.g. BorisIntegrators.h.
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
#include "Steppers/NumericalIntegrator.h"

#include "Distribution/EmittingSamplers.h"
#include "Elements/OpalBeamline.h"
#include "PartBunch/PartBunch.h"
#include "SpaceCharge/SpaceChargeSolver.h"
#include "Steppers/BorisIntegrators.h"
#include "Steppers/RungeKuttaIntegrator.h"
#include "Utilities/OpalException.h"

#include <utility>

NumericalIntegrator::NumericalIntegrator(Context context)
    : context_m(std::move(context)),
      pushTimer_m(IpplTimings::getTimer("Step: push")),
      spaceChargeTimer_m(IpplTimings::getTimer("Step: space charge")),
      emissionTimer_m(IpplTimings::getTimer("Step: emission")),
      externalFieldsTimer_m(IpplTimings::getTimer("Step: ext. fields")),
      spinTimer_m(IpplTimings::getTimer("Step: spin")),
      kickTimer_m(IpplTimings::getTimer("Step: kick")) {}

void NumericalIntegrator::halfDrift() {
    IpplTimings::startTimer(pushTimer_m);
    context_m.bunch.pushParticles();
    IpplTimings::stopTimer(pushTimer_m);
}

void NumericalIntegrator::kick() {
    // The two spin half rotations see the momentum before and after the kick.
    IpplTimings::startTimer(spinTimer_m);
    context_m.bunch.evolveSpinHalfStep();
    IpplTimings::stopTimer(spinTimer_m);

    IpplTimings::startTimer(kickTimer_m);
    context_m.bunch.kickParticles();
    IpplTimings::stopTimer(kickTimer_m);

    IpplTimings::startTimer(spinTimer_m);
    context_m.bunch.evolveSpinHalfStep();
    IpplTimings::stopTimer(spinTimer_m);
}

void NumericalIntegrator::solveSelfFields() {
    IpplTimings::startTimer(spaceChargeTimer_m);
    context_m.bunch.resetFields();
    context_m.spaceChargeSolver.solveStep(context_m.bunch, context_m.minBinEmitted, [this] {
        return context_m.emittingSamplers.spaceChargeEmissionProgress(
                context_m.bunch.getT(), context_m.aggressiveStateSync);
    });
    IpplTimings::stopTimer(spaceChargeTimer_m);
}

void NumericalIntegrator::emit(EmissionStart start) {
    IpplTimings::startTimer(emissionTimer_m);
    context_m.emittingSamplers.emit(
            context_m.bunch, context_m.beamline, context_m.bunch.getT(), context_m.bunch.getdT(),
            start);
    IpplTimings::stopTimer(emissionTimer_m);
}

void NumericalIntegrator::applySourcePlaneLoss() {
    IpplTimings::startTimer(externalFieldsTimer_m);
    const size_t nSourceMarked = context_m.beamline.markBackwardParticlesAtSourcePlane(
            context_m.bunch, context_m.dirichletPlane);
    if (nSourceMarked > 0) {
        Inform m("NumericalIntegrator::applySourcePlaneLoss");
        context_m.bunch.deleteInvalidParticles(true, m, "backward source-plane particles");
    }
    IpplTimings::stopTimer(externalFieldsTimer_m);
}

void NumericalIntegrator::addExternalFields(double fieldFraction, bool record) {
    IpplTimings::startTimer(externalFieldsTimer_m);
    context_m.beamline.applyExternalFields(context_m.bunch, fieldFraction, record);
    IpplTimings::stopTimer(externalFieldsTimer_m);
}

track_reference::State NumericalIntegrator::advanceReference(
        const track_reference::State& start, double dt, double endTime, double mass, double charge,
        const track_reference::FieldEvaluator& field) const {
    return track_reference::advance(start, dt, endTime, mass, charge, field);
}

std::unique_ptr<NumericalIntegrator> makeNumericalIntegrator(
        const std::string& name, NumericalIntegrator::Context context) {
    if (name == "BORIS") return std::make_unique<BorisMidpoint>(std::move(context));
    if (name == "BORIS-PRESTEP") return std::make_unique<BorisPrestep>(std::move(context));
    if (name == "RK4") return std::make_unique<RungeKutta4>(std::move(context));
    throw OpalException("makeNumericalIntegrator", "Unknown time integrator \"" + name + "\".");
}
