//
// Class NumericalIntegrator
//   The particle part of one time step, as an exchangeable integration scheme.
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
#include "Utilities/OpalException.h"

#include <utility>

namespace {

    class BorisMidpoint final : public NumericalIntegrator {
    public:
        explicit BorisMidpoint(Context context) : NumericalIntegrator(std::move(context)) {}

        void advance() override {
            halfDrift();
            solveSelfFields();
            emit();
            addExternalFields();
            kick();
            halfDrift();
        }
    };

    class BorisPrestep final : public NumericalIntegrator {
    public:
        explicit BorisPrestep(Context context) : NumericalIntegrator(std::move(context)) {}

        void advance() override {
            solveSelfFields();
            halfDrift();
            emit();
            addExternalFields();
            kick();
            halfDrift();
        }
    };

}  // namespace

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

void NumericalIntegrator::emit() {
    IpplTimings::startTimer(emissionTimer_m);
    context_m.emittingSamplers.emit(
            context_m.bunch, context_m.beamline, context_m.bunch.getT(), context_m.bunch.getdT());
    IpplTimings::stopTimer(emissionTimer_m);
}

void NumericalIntegrator::addExternalFields() {
    IpplTimings::startTimer(externalFieldsTimer_m);
    // The source plane is the X/Y plane at z = R0Z of the emitting EMISSIONSOURCE. Space-charge
    // kicks can push particles behind it; the loss cannot wait for the cleanup after the step.
    const size_t nSourceMarked = context_m.beamline.markBackwardParticlesAtSourcePlane(
            context_m.bunch, context_m.dirichletPlane);
    if (nSourceMarked > 0) {
        Inform m("NumericalIntegrator::addExternalFields");
        context_m.bunch.deleteInvalidParticles(
                true, m, "backward source-plane particles during external-field evaluation");
    }

    context_m.beamline.applyExternalFields(context_m.bunch);
    IpplTimings::stopTimer(externalFieldsTimer_m);
}

std::unique_ptr<NumericalIntegrator> makeNumericalIntegrator(
        const std::string& name, NumericalIntegrator::Context context) {
    if (name == "BORIS") return std::make_unique<BorisMidpoint>(std::move(context));
    if (name == "BORIS-PRESTEP") return std::make_unique<BorisPrestep>(std::move(context));
    throw OpalException("makeNumericalIntegrator", "Unknown time integrator \"" + name + "\".");
}
