/**
 * @file ParallelTracker.cpp
 * @brief Implementation of ParallelTracker (OPAL-T time-based parallel tracker).
 *
 * @copyright Copyright (c) 200x - 2014, Christof Kraus, Paul Scherrer Institut, Villigen PSI,
 * Switzerland
 * @copyright 2015 - 2016, Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
 * @copyright 2017 - 2020, Christof Metzger-Kraus
 *
 * All rights reserved. This file is part of OPAL.
 *
 * OPAL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * You should have received a copy of the GNU General Public License
 * along with OPAL. If not, see https://www.gnu.org/licenses/.
 */
#include "Algorithms/ParallelTracker.h"

#include <cmath>
#include <string>
#include <vector>

#include "Structure/EMFieldDumpWriter.h"

#include "SpaceCharge/SpaceChargeSolver.h"

#include "Structure/CheckpointFile.h"
#include "Structure/DataSink.h"
#include "Utilities/Timer.h"
#include "Utilities/Util.h"

extern Inform* gmsg;

// --- Constructors ---

/**
 * @brief Construct the tracker from what buildTracking() built for one RUN.
 */
ParallelTracker::ParallelTracker(TrackingSetup setup)
    : itsBunch_m(setup.bunch),
      options_m(std::move(setup.options)),
      spaceChargeFieldUpdate_m(setup.spaceChargeFieldUpdate),
      itsDataSink_m(setup.dataSink),
      spaceChargeSolver_m(setup.spaceChargeSolver),
      dirichletPlane_m(std::move(setup.dirichletPlane)),
      itsOpalBeamline_m(std::move(setup.beamline)),
      beamStops_m(std::move(setup.beamStops)),
      emittingSamplers_m(std::move(setup.emittingSamplers)),
      dtCurrentTrack_m(0.0),
      sStart_m(setup.sStart),
      globalPhaseShift_m(setup.globalPhaseShift),
      restarting_m(setup.restarting),
      restartGlobalStep_m(setup.restartGlobalStep),
      restartDt_m(setup.restartDt),
      restartPosition_m(setup.restartPosition),
      checkpointPath_m(std::move(setup.checkpointPath)),
      emFieldDumps_m(std::move(setup.emFieldDumps)),
      trackStartTimer_m(IpplTimings::getTimer("Track: start")),
      stepTimeStepTimer_m(IpplTimings::getTimer("Step: time step")),
      stepSpaceChargeTimer_m(IpplTimings::getTimer("Step: space charge")),
      stepFirstPushTimer_m(IpplTimings::getTimer("Step: first push")),
      stepEmissionTimer_m(IpplTimings::getTimer("Step: emission")),
      stepExternalFieldsTimer_m(IpplTimings::getTimer("Step: ext. fields")),
      stepSpinTimer_m(IpplTimings::getTimer("Step: spin")),
      stepKickAndPushTimer_m(IpplTimings::getTimer("Step: kick and push")),
      stepProcessesTimer_m(IpplTimings::getTimer("Step: processes")),
      stepReferenceTimer_m(IpplTimings::getTimer("Step: reference")),
      stepRemoveLostTimer_m(IpplTimings::getTimer("Step: remove lost")),
      stepOutputTimer_m(IpplTimings::getTimer("Step: output")),
      stepStopRulesTimer_m(IpplTimings::getTimer("Step: stop rules")),
      stepCheckpointTimer_m(IpplTimings::getTimer("Step: checkpoint")),
      trackFinishTimer_m(IpplTimings::getTimer("Track: finish")),
      references_m(
              *setup.bunch, itsOpalBeamline_m, std::move(setup.beamReferenceMomenta),
              options_m.autoPhase, std::move(setup.cavityPhases)) {
    for (unsigned int i = 0; i < setup.sStop.size(); ++i) {
        stepSizes_m.push_back(setup.dt[i], setup.sStop[i], setup.maxSteps[i]);
    }

    stepSizes_m.sortAscendingSStop();
    stepSizes_m.resetIterator();
}

// --- execute() ---

void ParallelTracker::execute() {
    IpplTimings::startTimer(trackStartTimer_m);
    const bool started = startTracking();
    IpplTimings::stopTimer(trackStartTimer_m);
    if (started) trackSteps();
}

bool ParallelTracker::startTracking() {
    Inform m("ParallelTracker::startTracking");
    // The reference start and pass may cross a MONITOR; only the step loop records.
    itsOpalBeamline_m.setMonitorsRecordReference(false);
    if (restarting_m) {
        stepSizes_m.resumeFromCheckpoint(restartPosition_m, restartDt_m, restartGlobalStep_m);
    }

    // PartBunch::resetPcActive() ran in the constructor while containers were still empty
    // (allocate-then-destroy for capacity). Initial particles are loaded later in TrackRun
    // without refreshing these flags, so reference updates must not skip all containers.
    if (!restarting_m || restartPosition_m.stepsCompletedInSegment == 0) {
        itsBunch_m->resetPcActive();
    }
    emittingSamplers_m.activateEmittingContainers(*itsBunch_m, itsBunch_m->getT());

    // Ensure the time step is positive for the setup phase (buildTracking() did so as well,
    // before the elements were initialised).
    itsBunch_m->setdT(std::abs(itsBunch_m->getdT()));

    double minTimeStep = stepSizes_m.getMinTimeStep();
    m << level3 << "Selected minimum time step from configuration: " << minTimeStep << endl;

    placeReferences();

    if (stepSizes_m.reachedEnd()) {
        *gmsg << level1 << "* Checkpoint is already at the end of the configured tracking "
              << "schedule; no integration steps remain." << endl;
        itsOpalBeamline_m.switchElementsOff();
        Kokkos::fence();
        return false;
    }

    references_m.tracePrimary(stepSizes_m, minTimeStep);

    // Legacy OPAL emission starts a fresh tracker before the RF reference time for centered
    // flat-top pulses. Checkpoint time already contains this shift and must not be shifted again.
    double time = itsBunch_m->getT();
    if (!restarting_m) {
        time -= globalPhaseShift_m;
        itsBunch_m->setT(time);
    }
    m << level4 << "Reset bunch time to " << time << "." << endl;

    printTrackStart(time);

    // Handle any dump field requests
    writeEMFieldDumps(emFieldDumps_m, itsOpalBeamline_m.getElements());

    // Only crossings in the step loop write MONITOR output.
    itsOpalBeamline_m.setMonitorsRecordReference(true);
    return true;
}

void ParallelTracker::placeReferences() {
    // A fresh run derives each reference pose from its initial distribution. A restart must keep
    // the per-container reference pose and reference-to-lab transform stored in the checkpoint.
    if (!restarting_m) {
        references_m.setInitialPoses(emittingSamplers_m);

        // Integrate reference orbits forward until all container path lengths reach sStart_m.
        const auto selectTimeStep = [this](double dt) {
            dtCurrentTrack_m = dt;
            itsBunch_m->setdT(selectDT());
        };
        if (references_m.findStartPositions(sStart_m, stepSizes_m, selectTimeStep)) changeDT();
        stepSizes_m.advanceToPos(sStart_m);
    } else {
        stepSizes_m.advanceToIndex(restartPosition_m.segment);
    }
    // Decay daughter containers use their root's reference particle, frame and path length.
    // Copy them after the start search, and after a restart in case the checkpoint differs.
    references_m.copyFollowers();
}

void ParallelTracker::printTrackStart(double time) {
    const unsigned long long step = itsBunch_m->getGlobalTrackStep();
    OPALTimer::Timer myt1;
    *gmsg << level1 << "* Track start at: " << myt1.time() << ", t= " << Util::getTimeString(time)
          << "; "
          << "sStart at: " << Util::getLengthString(itsBunch_m->getParticleContainer(0)->get_sPos())
          << endl
          << "* Initial dt = " << Util::getTimeString(itsBunch_m->getdT()) << endl
          << "* Max integration steps = " << stepSizes_m.getMaxSteps() << ", next step = " << step
          << endl;

    beamStops_m.announce(itsBunch_m->getNumParticleContainers());

    stepSizes_m.printDirect(*gmsg);
}

void ParallelTracker::trackSteps() {
    Inform m("ParallelTracker::trackSteps");

    // The position within the step-size schedule.
    std::size_t segmentIndex = stepSizes_m.getCurrentIndex();
    unsigned long long stepsCompletedInSegment =
            restarting_m ? restartPosition_m.stepsCompletedInSegment : 0;

    bool stopRun = false;
    while (!stepSizes_m.reachedEnd()) {
        IpplTimings::startTimer(stepTimeStepTimer_m);
        const unsigned long long stepsInSegment = stepSizes_m.getNumSteps();
        startSegment(stepsCompletedInSegment);

        while (stepsCompletedInSegment < stepsInSegment) {
            // In the first pass the timer still runs from startSegment(); starting a running
            // timer does nothing.
            IpplTimings::startTimer(stepTimeStepTimer_m);
            // Resolve emission-dependent dt before the step so the reference, particle
            // transport and self-field evaluation share the same duration.
            changeDT();
            if (!itsBunch_m->anyPcActive()) {
                m << level4 << "No active particle containers; ending inner track segment." << endl;
                IpplTimings::stopTimer(stepTimeStepTimer_m);
                break;
            }
            IpplTimings::stopTimer(stepTimeStepTimer_m);

            advanceParticles();

            IpplTimings::startTimer(stepReferenceTimer_m);
            advanceReference();
            IpplTimings::stopTimer(stepReferenceTimer_m);

            IpplTimings::startTimer(stepRemoveLostTimer_m);
            const auto due         = itsDataSink_m->outputsDue(itsBunch_m->getGlobalTrackStep());
            const auto stepMoments = removeLostParticles(due);
            IpplTimings::stopTimer(stepRemoveLostTimer_m);

            IpplTimings::startTimer(stepOutputTimer_m);
            itsDataSink_m->writeStep(
                    *itsBunch_m, due, stepMoments, referenceFields(),
                    spaceChargeSolver_m->reportedBinCount());
            IpplTimings::stopTimer(stepOutputTimer_m);

            IpplTimings::startTimer(stepStopRulesTimer_m);
            itsBunch_m->incTrackSteps();
            ++stepsCompletedInSegment;

            const auto stop =
                    beamStops_m.applyStopRules(*itsBunch_m, stepSizes_m, references_m.hitElement());
            const bool segmentComplete = !itsBunch_m->anyPcActive()
                                         || stepsCompletedInSegment == stepsInSegment
                                         || stop.switchSegment;
            IpplTimings::stopTimer(stepStopRulesTimer_m);
            writeCheckpointIfDue(segmentComplete, segmentIndex, stepsCompletedInSegment);

            if (stop.runFinished) {
                stopRun = true;
                break;
            }
            if (stop.switchSegment) break;
            if (!itsBunch_m->anyPcActive()) {
                m << level2
                  << "All active containers reached current sStop. Preparing to switch to "
                  << "next configuration." << endl;
                break;
            }
        }
        // Still running if the inner loop body never ran.
        IpplTimings::stopTimer(stepTimeStepTimer_m);

        if (stopRun) break;
        ++stepSizes_m;
        ++segmentIndex;
        stepsCompletedInSegment = 0;
    }
    IpplTimings::startTimer(trackFinishTimer_m);
    finishTracking();
    IpplTimings::stopTimer(trackFinishTimer_m);
}

void ParallelTracker::startSegment(unsigned long long stepsCompletedInSegment) {
    Inform m("ParallelTracker::startSegment");
    const unsigned long long step = itsBunch_m->getGlobalTrackStep();
    const unsigned long long trackSteps =
            step + (stepSizes_m.getNumSteps() - stepsCompletedInSegment);
    dtCurrentTrack_m = stepSizes_m.getdT();

    // Select global dt from dtCurrentTrack_m and copy to all container dt views.
    changeDT();
    if (stepsCompletedInSegment == 0) {
        itsBunch_m->resetPcActive();
    }
    emittingSamplers_m.activateEmittingContainers(*itsBunch_m, itsBunch_m->getT());
    // Beams that passed their ZSTOP plane stay stopped.
    beamStops_m.keepDoneBeamsStopped(*itsBunch_m);

    m << level2 << "Starting track with dt = " << Util::getTimeString(dtCurrentTrack_m)
      << ", track steps = " << step << " to " << trackSteps << "." << endl;
}

void ParallelTracker::advanceParticles() {
    Inform m("ParallelTracker::advanceParticles");
    const long long step = itsBunch_m->getGlobalTrackStep();

    // Historical OPAL evaluated the self-field at R_n and carried the gathered per-particle field
    // through the first half drift. Keep this as an explicit compatibility mode; the default
    // evaluates the field at R_{n+1/2} below.
    if (spaceChargeFieldUpdate_m == SpaceChargeFieldUpdate::PRESTEP) {
        IpplTimings::startTimer(stepSpaceChargeTimer_m);
        itsBunch_m->resetFields();
        computeSpaceChargeFields();
        IpplTimings::stopTimer(stepSpaceChargeTimer_m);
    }

    // First half of the time integration
    IpplTimings::startTimer(stepFirstPushTimer_m);
    timeIntegration1();
    m << level4 << "timeIntegration1 done at step " << step << "." << endl;
    IpplTimings::stopTimer(stepFirstPushTimer_m);

    if (spaceChargeFieldUpdate_m == SpaceChargeFieldUpdate::MIDPOINT) {
        IpplTimings::startTimer(stepSpaceChargeTimer_m);
        itsBunch_m->resetFields();
        computeSpaceChargeFields();
        IpplTimings::stopTimer(stepSpaceChargeTimer_m);
    }

    // Emit particles in the reference frame with their own fractional timesteps.
    // They experience external fields in their creation step, but no space charge.
    IpplTimings::startTimer(stepEmissionTimer_m);
    emittingSamplers_m.emit(
            *itsBunch_m, itsOpalBeamline_m, itsBunch_m->getT(), itsBunch_m->getdT());
    IpplTimings::stopTimer(stepEmissionTimer_m);

    IpplTimings::startTimer(stepExternalFieldsTimer_m);
    computeExternalFields();
    IpplTimings::stopTimer(stepExternalFieldsTimer_m);

    // First half of the Thomas-BMT spin rotation (momentum before the kick);
    // timeIntegration2 does the second half after the kick.
    // No-op for containers that have not registered the Pol attribute.
    IpplTimings::startTimer(stepSpinTimer_m);
    itsBunch_m->evolveSpinHalfStep();
    IpplTimings::stopTimer(stepSpinTimer_m);

    // Second half of the time integration
    IpplTimings::startTimer(stepKickAndPushTimer_m);
    timeIntegration2();
    IpplTimings::stopTimer(stepKickAndPushTimer_m);

    // Apply global processes (e.g. decay) and mark affected particles for deletion
    IpplTimings::startTimer(stepProcessesTimer_m);
    const size_t nProcessMarked = itsBunch_m->applyGlobalProcesses(itsBunch_m->getdT());
    if (nProcessMarked > 0) {
        // Decay-like global processes may have created daughters and marked parents
        // invalid. Delete the parents before time/reference updates so they do not
        // contribute to post-step moments.
        itsBunch_m->deleteInvalidParticles(false, m, "global processes");
    }

    // Small sanity check
    if (itsBunch_m->getTotalNumAllContainers() == 0) {
        m << level5 << "WARNING: No particles in the bunch at step " << step << " on rank "
          << ippl::Comm->rank() << ". This has no effect on the simulation." << endl;
    }
    IpplTimings::stopTimer(stepProcessesTimer_m);
}

void ParallelTracker::advanceReference() {
    Inform m("ParallelTracker::advanceReference");
    itsBunch_m->incrementT();
    m << level5 << "Incremented bunch time to " << Util::getTimeString(itsBunch_m->getT()) << "."
      << endl;

    references_m.advance();
}

std::vector<opalx::moments::BeamMoments> ParallelTracker::removeLostParticles(
        const DataSink::StepOutputs& due) {
    Inform m("ParallelTracker::removeLostParticles");

    // Every frame change of this step is closed: moments, output and the next step read
    // reference-frame data.
    itsBunch_m->requireReferenceFrame("ParallelTracker::removeLostParticles");

    // Enforce the transverse aperture of every nearby element (collective call).
    size_t nStepMarked = itsOpalBeamline_m.applyElementApertures(*itsBunch_m);
    if (nStepMarked > 0) {
        m << level2 << "Marked " << nStepMarked
          << " particles outside a transverse aperture for deletion." << endl;
    }

    // The outlier box and the outputs of this step read one computation per container,
    // unless particles are deleted in between.
    const double sigmas = options_m.boundpDestroy;
    const auto outputs  = [&](const ParticleContainer_t& pc) {
        return itsDataSink_m->requiredMoments(due, pc.hasSpin());
    };
    auto stepMoments = itsBunch_m->markOutliers(sigmas, outputs, nStepMarked);

    // Single flush of all marks accumulated this step.
    if (nStepMarked > 0) {
        itsBunch_m->deleteInvalidParticles(
                true, m, "transverse aperture / " + std::to_string(sigmas) + "-sigma boundary");
        stepMoments.assign(stepMoments.size(), {});
    }
    return stepMoments;
}

void ParallelTracker::writeCheckpointIfDue(
        bool segmentComplete, std::size_t segmentIndex,
        unsigned long long stepsCompletedInSegment) {
    if (options_m.checkpointFreq > 0
        && itsBunch_m->getGlobalTrackStep() % options_m.checkpointFreq == 0) {
        IpplTimings::startTimer(stepCheckpointTimer_m);
        Inform m("ParallelTracker::writeCheckpointIfDue");
        const std::size_t nextSegment = segmentComplete ? segmentIndex + 1 : segmentIndex;
        const unsigned long long nextSegmentSteps = segmentComplete ? 0 : stepsCompletedInSegment;
        CheckpointFile::write(
                checkpointPath_m, *itsBunch_m, nextSegment, nextSegmentSteps, beamStops_m.stops(),
                references_m.cavityPhases());
        m << level2 << "Wrote checkpoint '" << checkpointPath_m << "' after global step "
          << itsBunch_m->getGlobalTrackStep() << "." << endl;
        IpplTimings::stopTimer(stepCheckpointTimer_m);
    }
}

void ParallelTracker::finishTracking() {
    Inform m("ParallelTracker::finishTracking");
    itsDataSink_m->writeLastStep(
            *itsBunch_m, referenceFields(), spaceChargeSolver_m->reportedBinCount());

    itsOpalBeamline_m.switchElementsOff();

    // Ensure all Kokkos operations are complete
    Kokkos::fence();

    m << level2 << "Total FieldSolver calls: " << spaceChargeSolver_m->backendSolveCount() << endl;
    if (ippl::Comm->size() > 1) {
        m << level2 << "Total binary repartitions: " << spaceChargeSolver_m->redistributionCount()
          << endl;
    }

    OPALTimer::Timer myt3;
    *gmsg << level1 << endl
          << "* Done executing ParallelTracker at " << myt3.time() << endl
          << endl;
}

// --- PIC integration and fields ---

void ParallelTracker::timeIntegration1() {
    Inform m("ParallelTracker::timeIntegration1");
    itsBunch_m->pushParticles();
    m << level4 << "Push particles done for all containers." << endl;
}

void ParallelTracker::timeIntegration2() {
    Inform m("ParallelTracker::timeIntegration2");

    itsBunch_m->kickAndPushParticles();
    m << level4 << "Kick/push particles done for all containers." << endl;
}

void ParallelTracker::computeSpaceChargeFields() {
    spaceChargeSolver_m->solveStep(*itsBunch_m, options_m.minBinEmitted, [this] {
        return emittingSamplers_m.spaceChargeEmissionProgress(
                itsBunch_m->getT(), options_m.aggressiveStateSync);
    });
}

void ParallelTracker::computeExternalFields() {
    Inform msg("ParallelTracker ", *gmsg);

    // Source-plane loss has to be handled here, before external fields are evaluated. Deferring it
    // to the post-step invalid-particle cleanup would let particles behind the cathode receive one
    // more field kick.
    //
    // "Source-plane" refers to the X/Y plane at z=R0Z defined through the emitting EMISSIONSOURCE,
    // from where particles are emitted. It can happen that particles land behind the source plane
    // because e.g. of space-charge kicks. Particles with negative z momentum located behind R0Z
    // (with a certain threshold) are then removed from the simulation.
    const size_t nSourceMarked =
            itsOpalBeamline_m.markBackwardParticlesAtSourcePlane(*itsBunch_m, dirichletPlane_m);
    if (nSourceMarked > 0) {
        itsBunch_m->deleteInvalidParticles(
                true, msg, "backward source-plane particles during external-field evaluation");
    }

    itsOpalBeamline_m.applyExternalFields(*itsBunch_m);
}

// --- Time step ---

double ParallelTracker::selectDT() const {
    return emittingSamplers_m.emissionTimeStep(itsBunch_m->getT()).value_or(dtCurrentTrack_m);
}

void ParallelTracker::changeDT() {
    Inform m("ParallelTracker::changeDT");
    itsBunch_m->setTimeStep(selectDT());
    m << level5 << "Changed particle container time step to " << itsBunch_m->getdT() << "." << endl;
}

std::function<DataSink::ReferenceFields()> ParallelTracker::referenceFields() {
    return [this] {
        return itsOpalBeamline_m.getReferenceFields(*itsBunch_m);
    };
}
