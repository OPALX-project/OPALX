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

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

#include "BasicActions/DumpEMFields.h"

#include "Physics/Physics.h"
#include "Physics/Units.h"

#include "Processes/GlobalProcesses/GlobalProcess.h"

#include "SpaceCharge/SpaceChargeSolver.h"

#include "Structure/CheckpointFile.h"
#include "Structure/DataSink.h"
#include "Utilities/OpalException.h"
#include "Utilities/Timer.h"
#include "Utilities/Util.h"

extern Inform* gmsg;

using opalx::moments::Moment;

namespace {
    /// Whether an output with frequency @p freq (0: never) is due at global step @p step: the
    /// last step of every cycle of @p freq steps.
    bool outputDue(int freq, long long step) { return freq > 0 && (step % freq) + 1 == freq; }
}  // namespace

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
      sStart_m(setup.sStart),
      globalPhaseShift_m(setup.globalPhaseShift),
      dtCurrentTrack_m(0.0),
      emittingSamplers_m(std::move(setup.emittingSamplers)),
      restarting_m(setup.restarting),
      restartGlobalStep_m(setup.restartGlobalStep),
      restartDt_m(setup.restartDt),
      restartPosition_m(setup.restartPosition),
      checkpointPath_m(std::move(setup.checkpointPath)),
      timeIntegrationTimer1_m(IpplTimings::getTimer("TIntegration1")),
      timeIntegrationTimer2_m(IpplTimings::getTimer("TIntegration2")),
      fieldEvaluationTimer_m(IpplTimings::getTimer("External field eval")),
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
    if (startTracking()) trackSteps();
}

bool ParallelTracker::startTracking() {
    Inform m("ParallelTracker::startTracking");
    // The reference start and pass may cross a MONITOR; only the step loop records.
    itsOpalBeamline_m.setMonitorsRecordReference(false);
    if (restarting_m) {
        resumeRestartSchedule();
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

    // A fresh run derives each reference pose from its initial distribution. A restart must keep
    // the per-container reference pose and reference-to-lab transform stored in the checkpoint.
    if (!restarting_m) {
        references_m.setInitialPoses(emittingSamplers_m);

        // Integrate reference orbits forward until all container path lengths reach sStart_m.
        const auto selectTimeStep = [this](double dt) {
            dtCurrentTrack_m = dt;
            selectDT();
        };
        if (references_m.findStartPositions(sStart_m, stepSizes_m, selectTimeStep)) changeDT();
        stepSizes_m.advanceToPos(sStart_m);
    } else {
        stepSizes_m.advanceToIndex(restartPosition_m.segment);
    }
    // Decay daughter containers use their root's reference particle, frame and path length.
    // Copy them after the start search, and after a restart in case the checkpoint differs.
    references_m.copyFollowers();

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

    const unsigned long long step = itsBunch_m->getGlobalTrackStep();
    OPALTimer::Timer myt1;
    *gmsg << level1 << "* Track start at: " << myt1.time() << ", t= " << Util::getTimeString(time)
          << "; "
          << "sStart at: " << Util::getLengthString(itsBunch_m->getParticleContainer(0)->get_sPos())
          << endl
          << "* Initial dt = " << Util::getTimeString(itsBunch_m->getdT()) << endl
          << "* Max integration steps = " << stepSizes_m.getMaxSteps() << ", next step = " << step
          << endl;

    beamStops_m.resize(itsBunch_m->getNumParticleContainers());
    for (size_t i = 0; i < beamStops_m.size(); ++i) {
        const BeamStop& stop = beamStops_m[i];
        if (stop.zStop) {
            *gmsg << level1 << "* Beam " << i << " stops when its last particle passes ZSTOP = "
                  << Util::getLengthString(*stop.zStop) << "." << endl;
        } else if (stop.untilMaxSteps) {
            *gmsg << level1 << "* Beam " << i << " is tracked until MAXSTEPS." << endl;
        }
    }

    stepSizes_m.printDirect(*gmsg);

    // Handle any dump field requests
    DumpEMFields::writeFields(itsOpalBeamline_m.getElements());

    // Only crossings in the step loop write MONITOR output.
    itsOpalBeamline_m.setMonitorsRecordReference(true);
    return true;
}

void ParallelTracker::resumeRestartSchedule() {
    Inform m("ParallelTracker::resumeRestartSchedule");
    stepSizes_m.advanceToResumePosition(restartPosition_m);
    if (!stepSizes_m.reachedEnd()) {
        const double configuredDt = stepSizes_m.getdT();
        const double tolerance    = 100.0 * std::numeric_limits<double>::epsilon()
                                 * std::max(
                                         {std::abs(configuredDt), std::abs(restartDt_m),
                                          std::numeric_limits<double>::min()});
        if (std::abs(configuredDt - restartDt_m) > tolerance) {
            throw OpalException(
                    "ParallelTracker::resumeRestartSchedule",
                    "restart input changes the time-step schedule at the checkpointed global "
                    "step");
        }
    }
    m << level2 << "Restart global step " << restartGlobalStep_m << " resumes step-size "
      << "segment " << restartPosition_m.segment << " with "
      << restartPosition_m.stepsCompletedInSegment << " completed steps in that segment." << endl;
}

void ParallelTracker::trackSteps() {
    Inform m("ParallelTracker::trackSteps");

    // The position within the step-size schedule.
    std::size_t segmentIndex = stepSizes_m.getCurrentIndex();
    unsigned long long stepsCompletedInSegment =
            restarting_m ? restartPosition_m.stepsCompletedInSegment : 0;

    bool stopRun = false;
    while (!stepSizes_m.reachedEnd()) {
        const unsigned long long stepsInSegment = stepSizes_m.getNumSteps();
        startSegment(stepsCompletedInSegment);

        while (stepsCompletedInSegment < stepsInSegment) {
            // Resolve emission-dependent dt before the step so the reference, particle
            // transport and self-field evaluation share the same duration.
            changeDT();
            if (!itsBunch_m->anyPcActive()) {
                m << level4 << "No active particle containers; ending inner track segment." << endl;
                break;
            }

            advanceParticles();
            advanceReference();

            const long long step   = itsBunch_m->getGlobalTrackStep();
            const bool psDump      = outputDue(options_m.psDumpFreq, step);
            const bool statDump    = outputDue(options_m.statDumpFreq, step);
            const auto stepMoments = removeLostParticles(psDump, statDump);
            dumpStats(psDump, statDump, stepMoments);

            itsBunch_m->incTrackSteps();
            ++stepsCompletedInSegment;

            const StopRules stop       = applyStopRules();
            const bool segmentComplete = !itsBunch_m->anyPcActive()
                                         || stepsCompletedInSegment == stepsInSegment
                                         || stop.switchSegment;
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

        if (stopRun) break;
        ++stepSizes_m;
        ++segmentIndex;
        stepsCompletedInSegment = 0;
    }
    finishTracking();
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
    for (size_t i = 0; i < beamStops_m.size(); ++i) {
        if (beamStops_m[i].done) itsBunch_m->setPcAtSStop(i);
    }

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
        itsBunch_m->resetFields();
        computeSpaceChargeFields();
    }

    // First half of the time integration
    timeIntegration1();
    m << level4 << "timeIntegration1 done at step " << step << "." << endl;

    if (spaceChargeFieldUpdate_m == SpaceChargeFieldUpdate::MIDPOINT) {
        itsBunch_m->resetFields();
        computeSpaceChargeFields();
    }

    // Emit particles in the reference frame with their own fractional timesteps.
    // They experience external fields in their creation step, but no space charge.
    emittingSamplers_m.emit(
            *itsBunch_m, itsOpalBeamline_m, itsBunch_m->getT(), itsBunch_m->getdT());

    computeExternalFields();

    // First half of the Thomas-BMT spin rotation (momentum before the kick);
    // timeIntegration2 does the second half after the kick.
    // No-op for containers that have not registered the Pol attribute.
    itsBunch_m->evolveSpinHalfStep();

    // Second half of the time integration
    timeIntegration2();

    // Apply global processes (e.g. decay) and mark affected particles for deletion
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
}

void ParallelTracker::advanceReference() {
    Inform m("ParallelTracker::advanceReference");
    itsBunch_m->incrementT();
    m << level5 << "Incremented bunch time to " << Util::getTimeString(itsBunch_m->getT()) << "."
      << endl;

    references_m.advance();
}

std::vector<opalx::moments::BeamMoments> ParallelTracker::removeLostParticles(
        bool psDump, bool statDump) {
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
        return outputMomentsDue(pc, psDump, statDump, stepLogDue());
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

ParallelTracker::StopRules ParallelTracker::applyStopRules() {
    Inform m("ParallelTracker::applyStopRules");
    const auto& particleContainers = itsBunch_m->getParticleContainers();

    // The primary reference switches to the next DT segment at each ZSTOP entry before the
    // last, and a beam with a ZSTOP stops once its last particle has passed its stop plane.
    StopRules stop;
    for (size_t i = 0; i < particleContainers.size(); ++i) {
        const auto& pc = particleContainers[i];
        if (!pc || !itsBunch_m->isPcActive(i)) {
            continue;
        }
        ippl::Vector<double, 3> pdivg = pc->getRefPartP() / Util::getGamma(pc->getRefPartP());
        double beta                   = euclidean_norm(pdivg);
        double driftPerTimeStep       = std::abs(itsBunch_m->getdT()) * Physics::c * beta;
        m << level4 << "Calculated drift per time step (container " << i
          << "): " << Util::getLengthString(driftPerTimeStep) << "." << endl;

        if (i == 0 && !stepSizes_m.isLastSegment()
            && beam_stop::reached(pc->get_sPos(), stepSizes_m.getSStop(), driftPerTimeStep)) {
            m << level2 << "Primary reference reached sStop = "
              << Util::getLengthString(stepSizes_m.getSStop())
              << "; switching to the next DT segment." << endl;
            stop.switchSegment = true;
        }
        if (i < beamStops_m.size() && !itsBunch_m->isFollower(i)) {
            BeamStop& beamStop = beamStops_m[i];
            if (beamStop.zStop && !beamStop.planeSet
                && beam_stop::reached(pc->get_sPos(), *beamStop.zStop, driftPerTimeStep)) {
                beam_stop::setPlane(beamStop, pc->getRefPartR(), pc->getRefPartP(), pc->get_sPos());
                *gmsg << level1 << "* Beam " << i
                      << ": reference reached ZSTOP = " << Util::getLengthString(*beamStop.zStop)
                      << " at t = " << Util::getTimeString(itsBunch_m->getT())
                      << "; tracking until its last particle passes." << endl;
            }
        }
    }
    updateBeamStops();

    bool referenceHit = references_m.hitElement();
    ippl::Comm->allreduce(referenceHit, 1, std::logical_or<bool>());
    if (referenceHit) {
        *gmsg << level1 << "* A reference particle hit an element; the run ends." << endl;
    }
    const bool beamsDone = beam_stop::runFinished(beamStops_m);
    if (beamsDone) {
        *gmsg << level1 << "* Every beam with a ZSTOP is done; the run ends." << endl;
    }
    stop.runFinished = referenceHit || beamsDone;
    return stop;
}

void ParallelTracker::writeCheckpointIfDue(
        bool segmentComplete, std::size_t segmentIndex,
        unsigned long long stepsCompletedInSegment) {
    if (options_m.checkpointFreq > 0
        && itsBunch_m->getGlobalTrackStep() % options_m.checkpointFreq == 0) {
        Inform m("ParallelTracker::writeCheckpointIfDue");
        const std::size_t nextSegment = segmentComplete ? segmentIndex + 1 : segmentIndex;
        const unsigned long long nextSegmentSteps = segmentComplete ? 0 : stepsCompletedInSegment;
        CheckpointFile::write(
                checkpointPath_m, *itsBunch_m, nextSegment, nextSegmentSteps, beamStops_m,
                references_m.cavityPhases());
        m << level2 << "Wrote checkpoint '" << checkpointPath_m << "' after global step "
          << itsBunch_m->getGlobalTrackStep() << "." << endl;
    }
}

void ParallelTracker::finishTracking() {
    Inform m("ParallelTracker::finishTracking");
    // The step loop wrote the last step already if it was due there.
    const long long lastStep = itsBunch_m->getGlobalTrackStep() - 1;
    bool const psDump = options_m.psDumpFreq > 0 && !outputDue(options_m.psDumpFreq, lastStep);
    bool const statDump =
            options_m.statDumpFreq > 0 && !outputDue(options_m.statDumpFreq, lastStep);

    writePhaseSpace(psDump, statDump);

    if (psDump) {
        *gmsg << level2 << "* Dump phase space of last step" << endl;
    }

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
    IpplTimings::startTimer(timeIntegrationTimer1_m);
    itsBunch_m->pushParticles();
    IpplTimings::stopTimer(timeIntegrationTimer1_m);
    m << level4 << "Push particles done for all containers." << endl;
}

void ParallelTracker::timeIntegration2() {
    Inform m("ParallelTracker::timeIntegration2");

    IpplTimings::startTimer(timeIntegrationTimer2_m);
    itsBunch_m->kickAndPushParticles();
    m << level4 << "Kick/push particles done for all containers." << endl;

    IpplTimings::stopTimer(timeIntegrationTimer2_m);
}

void ParallelTracker::computeSpaceChargeFields() {
    if (spaceChargeSolver_m == nullptr) {
        throw OpalException(
                "ParallelTracker::computeSpaceChargeFields",
                "No space-charge solver is available. Use TYPE=NONE for a configured no-op "
                "solver.");
    }
    spaceChargeSolver_m->solveStep(*itsBunch_m, options_m.minBinEmitted, [this] {
        return emittingSamplers_m.spaceChargeEmissionProgress(
                itsBunch_m->getT(), options_m.aggressiveStateSync);
    });
}

void ParallelTracker::computeExternalFields() {
    IpplTimings::startTimer(fieldEvaluationTimer_m);
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

    IpplTimings::stopTimer(fieldEvaluationTimer_m);
}

void ParallelTracker::updateBeamStops() {
    for (size_t i = 0; i < beamStops_m.size(); ++i) {
        BeamStop& stop = beamStops_m[i];
        if (!stop.planeSet || stop.done) continue;

        // Particles are in the bunch frame: move the lab plane into it.
        double minimum = std::numeric_limits<double>::max();
        if (const auto pc = itsBunch_m->getParticleContainer(i)) {
            const CoordinateSystemTrafo& toLab = pc->getToLabTrafo();
            minimum                            = beam_stop::minDistanceToPlane(
                    pc->R.getView(), pc->getLocalNum(), toLab.transformFrom(stop.planeOrigin),
                    toLab.rotateFrom(stop.planeNormal));
        }
        ippl::Comm->allreduce(minimum, 1, std::less<double>());
        if (minimum >= 0.0) {
            stop.done = true;
            itsBunch_m->setPcAtSStop(i);
            *gmsg << level1 << "* Beam " << i << ": its last particle passed ZSTOP at t = "
                  << Util::getTimeString(itsBunch_m->getT()) << "; it is no longer tracked."
                  << endl;
        }
    }
}

// --- Helpers (beamline, dt, bounds, I/O) ---

void ParallelTracker::selectDT() {
    itsBunch_m->setdT(
            emittingSamplers_m.emissionTimeStep(itsBunch_m->getT()).value_or(dtCurrentTrack_m));
}

void ParallelTracker::changeDT() {
    Inform m("ParallelTracker::changeDT");
    selectDT();
    double newdT = itsBunch_m->getdT();
    for (const auto& pc : itsBunch_m->getParticleContainers()) {
        if (pc) {
            pc->dt = newdT;
        }
    }
    m << level5 << "Changed particle container time step to " << newdT << "." << endl;
}

bool ParallelTracker::stepLogDue() const {
    return options_m.stepInfoFreq > 0
           && itsBunch_m->getGlobalTrackStep() % options_m.stepInfoFreq == 0;
}

opalx::moments::Request ParallelTracker::outputMomentsDue(
        const ParticleContainer_t& pc, bool psDump, bool statDump, bool stepLog) const {
    const opalx::moments::Request log =
            stepLog ? opalx::moments::Request(Moment::MeanMomentum) : opalx::moments::Request{};
    return log | DataSink::requiredMoments(psDump, statDump, pc.hasSpin());
}

std::vector<opalx::moments::BeamMoments> ParallelTracker::momentsOfContainers(
        bool psDump, bool statDump, bool stepLog,
        std::span<const opalx::moments::BeamMoments> computed) const {
    const auto& containers = itsBunch_m->getParticleContainers();
    std::vector<opalx::moments::BeamMoments> moments(containers.size());
    for (size_t i = 0; i < containers.size(); ++i) {
        if (!containers[i] || containers[i]->getTotalNum() == 0) {
            continue;
        }
        const opalx::moments::Request request =
                outputMomentsDue(*containers[i], psDump, statDump, stepLog);
        if (request.empty()) {
            continue;
        }
        const bool reuse = i < computed.size() && computed[i].available().contains(request);
        moments[i]       = reuse ? computed[i] : opalx::moments::compute(*containers[i], request);
    }
    return moments;
}

void ParallelTracker::dumpStats(
        bool psDump, bool statDump, std::span<const opalx::moments::BeamMoments> computed) {
    OPALTimer::Timer myt2;
    const size_t totalAll      = itsBunch_m->getTotalNumAllContainers();
    const long long globalStep = itsBunch_m->getGlobalTrackStep();
    const bool printStepInfo   = stepLogDue();

    if (totalAll == 0 && printStepInfo) {
        *gmsg << level1 << "* " << myt2.time() << " "
              << "Step " << std::setw(6) << globalStep << "; "
              << "   -- no emission yet --     "
              << "t= " << Util::getTimeString(itsBunch_m->getT()) << endl;
    }

    // One request per container covers the step log and every output due now.
    const auto moments = momentsOfContainers(psDump, statDump, printStepInfo, computed);

    bool anyLogged         = false;
    const auto& containers = itsBunch_m->getParticleContainers();
    for (size_t ci = 0; ci < containers.size(); ++ci) {
        const auto& pc = containers[ci];
        if (!pc || pc->getTotalNum() == 0) {
            continue;
        }
        const double sPos = pc->get_sPos();
        if (std::isnan(sPos) || std::isinf(sPos)) {
            throw OpalException(
                    "ParallelTracker::dumpStats()",
                    "invalid path length s for particle container " + std::to_string(ci));
        }
        if (printStepInfo) {
            *gmsg << level1 << "* " << myt2.time() << " "
                  << "Step " << std::setw(6) << globalStep << " "
                  << "container[" << ci << "] ";
            *gmsg << "at " << Util::getLengthString(sPos) << ", ";
            *gmsg << "t= " << Util::getTimeString(itsBunch_m->getT()) << ", "
                  << "E=" << Util::getEnergyString(moments[ci].meanKineticEnergy()) << endl;
        }
        anyLogged = true;
    }

    if (anyLogged || statDump) {
        writePhaseSpace(psDump, statDump, moments);
    }

    if (outputDue(options_m.printRankDistrFreq, itsBunch_m->getGlobalTrackStep())) {
        for (size_t i = 0; i < containers.size(); ++i) {
            const auto& pc = containers[i];
            if (!pc || !itsBunch_m->isPcActive(i)) {
                continue;
            }
            pc->printRankLoadInfo("container[" + std::to_string(i) + "]");
        }
    }
}

void ParallelTracker::writePhaseSpace(
        bool psDump, bool statDump, std::span<const opalx::moments::BeamMoments> computed) {
    Inform m("ParallelTracker::writePhaseSpace");
    Vector_t<double, 3> externalE, externalB;

    const size_t nContainers = itsBunch_m->getNumParticleContainers();
    std::vector<std::array<Vector_t<double, 3>, 2>> fdByContainer(nContainers);
    const auto moments = momentsOfContainers(psDump, statDump, false, computed);

    if (psDump || statDump) {
        for (size_t i = 0; i < nContainers; ++i) {
            auto pc = itsBunch_m->getParticleContainer(i);
            if (!pc || pc->getTotalNum() == 0) {
                fdByContainer[i][0] = Vector_t<double, 3>(0.0);
                fdByContainer[i][1] = Vector_t<double, 3>(0.0);
                continue;
            }
            externalB = Vector_t<double, 3>(0.0);
            externalE = Vector_t<double, 3>(0.0);
            itsOpalBeamline_m.getFieldAt(
                    pc->getRefPartR(), pc->getRefPartP(),
                    itsBunch_m->getT() - 0.5 * itsBunch_m->getdT(), externalE, externalB);
            fdByContainer[i][0] = externalB;
            fdByContainer[i][1] = externalE * Units::Vpm2MVpm;
            m << level5 << "External fields (container " << i << "): externalE = " << externalE
              << ", externalB = " << externalB << endl;
        }
    }

    if (statDump) {
        const int reportedBinCount =
                spaceChargeSolver_m == nullptr ? 1 : spaceChargeSolver_m->reportedBinCount();
        itsDataSink_m->dumpSDDS(*itsBunch_m, fdByContainer, moments, reportedBinCount);
        *gmsg << level3 << "* Wrote beam statistics." << endl;
    }

    if (psDump && itsBunch_m->getTotalNumAllContainers() > 0) {
        itsDataSink_m->dumpH5(*itsBunch_m, fdByContainer, moments);
        *gmsg << level3 << "* Wrote beam phase space." << endl;
    }
}
