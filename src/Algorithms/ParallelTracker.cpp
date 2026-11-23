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

#include "Algorithms/TrackReferenceStep.h"
#include "BasicActions/DumpEMFields.h"

#include "AbsBeamline/RFCavity.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/CavityAutophaser.h"
#include "Algorithms/ReferenceTracer.h"

#include "Physics/Physics.h"
#include "Physics/Units.h"

#include "Processes/GlobalProcesses/GlobalProcess.h"

#include "Steppers/SpinTBMTPusher.h"

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
      beamReferenceMomenta_m(std::move(setup.beamReferenceMomenta)),
      spaceChargeFieldUpdate_m(setup.spaceChargeFieldUpdate),
      itsDataSink_m(setup.dataSink),
      spaceChargeSolver_m(setup.spaceChargeSolver),
      dirichletPlane_m(std::move(setup.dirichletPlane)),
      spaceChargeContainerActivity_m(),
      itsOpalBeamline_m(std::move(setup.beamline)),
      referenceHitMaterial_m(false),
      beamStops_m(std::move(setup.beamStops)),
      sStart_m(setup.sStart),
      globalPhaseShift_m(setup.globalPhaseShift),
      dtCurrentTrack_m(0.0),
      emittingSamplers_m(std::move(setup.emittingSamplers)),
      restarting_m(setup.restarting),
      restartGlobalStep_m(setup.restartGlobalStep),
      restartDt_m(setup.restartDt),
      restartPosition_m(setup.restartPosition),
      cavityPhases_m(std::move(setup.cavityPhases)),
      timeIntegrationTimer1_m(IpplTimings::getTimer("TIntegration1")),
      timeIntegrationTimer2_m(IpplTimings::getTimer("TIntegration2")),
      fieldEvaluationTimer_m(IpplTimings::getTimer("External field eval")),
      referencePassTimer_m(IpplTimings::getTimer("ReferencePass")) {
    for (unsigned int i = 0; i < setup.sStop.size(); ++i) {
        stepSizes_m.push_back(setup.dt[i], setup.sStop[i], setup.maxSteps[i]);
    }

    stepSizes_m.sortAscendingSStop();
    stepSizes_m.resetIterator();
    initializeSpaceChargeContainerActivity();
}

void ParallelTracker::initializeSpaceChargeContainerActivity() {
    const auto& particleContainers = itsBunch_m->getParticleContainers();
    spaceChargeContainerActivity_m.assign(particleContainers.size(), 0);
    for (const auto& container : particleContainers) {
        if (!container) {
            throw OpalException(
                    "ParallelTracker::initializeSpaceChargeContainerActivity",
                    "Cannot track space charge for a null particle container.");
        }
    }
}
// --- execute() ---

/**
 * @copybrief ParallelTracker::execute
 */
void ParallelTracker::execute() {
    if (startTracking()) trackSteps();
}

bool ParallelTracker::startTracking() {
    Inform m("ParallelTracker::startTracking");
    // The reference start and pass may cross a MONITOR; only the step loop records.
    itsOpalBeamline_m.setMonitorsRecordReference(false);
    StepSizeConfig::ResumePosition restartPosition{0, 0};
    if (restarting_m) {
        restartPosition = restartPosition_m;
        stepSizes_m.advanceToResumePosition(restartPosition);
        if (!stepSizes_m.reachedEnd()) {
            const double configuredDt = stepSizes_m.getdT();
            const double tolerance    = 100.0 * std::numeric_limits<double>::epsilon()
                                     * std::max(
                                             {std::abs(configuredDt), std::abs(restartDt_m),
                                              std::numeric_limits<double>::min()});
            if (std::abs(configuredDt - restartDt_m) > tolerance) {
                throw OpalException(
                        "ParallelTracker::startTracking",
                        "restart input changes the time-step schedule at the checkpointed global "
                        "step");
            }
        }
        m << level2 << "Restart global step " << restartGlobalStep_m << " resumes step-size "
          << "segment " << restartPosition.segment << " with "
          << restartPosition.stepsCompletedInSegment << " completed steps in that segment." << endl;
    }

    // PartBunch::resetPcActive() ran in the constructor while containers were still empty
    // (allocate-then-destroy for capacity). Initial particles are loaded later in TrackRun
    // without refreshing these flags, so reference updates must not skip all containers.
    if (!restarting_m || restartPosition.stepsCompletedInSegment == 0) {
        itsBunch_m->resetPcActive();
    }
    activateEmittingContainers(itsBunch_m->getT());

    // Initialize the Boris particle pusher
    BorisPusher pusher;
    m << level3 << "Initialized Boris pusher." << endl;

    // Ensure the time step is positive for the setup phase (buildTracking() did so as well,
    // before the elements were initialised).
    itsBunch_m->setdT(std::abs(itsBunch_m->getdT()));

    // Select the minimal time step from the configuration
    double minTimeStep = stepSizes_m.getMinTimeStep();
    m << level3 << "Selected minimum time step from configuration: " << minTimeStep << endl;

    // A fresh run derives each reference pose from its initial distribution. A restart must keep
    // the per-container reference pose and reference-to-lab transform stored in the checkpoint.
    const auto& particleContainers = itsBunch_m->getParticleContainers();
    if (!restarting_m) {
        CoordinateSystemTrafo beamlineToLab = itsOpalBeamline_m.getCSTrafoLab2Local().inverted();
        for (size_t ci = 0; ci < particleContainers.size(); ++ci) {
            const auto& pc = particleContainers[ci];
            if (!pc || itsBunch_m->isFollower(ci)) {
                continue;  // a decay daughter container takes its root's reference below
            }
            if (!pc->getReference()) {
                throw OpalException(
                        "ParallelTracker::startTracking",
                        "Particle container has null PartData reference during lab-frame init.");
            }
            pc->setToLabTrafo(beamlineToLab);

            // Resolve the reference particle's pose in the beamline frame. Position and
            // momentum follow the same 3-tier rule: the bunch mean when particles already
            // exist; otherwise the input-specified emission offsets (R0, P0) reported by
            // the sampler; otherwise the design pose (lattice origin, beta*gamma along +z).
            Vector_t<double, 3> refR = 0.0;
            Vector_t<double, 3> refP = 0.0;
            if (pc->getTotalNum() > 0) {
                const auto moments =
                        opalx::moments::compute(*pc, Moment::MeanPosition | Moment::MeanMomentum);
                refR = moments.meanR();
                refP = moments.meanP();
            } else {
                // Empty container (e.g. an emitted distribution before its first emission):
                // the bunch mean is undefined, so take the emission offsets from the sampler.
                bool useSamplerPosition         = false;
                bool useSamplerMomentum         = false;
                Vector_t<double, 3> samplerRefR = 0.0;
                Vector_t<double, 3> samplerRefP = 0.0;
                if (ci < emittingSamplers_m.size()) {
                    for (const auto& sampler : emittingSamplers_m[ci]) {
                        if (!sampler) {
                            continue;
                        }
                        if (!useSamplerMomentum && sampler->hasInitialReferenceMomentum()) {
                            samplerRefP        = sampler->getInitialReferenceMomentum();
                            useSamplerMomentum = true;
                        }
                        if (!useSamplerPosition && sampler->hasInitialReferencePosition()) {
                            samplerRefR        = sampler->getInitialReferencePosition();
                            useSamplerPosition = true;
                        }
                    }
                }

                // Position: emission offset R0, else the lattice origin.
                refR = useSamplerPosition ? samplerRefR : Vector_t<double, 3>(0.0);

                // Momentum: emission offset P0, else the design beta*gamma along +z.
                if (useSamplerMomentum) {
                    if (dot(samplerRefP, samplerRefP) <= 0.0) {
                        throw OpalException(
                                "ParallelTracker::startTracking",
                                "Sampler-provided initial reference momentum is zero.");
                    }
                    refP = samplerRefP;
                } else {
                    const PartData& pref = *pc->getReference();
                    const double P0      = pref.getP() / pref.getM();  // beta*gamma from BEAM pc
                    refP                 = Vector_t<double, 3>(0.0, 0.0, P0);
                }
            }

            // A file distribution with PC, ENERGY or GAMMA on BEAM takes that momentum along
            // REFDIRECTION. The mean of a wide angular spread is much shorter than the particle
            // momenta.
            if (ci < beamReferenceMomenta_m.size() && beamReferenceMomenta_m[ci]) {
                refP = *beamReferenceMomenta_m[ci];
            }

            pc->getRefPartR() = beamlineToLab.transformTo(refR);
            pc->getRefPartP() = beamlineToLab.rotateTo(refP);
        }

        m << level4
          << "Transformed reference particle position and momentum to lab frame (all containers)."
          << endl;

        // Integrate reference orbits forward until all container path lengths reach sStart_m.
        findStartPositions(pusher);
        stepSizes_m.advanceToPos(sStart_m);
    } else {
        stepSizes_m.advanceToIndex(restartPosition.segment);
    }
    // Decay daughter containers use their root's reference particle, frame and path length.
    // Copy them after the start search, and after a restart in case the checkpoint differs.
    copyFollowerReferences();

    if (stepSizes_m.reachedEnd()) {
        *gmsg << level1 << "* Checkpoint is already at the end of the configured tracking "
              << "schedule; no integration steps remain." << endl;
        itsOpalBeamline_m.switchElementsOff();
        Kokkos::fence();
        return false;
    }

    // Global spatial bounds: union over all containers
    Vector_t<double, 3> rmin(0.0), rmax(0.0);
    if (itsBunch_m->getTotalNumAllContainers() > 0) {
        computeInitialBounds(rmin, rmax);
        m << level4 << "Initial bunch bounds:\n"
          << "  rmin = " << rmin << "\n"
          << "  rmax = " << rmax << endl;
    }

    // Print reference particle information for all containers
    printInitialContainerRefs(m);

    // Start timing for the reference pass
    IpplTimings::startTimer(referencePassTimer_m);

    // The primary beam (the first container with a reference) traces its reference particle once
    // through the lattice: autophasing, element-length and reach checks.
    const size_t nContainers = itsBunch_m->getNumParticleContainers();
    for (size_t ci = 0; ci < nContainers; ++ci) {
        const auto& pc = itsBunch_m->getParticleContainer(ci);
        if (!pc || !pc->getReference() || itsBunch_m->isFollower(ci)) {
            continue;  // a decay daughter container follows its root's reference
        }
        ReferenceTracer tracer(
                *pc->getReference(), itsOpalBeamline_m, stepSizes_m, minTimeStep, &cavityPhases_m);
        ReferenceTracer::RayState start;
        start.position   = pc->getRefPartR();
        start.momentum   = pc->getRefPartP();
        start.time       = itsBunch_m->getT();
        start.pathLength = pc->get_sPos();
        tracer.execute(start, -rmin(2));
        break;  // only the primary beam
    }
    m << level4 << "Reference pass done." << endl;
    // Stop timing for the reference pass
    IpplTimings::stopTimer(referencePassTimer_m);

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

    referenceHitMaterial_m = false;
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

void ParallelTracker::trackSteps() {
    Inform m("ParallelTracker::trackSteps");
    BorisPusher pusher;

    // The position within the step-size schedule.
    std::size_t segmentIndex = stepSizes_m.getCurrentIndex();
    unsigned long long stepsCompletedInSegment =
            restarting_m ? restartPosition_m.stepsCompletedInSegment : 0;

    bool stopRun = false;
    while (!stepSizes_m.reachedEnd()) {
        const unsigned long long stepsInSegment = stepSizes_m.getNumSteps();
        // Only the first segment of a restart can start with completed steps.
        startSegment(stepsCompletedInSegment > 0, stepsCompletedInSegment);

        while (stepsCompletedInSegment < stepsInSegment) {
            // Resolve emission-dependent dt before the step so the reference, particle
            // transport and self-field evaluation share the same duration.
            changeDT();
            if (!itsBunch_m->anyPcActive()) {
                m << level4 << "No active particle containers; ending inner track segment." << endl;
                break;
            }

            advanceParticles(pusher);
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
                *gmsg << level1 << "* Every beam with a ZSTOP is done; the run ends." << endl;
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

void ParallelTracker::startSegment(
        bool resumesInSegment, unsigned long long stepsCompletedInSegment) {
    Inform m("ParallelTracker::startSegment");
    const unsigned long long step = itsBunch_m->getGlobalTrackStep();
    const unsigned long long trackSteps =
            step + (stepSizes_m.getNumSteps() - stepsCompletedInSegment);
    dtCurrentTrack_m = stepSizes_m.getdT();

    // Select global dt from dtCurrentTrack_m and copy to all container dt views.
    changeDT();
    if (!resumesInSegment) {
        itsBunch_m->resetPcActive();
    }
    activateEmittingContainers(itsBunch_m->getT());
    // Beams that passed their ZSTOP plane stay stopped.
    for (size_t i = 0; i < beamStops_m.size(); ++i) {
        if (beamStops_m[i].done) itsBunch_m->setPcAtSStop(i);
    }

    m << level2 << "Starting track with dt = " << Util::getTimeString(dtCurrentTrack_m)
      << ", track steps = " << step << " to " << trackSteps << "." << endl;
}

/**
 * @copybrief ParallelTracker::advanceParticles
 *
 * @par Frame of reference
 * Particle R and mesh are in the reference frame for the whole step except inside
 * computeSpaceChargeFields() (the solve frame) and computeExternalFields() (each element's local
 * frame, through a FrameScope).
 */
void ParallelTracker::advanceParticles(BorisPusher& pusher) {
    Inform m("ParallelTracker::advanceParticles");
    const long long step = itsBunch_m->getGlobalTrackStep();

    // Historical OPAL evaluated the self-field at R_n and carried the gathered per-particle field
    // through the first half drift. Keep this as an explicit compatibility mode; the default
    // evaluates the field at R_{n+1/2} below.
    if (spaceChargeFieldUpdate_m == SpaceChargeFieldUpdate::PRESTEP) {
        resetFields();
        computeSpaceChargeFields();
    }

    // First half of the time integration
    timeIntegration1(pusher);
    m << level4 << "timeIntegration1 done at step " << step << "." << endl;

    if (spaceChargeFieldUpdate_m == SpaceChargeFieldUpdate::MIDPOINT) {
        resetFields();
        computeSpaceChargeFields();
    }

    // Emit particles in the reference frame with their own fractional timesteps.
    // They experience external fields in their creation step, but no space charge.
    emitFromEmissionSources(itsBunch_m->getT(), itsBunch_m->getdT());

    computeExternalFields();

    // First half of the Thomas-BMT spin rotation (momentum before the kick);
    // timeIntegration2 does the second half after the kick.
    // No-op for containers that have not registered the Pol attribute.
    evolveSpinTBMT();

    // Second half of the time integration
    timeIntegration2(pusher);

    // Apply global processes (e.g. decay) and mark affected particles for deletion
    const size_t nProcessMarked = applyGlobalProcesses(itsBunch_m->getdT());
    if (nProcessMarked > 0) {
        // Decay-like global processes may have created daughters and marked parents
        // invalid. Delete the parents before time/reference updates so they do not
        // contribute to post-step moments.
        deleteInvalidParticles(false, m, "global processes");
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

    if (itsBunch_m->getT() > 0.0) {
        updateReferenceParticles();
        updateRefToLabCSTrafo();
    }

    const auto& particleContainers = itsBunch_m->getParticleContainers();
    for (size_t i = 0; i < particleContainers.size(); ++i) {
        const auto& pc = particleContainers[i];
        if (!pc || !itsBunch_m->isPcActive(i)) {
            continue;
        }
        m << level4 << "Current path length (container " << i << ") is " << pc->get_sPos() << "."
          << endl;
    }
}

std::vector<opalx::moments::BeamMoments> ParallelTracker::removeLostParticles(
        bool psDump, bool statDump) {
    Inform m("ParallelTracker::removeLostParticles");
    const auto& particleContainers = itsBunch_m->getParticleContainers();

    // Every frame change of this step is closed: moments, output and the next step read
    // reference-frame data.
    for (const auto& pc : particleContainers)
        if (pc) pc->requireReferenceFrame("ParallelTracker::removeLostParticles");

    // Enforce the transverse aperture of every nearby element (collective call).
    size_t nStepMarked = applyElementApertures();
    if (nStepMarked > 0) {
        m << level2 << "Marked " << nStepMarked
          << " particles outside a transverse aperture for deletion." << endl;
    }

    // The outlier box and the outputs of this step read one computation per container,
    // unless particles are deleted in between.
    const double sigmas = options_m.boundpDestroy;
    std::vector<opalx::moments::BeamMoments> stepMoments(particleContainers.size());
    for (size_t i = 0; i < particleContainers.size(); ++i) {
        const auto& pc = particleContainers[i];
        if (!pc || !itsBunch_m->isPcActive(i) || sigmas <= 0.0 || pc->getTotalNum() == 0) {
            continue;
        }
        stepMoments[i] = opalx::moments::compute(
                *pc, Moment::MeanPosition | Moment::Spread
                             | outputMomentsDue(*pc, psDump, statDump, stepLogDue()));
        const Vector_t<double, 3> halfWidth = sigmas * stepMoments[i].rmsR();
        const size_t nMarked                = pc->markParticlesOutside(
                stepMoments[i].meanR() - halfWidth, stepMoments[i].meanR() + halfWidth);
        nStepMarked += nMarked;
        if (nMarked > 0) {
            m << level2 << "Marked " << nMarked << " particles outside " << sigmas
              << "-sigma boundary for deletion (container " << i << ")." << endl;
        }
    }

    // Single flush of all marks accumulated this step.
    if (nStepMarked > 0) {
        deleteInvalidParticles(
                true, m, "transverse aperture / " + std::to_string(sigmas) + "-sigma boundary");
        stepMoments.assign(particleContainers.size(), {});
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

    ippl::Comm->allreduce(referenceHitMaterial_m, 1, std::logical_or<bool>());
    if (referenceHitMaterial_m) {
        *gmsg << level1 << "* A reference particle hit an element; the run ends." << endl;
    }
    stop.runFinished = referenceHitMaterial_m || beam_stop::runFinished(beamStops_m);
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
        const std::string checkpointPath =
                CheckpointFile::defaultPath(OpalData::getInstance()->getInputBasename());
        CheckpointFile::write(
                checkpointPath, *itsBunch_m, nextSegment, nextSegmentSteps, beamStops_m,
                cavityPhases_m);
        m << level2 << "Wrote checkpoint '" << checkpointPath << "' after global step "
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

    if (spaceChargeSolver_m != nullptr) {
        m << level2 << "Total FieldSolver calls: " << spaceChargeSolver_m->backendSolveCount()
          << endl;
    }
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

/**
 * @copybrief ParallelTracker::timeIntegration1
 */
void ParallelTracker::timeIntegration1(BorisPusher& pusher) {
    Inform m("ParallelTracker::timeIntegration1");
    IpplTimings::startTimer(timeIntegrationTimer1_m);
    const size_t n = itsBunch_m->getNumParticleContainers();
    for (size_t i = 0; i < n; ++i) {
        if (!itsBunch_m->isPcActive(i)) {
            continue;
        }
        auto pc = itsBunch_m->getParticleContainer(i);
        if (!pc) {
            continue;
        }
        pushParticles(pusher, *pc);
    }
    IpplTimings::stopTimer(timeIntegrationTimer1_m);
    m << level4 << "Push particles done for all containers." << endl;
}

/**
 * @copybrief ParallelTracker::timeIntegration2
 */
void ParallelTracker::timeIntegration2(BorisPusher& pusher) {
    Inform m("ParallelTracker::timeIntegration2");

    IpplTimings::startTimer(timeIntegrationTimer2_m);
    const size_t n = itsBunch_m->getNumParticleContainers();
    for (size_t i = 0; i < n; ++i) {
        if (!itsBunch_m->isPcActive(i)) {
            continue;
        }
        auto pc = itsBunch_m->getParticleContainer(i);
        if (!pc) {
            continue;
        }
        kickParticles(pusher, *pc);
        if (pc->hasSpin()) {
            evolveSpinHalfStep(*pc);  // second half of the spin rotation, momentum after the kick
        }
        pushParticles(pusher, *pc);
        pc->dt = itsBunch_m->getdT();
    }
    m << level4 << "Kick/push particles done for all containers." << endl;

    IpplTimings::stopTimer(timeIntegrationTimer2_m);
}

ParticleFrame ParallelTracker::makeSpaceChargeFrame(opalx::spacecharge::SolveFrame request) const {
    using opalx::spacecharge::SolveFrame;
    if (request == SolveFrame::Reference) return ParticleFrame();
    const auto primary = itsBunch_m->getParticleContainer();
    primary->requireReferenceFrame("ParallelTracker::makeSpaceChargeFrame");
    const bool centroid = request == SolveFrame::Centroid;
    const auto moments  = opalx::moments::compute(
            *primary, centroid ? Moment::MeanPosition | Moment::MeanMomentum
                                : opalx::moments::Request(Moment::MeanMomentum));
    return opalx::spacecharge::makeSolveFrame(
            request, centroid ? moments.meanR() : Vector_t<double, 3>(0.0), moments.meanP(),
            primary->getRefPartP(), primary->get_sPos());
}

ParallelTracker::SpaceChargeEmissionProgress ParallelTracker::spaceChargeEmissionProgress() const {
    SpaceChargeEmissionProgress progress;
    const double currentTime = itsBunch_m->getT();

    // While emission is active, CartesianPIC3D stretches its beam-frame mesh over the full source
    // pulse length, matching old OPAL. For example, 5 percent emission produces approximately a
    // factor-20 longitudinal stretch. Without it, early charge is compressed onto an artificially
    // short mesh and can receive excessive self-field kicks, including kicks back into the source.
    // The least-complete active source determines the required stretch.
    for (const auto& samplers : emittingSamplers_m) {
        for (const auto& sampler : samplers) {
            if (!sampler || sampler->isEmissionDone(currentTime)) {
                continue;
            }
            const double samplerFraction = sampler->getEmittedFraction();
            if (samplerFraction < 1.0) {
                progress.active   = true;
                progress.fraction = std::min(progress.fraction, samplerFraction);
            }
        }
    }

    if (options_m.aggressiveStateSync) {
        // Preserve legacy state convergence across ranks: any active source keeps stretching
        // enabled, and the minimum emitted fraction supplies the most conservative domain length.
        bool globallyActive = progress.active;
        ippl::Comm->allreduce(progress.active, globallyActive, 1, std::logical_or<bool>());
        const double localFraction =
                progress.active ? std::clamp(progress.fraction, 0.0, 1.0) : 1.0;
        double globalFraction = localFraction;
        ippl::Comm->allreduce(localFraction, globalFraction, 1, std::less<double>());
        progress.active   = globallyActive;
        progress.fraction = globallyActive ? globalFraction : 1.0;
    } else {
        progress.fraction = progress.active ? std::clamp(progress.fraction, 0.0, 1.0) : 1.0;
    }
    return progress;
}

/**
 * @copybrief ParallelTracker::computeSpaceChargeFields
 *
 * @par Frame of reference
 * The algorithm requests its solve frame (SpaceChargeAlgorithm::solveFrame), which is built here
 * from the primary container's moments. The particles are in their reference frame on entry and
 * again on return; the algorithm enters the solve frame and leaves it.
 */
void ParallelTracker::computeSpaceChargeFields() {
    Inform m("ParallelTracker::computeSpaceChargeFields");
    if (spaceChargeSolver_m == nullptr) {
        throw OpalException(
                "ParallelTracker::computeSpaceChargeFields",
                "No space-charge solver is available. Use TYPE=NONE for a configured no-op "
                "solver.");
    }

    const size_t totalParticles = itsBunch_m->getTotalNumAllContainers();
    if (totalParticles <= static_cast<size_t>(options_m.minBinEmitted)) {
        m << level4
          << "Skipping space charge until more than MINBINEMITTED=" << options_m.minBinEmitted
          << " particles are present (total=" << totalParticles << ")." << endl;
        return;
    }

    const auto frame    = makeSpaceChargeFrame(spaceChargeSolver_m->solveFrame());
    const auto emission = spaceChargeEmissionProgress();

    for (std::size_t index = 0; index < spaceChargeContainerActivity_m.size(); ++index) {
        spaceChargeContainerActivity_m[index] = itsBunch_m->isPcActive(index) ? 1 : 0;
    }

    using namespace opalx::spacecharge;
    const std::size_t step = static_cast<std::size_t>(itsBunch_m->getGlobalTrackStep());

    // Algorithms borrow the stable containers at construction; this context supplies step state.
    SpaceChargeStepState stepState{
            step,
            itsBunch_m->getT(),
            itsBunch_m->getdT(),
            emission.active,
            emission.fraction,
            ippl::Comm->size(),
            frame};
    SpaceChargeSolveContext context(spaceChargeContainerActivity_m, std::move(stepState));
    spaceChargeSolver_m->solve(context);
    for (const auto& pc : itsBunch_m->getParticleContainers())
        if (pc) pc->requireReferenceFrame("ParallelTracker::computeSpaceChargeFields");
    m << level3 << "Compute space-charge fields done." << endl;
}

/**
 * @copybrief ParallelTracker::computeExternalFields
 */
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
    const size_t nSourceMarked = markBackwardParticlesAtSourcePlane();
    if (nSourceMarked > 0) {
        deleteInvalidParticles(
                true, msg, "backward source-plane particles during external-field evaluation");
    }

    forEachElementInBunchFrame(
            [](const ElementBase& element) {
                return element.fieldFrame() == ElementFrame::Local;
            },
            [this](const std::shared_ptr<ElementBase>& element,
                   const std::shared_ptr<ParticleContainer_t>& pc) {
                element->apply(pc, stepTime());
            });

    IpplTimings::stopTimer(fieldEvaluationTimer_m);
}

/**
 * @copybrief ParallelTracker::forEachElementInBunchFrame
 */
void ParallelTracker::forEachElementInBunchFrame(
        const std::function<bool(const ElementBase&)>& selects,
        const std::function<void(
                const std::shared_ptr<ElementBase>&, const std::shared_ptr<ParticleContainer_t>&)>&
                func) {
    for (size_t ci = 0; ci < itsBunch_m->getNumParticleContainers(); ++ci) {
        const auto pc = itsBunch_m->getParticleContainer(ci);
        if (!pc) continue;
        const ParticleContainer_t::FrameScope scope(*pc, pc->frame());
        for (const auto& element : candidateElements(ci)) {
            if (!selects(*element)) continue;
            pc->enterFrame(elementFrame(element, *pc));
            func(element, pc);
        }
    }
}

/**
 * @copybrief ParallelTracker::candidateElements
 */
std::vector<std::shared_ptr<ElementBase>> ParallelTracker::candidateElements(size_t ci) {
    const auto pc = itsBunch_m->getParticleContainer(ci);
    if (!itsBunch_m->isPcActive(ci) || !pc) {
        return {};
    }
    // A container without its own or its root's reference particle has no bunch frame.
    const auto root = itsBunch_m->getParticleContainer(
            itsBunch_m->isFollower(ci) ? itsBunch_m->getPcRoot(ci) : ci);
    if (!root || !root->getReference()) {
        return {};
    }
    pc->requireReferenceFrame("ParallelTracker::candidateElements");

    // Bunch bounds for this container; an empty container queries the point at the reference
    // particle.
    Vector_t<double, 3> rmin(0.0), rmax(0.0);
    if (pc->getTotalNum() > 0) {
        const auto bounds = opalx::moments::compute(*pc, Moment::Bounds);
        rmin              = bounds.minR();
        rmax              = bounds.maxR();
    }

    // The bunch box, extended forward by the distance one step can cover, c dt: a MONITOR
    // records particles that cross its plane during the next step. Elements come in lattice
    // order.
    Vector_t<double, 3> upper = rmax;
    upper(2) += Physics::c * itsBunch_m->getdT();
    return itsOpalBeamline_m.getElementMapping().query(rmin, upper, pc->getToLabTrafo());
}

ParticleFrame ParallelTracker::beamlineFrame(const ParticleContainer_t& pc) const {
    return ParticleFrame(
            ParticleFrame::Kind::Beamline,
            itsOpalBeamline_m.getCSTrafoLab2Local() * pc.getToLabTrafo());
}

ParticleFrame ParallelTracker::elementFrame(
        const std::shared_ptr<ElementBase>& element, const ParticleContainer_t& pc) const {
    return ParticleFrame(
            ParticleFrame::Kind::Element,
            itsOpalBeamline_m.getMisalignment(element)
                    * (itsOpalBeamline_m.getCSTrafoLab2Local(element) * pc.getToLabTrafo()),
            element.get());
}

/**
 * @copybrief ParallelTracker::applyElementApertures
 */
size_t ParallelTracker::applyElementApertures() {
    size_t localMarked = 0;
    for (size_t ci = 0; ci < itsBunch_m->getNumParticleContainers(); ++ci) {
        const auto pc = itsBunch_m->getParticleContainer(ci);
        for (const auto& element : candidateElements(ci)) {
            // The check reads local positions in its kernel; the particles stay where they are.
            // Marks are flushed by the single deletion at the end of the step.
            const auto toLocal = pc->frame().to(elementFrame(element, *pc)).rigid();
            localMarked += element->markOutsideAperture(pc, toLocal);
        }
    }

    size_t globalMarked = 0;
    ippl::Comm->allreduce(localMarked, globalMarked, 1, std::plus<size_t>());
    return globalMarked;
}

/**
 * @copybrief ParallelTracker::emitFromEmissionSources
 */
void ParallelTracker::emitFromEmissionSources(double t, double dt) {
    const auto& containers = itsBunch_m->getParticleContainers();
    for (size_t ci = 0; ci < containers.size(); ++ci) {
        const auto& pc = containers[ci];
        if (!pc) {
            continue;
        }
        if (itsBunch_m->pcAtSStop(ci)) {
            continue;
        }

        // Record the extent of the position array before emission. If an internal resize
        // (Kokkos::realloc) happens during emission, the extent of R will change and we can
        // flag this as an error.
        const size_t extentBeforeEmission = pc->R.size();

        // Samplers append new particles in the beamline frame and do not touch the others, so
        // only the new particles change frame.
        const size_t nBefore = pc->getLocalNum();
        if (ci < emittingSamplers_m.size()) {
            for (const auto& sampler : emittingSamplers_m[ci]) {
                if (sampler) {
                    sampler->emitParticles(t, dt);
                }
            }
        }
        pc->adoptParticles(beamlineFrame(*pc), nBefore);

        pc->setM(pc->getMassPerParticle());
        pc->setQ(pc->getChargePerParticle());

        // Sanity guard: the total number of macroparticles in the bunch must
        // never exceed the globally configured BEAM::NALLOC value. Overshooting
        // this limit would trigger internal reallocations in the particle
        // container and silently drop already-tracked particles/delete their data in the particle
        // attributes. This is only a check for the number of local particles.
        const size_t extentAfterEmission = pc->R.size();
        if (extentAfterEmission != extentBeforeEmission) {
            throw OpalException(
                    "ParallelTracker::emitFromEmissionSources",
                    "Local particle storage was resized during emission (likely due to "
                    "over-emission causing a Kokkos::realloc). This leads to loss of "
                    "previously tracked particles. Please increase the total number of "
                    "macroparticles or adjust NPARTDIST / the emission profile. If you are "
                    "using emission sources, please check the emission profile and adjust "
                    "the number of particles emitted.");
        }
    }
    itsBunch_m->refreshPcActiveAfterEmit();
}

size_t ParallelTracker::applyGlobalProcesses(double dt) {
    const size_t nContainers        = itsBunch_m->getNumParticleContainers();
    const long long globalTrackStep = itsBunch_m->getGlobalTrackStep();
    size_t nMarked                  = 0;

    for (size_t ci = 0; ci < nContainers; ++ci) {
        auto pc = itsBunch_m->getParticleContainer(ci);
        if (!pc) {
            continue;
        }
        // Decay writes daughters in the shared reference frame of parent and daughter.
        pc->requireReferenceFrame("ParallelTracker::applyGlobalProcesses");
        const auto& processes = pc->getGlobalProcesses();
        if (processes.empty()) {
            continue;
        }
        for (const auto& process : processes) {
            if (process) {
                nMarked += process->apply(*pc, dt, globalTrackStep, ci);
            }
        }
    }
    // A decay may have put the first particles into an empty daughter container. Switch it on
    // now, so its particles take part in the next step's first push. nMarked is the global
    // count, so every rank does the same.
    if (nMarked > 0) {
        itsBunch_m->refreshPcActiveAfterEmit();
    }
    return nMarked;
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

void ParallelTracker::deleteInvalidParticles(
        bool activeOnly, Inform& m, const std::string& reason) {
    const size_t nContainers = itsBunch_m->getNumParticleContainers();

    for (size_t ci = 0; ci < nContainers; ++ci) {
        if (activeOnly && !itsBunch_m->isPcActive(ci)) {
            continue;
        }
        auto pc = itsBunch_m->getParticleContainer(ci);
        if (!pc) {
            continue;
        }

        const size_t nContainerDeleted = pc->deleteInvalidParticles();
        if (nContainerDeleted > 0) {
            m << level2 << "Deleted " << nContainerDeleted << " particles marked by " << reason
              << ", " << pc->getTotalNum() << " remaining (container " << ci << ")." << endl;
        }
    }
}

size_t ParallelTracker::markBackwardParticlesAtSourcePlane() {
    /// \todo this function should probably be integrated as a GunSource element similar to old
    /// OPAL.
    if (spaceChargeSolver_m == nullptr) {
        return 0;
    }

    using namespace opalx::spacecharge;
    if (dirichletPlane_m.kind == DirichletPlaneType::None) {
        return 0;
    }

    const double sourcePlaneZ = dirichletPlane_m.planeZ;
    // Legacy OPAL's SOURCE element is 5 cm long and is shifted upstream from ELEMEDGE.
    // Source::apply deletes only once a particle crosses the element-local entrance plane
    // (Rz <= 0), not when it crosses the cathode/image plane at ELEMEDGE.
    constexpr double legacySourceLength = 0.05;
    const double sourceLossPlaneZ       = sourcePlaneZ - legacySourceLength;

    size_type localTotalMarked = 0;
    const size_t nContainers   = itsBunch_m->getNumParticleContainers();
    for (size_t ci = 0; ci < nContainers; ++ci) {
        if (!itsBunch_m->isPcActive(ci)) {
            continue;
        }

        auto pc = itsBunch_m->getParticleContainer(ci);
        if (!pc || pc->getLocalNum() == 0) {
            continue;
        }

        // Read source-frame coordinates in the kernel; the particles stay where they are.
        const auto toSource = pc->frame().to(beamlineFrame(*pc)).rigid();

        auto Rview   = pc->R.getView();
        auto Pview   = pc->P.getView();
        auto invalid = pc->InvalidMask.getView();

        size_type localMarked  = 0;
        const size_type nLocal = static_cast<size_type>(pc->getLocalNum());
        Kokkos::parallel_reduce(
                "ParallelTracker::markBackwardParticlesAtSourcePlane", nLocal,
                KOKKOS_LAMBDA(const size_type i, size_type& count) {
                    const Vector_t<double, 3> localR = toSource.pointTo(Rview(i));
                    const Vector_t<double, 3> localP = toSource.vectorTo(Pview(i));
                    const bool backwards   = localR[2] <= sourceLossPlaneZ && localP[2] < 0.0;
                    const bool newlyMarked = backwards && !invalid(i);
                    invalid(i)             = invalid(i) || backwards;
                    count += newlyMarked ? 1 : 0;
                },
                localMarked);
        Kokkos::fence();

        localTotalMarked += localMarked;
    }

    size_type globalTotalMarked = 0;
    ippl::Comm->allreduce(localTotalMarked, globalTotalMarked, 1, std::plus<size_type>());
    return static_cast<size_t>(globalTotalMarked);
}

/**
 * @copybrief ParallelTracker::resetFields
 */
void ParallelTracker::resetFields() {
    const size_t n = itsBunch_m->getNumParticleContainers();
    for (size_t i = 0; i < n; ++i) {
        if (!itsBunch_m->isPcActive(i)) {
            continue;
        }
        auto pc = itsBunch_m->getParticleContainer(i);
        if (!pc) {
            continue;
        }
        pc->E = 0;
        pc->B = 0;
    }
}

/**
 * @brief Boris position push in unitless coordinates (per-particle dt via pusher internals).
 * @param pusher Boris pusher.
 * @param pc     Non-null target particle container.
 */
void ParallelTracker::pushParticles(
        const BorisPusher& pusher, PartBunch_t::ParticleContainer_t& pc) {
    // Per-particle dt is used so that newly emitted particles with fractional dt
    // (sampled during emission) are pushed proportionally to their sub-timestep fraction.
    pc.switchToUnitlessPositions();

    auto Rview = pc.R.getView();
    auto Pview = pc.P.getView();

    Kokkos::parallel_for(
            "pushParticles", pc.getLocalNum(), KOKKOS_LAMBDA(const size_t i) {
                // Half drift: x_{n+1/2} = x_n + (dt/2) * v; unitless form via pusher.push(..., 0).
                // TODO: verify sign convention for half-step push.
                Vector_t<double, 3> x = Rview(i);
                pusher.push(
                        x, Pview(i),
                        0);  // this 0 is "dt" that is not used with unitless positions!
                Rview(i) = x;
            });

    pc.switchOffUnitlessPositions();
    Kokkos::fence();
    ippl::Comm->barrier();
}

/**
 * @brief Boris velocity kick from E and B using per-particle dt.
 * @param pusher Boris pusher.
 * @param pc     Non-null target particle container.
 */
void ParallelTracker::kickParticles(
        const BorisPusher& pusher, PartBunch_t::ParticleContainer_t& pc) {
    Inform m("ParallelTracker::kickParticles");

    auto Pview  = pc.P.getView();
    auto dtview = pc.dt.getView();
    auto Efview = pc.E.getView();
    auto Bfview = pc.B.getView();
    m << level5 << "Got particle views for kick operation." << endl;
    // Mass (eV) and charge (proton charges) from this container's reference particle,
    // passed explicitly into BorisPusher::kick for GPU-safe kernels.

    const PartData& ref = *pc.getReference();
    const double mass   = ref.getM();
    const double charge = ref.getQ();
    Kokkos::parallel_for(
            "kickParticles", pc.getLocalNum(), KOKKOS_LAMBDA(const size_t i) {
                // Boris kick: Birdsall & Langdon (1985), ch. 4-4 (non-relativistic) and ch. 15-4
                // (relativistic); scaled P = (v/c)*gamma, mass in rest energy units.
                Vector_t<double, 3> p = Pview(i);
                // TODO: consider dropping unused R/dt arguments from kick when API allows.
                pusher.kick(0, p, Efview(i), Bfview(i), dtview(i), mass, charge);
                Pview(i) = p;
            });
    /*
    Wait until everyone completed the kick operation before proceeding. For now,
    this is just a precaution and could be removed for a small performance gain
    (technically, these shouldn't be necessary!).
    */
    Kokkos::fence();
    ippl::Comm->barrier();

    m << level5 << "Completed parallel kick operation." << endl;
}

void ParallelTracker::evolveSpinTBMT() {
    const size_t n = itsBunch_m->getNumParticleContainers();
    for (size_t i = 0; i < n; ++i) {
        if (!itsBunch_m->isPcActive(i)) {
            continue;
        }
        auto pc = itsBunch_m->getParticleContainer(i);
        if (!pc || !pc->hasSpin()) {
            continue;
        }
        evolveSpinHalfStep(*pc);
    }
}

void ParallelTracker::evolveSpinHalfStep(PartBunch_t::ParticleContainer_t& pc) {
    const PartData& ref = *pc.getReference();
    const double mass   = ref.getM();
    const double charge = ref.getQ();
    const double anom   = ref.getAnomaly();

    auto Polview = pc.Pol.getView();
    auto Pview   = pc.P.getView();
    auto Efview  = pc.E.getView();
    auto Bfview  = pc.B.getView();
    auto dtview  = pc.dt.getView();

    SpinTBMTPusher spinPusher;
    Kokkos::parallel_for(
            "evolveSpinHalfStep", pc.getLocalNum(), KOKKOS_LAMBDA(const size_t j) {
                spinPusher.evolve(
                        Polview(j), Pview(j), Efview(j), Bfview(j), 0.5 * dtview(j), mass, charge,
                        anom);
            });
    Kokkos::fence();
}

// --- Helpers (beamline, dt, bounds, I/O) ---

/**
 * @copybrief ParallelTracker::selectDT
 */
void ParallelTracker::selectDT() {
    double selectedDt        = dtCurrentTrack_m;
    double emissionDt        = std::numeric_limits<double>::max();
    bool hasEmissionDt       = false;
    const double currentTime = itsBunch_m->getT();

    for (const auto& samplers : emittingSamplers_m) {
        for (const auto& sampler : samplers) {
            if (!sampler || sampler->isEmissionDone(currentTime)) {
                continue;
            }
            const double samplerDt = sampler->getEmissionTimeStep();
            if (samplerDt > 0.0) {
                emissionDt    = std::min(emissionDt, samplerDt);
                hasEmissionDt = true;
            }
        }
    }

    if (hasEmissionDt) {
        selectedDt = emissionDt;
    }
    itsBunch_m->setdT(selectedDt);
}

/**
 * @copybrief ParallelTracker::changeDT
 */
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

/**
 * @copybrief ParallelTracker::activateEmittingContainers
 */
void ParallelTracker::activateEmittingContainers(double t) {
    for (size_t i = 0; i < emittingSamplers_m.size(); ++i) {
        for (const auto& sampler : emittingSamplers_m[i]) {
            if (sampler && !sampler->isEmissionDone(t)) {
                itsBunch_m->setPcActive(i);
                break;
            }
        }
    }
}

/**
 * @copybrief ParallelTracker::computeInitialBounds
 */
void ParallelTracker::computeInitialBounds(Vector_t<double, 3>& rmin, Vector_t<double, 3>& rmax) {
    std::vector<const ParticleContainer_t*> containers;
    for (const auto& pc : itsBunch_m->getParticleContainers()) {
        if (pc) {
            containers.push_back(pc.get());
        }
    }
    // Empty containers contribute nothing; without any particle the box is a point at the origin.
    const auto bounds = opalx::moments::compute(containers, Moment::Bounds);
    rmin              = bounds.minR();
    rmax              = bounds.maxR();
}

/**
 * @copybrief ParallelTracker::printInitialContainerRefs
 */
void ParallelTracker::printInitialContainerRefs(Inform& m) const {
    const auto& particleContainers = itsBunch_m->getParticleContainers();
    for (size_t i = 0; i < particleContainers.size(); ++i) {
        const auto& pc = particleContainers[i];
        if (!pc) {
            m << level3 << "ParallelTrack: container " << i << " is null." << endl;
            continue;
        }
        if (pc->getTotalNum() > 0) {
            // Collective; evaluated on every rank before the rank-filtered log line.
            const double meanPz = opalx::moments::compute(*pc, Moment::MeanMomentum).meanP()(2);
            m << level3 << "ParallelTrack (container " << i << "): momentum z = " << meanPz << "\n"
              << "RefPartR (container " << i << ") = " << pc->getRefPartR() << "\n"
              << "RefPartP (container " << i << ") = " << pc->getRefPartP() << endl;
        } else {
            m << level3 << "ParallelTrack: container " << i
              << " empty; total particles = " << itsBunch_m->getTotalNumAllContainers() << "\n"
              << "RefPartR (container " << i << ") = " << pc->getRefPartR() << "\n"
              << "RefPartP (container " << i << ") = " << pc->getRefPartP() << endl;
        }
    }
}

/**
 * @copybrief ParallelTracker::updateReferenceParticles
 */
void ParallelTracker::updateReferenceParticles() {
    const double dt = std::min(itsBunch_m->getT(), itsBunch_m->getdT());

    const size_t n = itsBunch_m->getNumParticleContainers();
    for (size_t i = 0; i < n; ++i) {
        // A decay daughter container has no reference particle of its own. Its root keeps
        // moving while any of its daughter containers is tracked, even when the root is empty.
        if (itsBunch_m->isFollower(i) || !itsBunch_m->referenceMoves(i)) {
            continue;
        }
        auto pcPtr = itsBunch_m->getParticleContainer(i);
        if (!pcPtr) {
            continue;
        }
        auto& pc                = *pcPtr;
        const PartData& refKick = *pc.getReference();
        const auto advance      = track_reference::advanceInBeamline;
        const auto end =
                advance(itsOpalBeamline_m, refKick, {pc.getRefPartR(), pc.getRefPartP()}, dt,
                        itsBunch_m->getT(), true);
        pc.getRefPartR() = end.position;
        pc.getRefPartP() = end.momentum;
        if (end.hitMaterial) {
            *gmsg << level1 << "The reference particle hit an element" << endl;
            referenceHitMaterial_m = true;
        }
    }
}

/**
 * @copybrief ParallelTracker::updateRefToLabCSTrafo
 */
void ParallelTracker::updateRefToLabCSTrafo() {
    // Map the lab-frame reference position and momentum into the current bunch frame, then
    // move and turn the bunch frame onto them. Momentum is only rotated, not translated.
    const double bunchDT = itsBunch_m->getdT();
    const size_t n       = itsBunch_m->getNumParticleContainers();
    std::vector<CoordinateSystemTrafo> updates(n);  // identity unless the container moved
    for (size_t i = 0; i < n; ++i) {
        if (itsBunch_m->isFollower(i) || !itsBunch_m->referenceMoves(i)) {
            continue;
        }
        auto pc = itsBunch_m->getParticleContainer(i);
        if (pc) {
            updates[i] = pc->updateRefToLabCSTrafo(bunchDT);
        }
    }
    // Decay daughter containers take the same frame change as their root, active or not, so
    // daughters created in this step end up in the root's new frame too.
    for (size_t i = 0; i < n; ++i) {
        if (!itsBunch_m->isFollower(i)) {
            continue;
        }
        const size_t root = itsBunch_m->getPcRoot(i);
        auto pc           = itsBunch_m->getParticleContainer(i);
        auto rootPc       = itsBunch_m->getParticleContainer(root);
        if (pc && rootPc) {
            pc->followFrameOf(*rootPc, updates[root]);
        }
    }
}

/**
 * @copybrief ParallelTracker::copyFollowerReferences
 */
void ParallelTracker::copyFollowerReferences() {
    const size_t n = itsBunch_m->getNumParticleContainers();
    for (size_t i = 0; i < n; ++i) {
        if (!itsBunch_m->isFollower(i)) {
            continue;
        }
        auto pc     = itsBunch_m->getParticleContainer(i);
        auto rootPc = itsBunch_m->getParticleContainer(itsBunch_m->getPcRoot(i));
        if (pc && rootPc) {
            pc->copyReferenceOf(*rootPc);
        }
    }
}

/**
 * @copybrief ParallelTracker::findStartPositions
 */
void ParallelTracker::findStartPositions(const BorisPusher& pusher) {
    // Primary container (index 0) drives segment advances.
    auto primary = itsBunch_m->getParticleContainer(0);

    if (sStart_m <= primary->get_sPos()) {
        return;
    }

    StepSizeConfig stepSizesCopy(stepSizes_m);

    double t = 0.0;
    itsBunch_m->setT(t);

    dtCurrentTrack_m = stepSizesCopy.getdT();
    selectDT();

    const auto& containers = itsBunch_m->getParticleContainers();

    // The last, fractional step takes every reference to sTarget; the bunch time advances by the
    // primary's share.
    const auto applyLastStep = [&](double sTarget, double speed) {
        const double tauPrimary = (sTarget - primary->get_sPos()) / speed;
        itsBunch_m->setT(itsBunch_m->getT() + tauPrimary);
        for (size_t i = 0; i < containers.size(); ++i) {
            const auto& pc = containers[i];
            if (!pc || itsBunch_m->isFollower(i)) {
                continue;
            }
            Vector_t<double, 3> pv =
                    pc->getRefPartP() * Physics::c / Util::getGamma(pc->getRefPartP());
            double speed_i = euclidean_norm(pv);
            double tau_i   = (sTarget - pc->get_sPos()) / speed_i;
            pc->applyFractionalStep(pusher, tau_i, sStart_m);
        }
    };

    while (true) {
        autophaseCavities(pusher);

        t += itsBunch_m->getdT();
        itsBunch_m->setT(t);

        std::vector<Vector_t<double, 3>> oldRs(containers.size());
        for (size_t i = 0; i < containers.size(); ++i) {
            if (containers[i]) {
                oldRs[i] = containers[i]->getRefPartR();
            }
        }

        updateReferenceParticles();

        for (size_t i = 0; i < containers.size(); ++i) {
            if (!containers[i] || itsBunch_m->isFollower(i)) {
                continue;
            }
            Vector_t<double, 3> dR = containers[i]->getRefPartR() - oldRs[i];
            containers[i]->set_sPos(containers[i]->get_sPos() + euclidean_norm(dR));
        }

        Vector_t<double, 3> tmp =
                primary->getRefPartP() * Physics::c / Util::getGamma(primary->getRefPartP());
        double speed = euclidean_norm(tmp);

        if (primary->get_sPos() > stepSizesCopy.getSStop()) {
            ++stepSizesCopy;

            if (stepSizesCopy.reachedEnd()) {
                --stepSizesCopy;
                applyLastStep(stepSizesCopy.getSStop(), speed);
                break;
            }

            dtCurrentTrack_m = stepSizesCopy.getdT();
            selectDT();
        }

        if (std::abs(primary->get_sPos() - sStart_m) <= 0.5 * itsBunch_m->getdT() * speed) {
            applyLastStep(sStart_m, speed);
            break;
        }
    }

    changeDT();
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

/**
 * @copybrief ParallelTracker::dumpStats
 */
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

/**
 * @copybrief ParallelTracker::writePhaseSpace
 */
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

// --- Autophasing (RF and traveling-wave cavities) ---

/**
 * @copybrief ParallelTracker::autophaseCavities
 */
void ParallelTracker::autophaseCavities(const BorisPusher& pusher) {
    const PartData& ref = *itsBunch_m->getParticleContainer()->getReference();
    double t            = itsBunch_m->getT();
    Vector_t<double, 3> nextR =
            itsBunch_m->getParticleContainer()->getRefPartR() / (Physics::c * itsBunch_m->getdT());
    pusher.push(nextR, itsBunch_m->getParticleContainer()->getRefPartP(), itsBunch_m->getdT());
    nextR *= Physics::c * itsBunch_m->getdT();

    auto elementSet = itsOpalBeamline_m.getElements(nextR);
    for (auto element : elementSet) {
        // A TravelingWave is an RFCavity.
        if (element->getType() != ElementType::TRAVELINGWAVE
            && element->getType() != ElementType::RFCAVITY) {
            continue;
        }
        const RFCavity* cavity = static_cast<const RFCavity*>(element.get());
        if (!cavity->getAutophaseVeto()) {
            CavityAutophaser ap(ref, element, &cavityPhases_m);
            ap.getPhaseAtMaxEnergy(
                    itsOpalBeamline_m.transformToLocalCS(
                            element, itsBunch_m->getParticleContainer()->getRefPartR()),
                    itsOpalBeamline_m.rotateToLocalCS(
                            element, itsBunch_m->getParticleContainer()->getRefPartP()),
                    t, itsBunch_m->getdT());
        }
    }
}
