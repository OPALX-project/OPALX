/**
 * @file ParallelTracker.h
 * @brief Parallel tracker with time as the independent variable.
 *
 * @copyright Copyright (c) 200x - 2014, Christof Kraus, Paul Scherrer Institut, Villigen PSI,
 * Switzerland
 * @copyright 2015 - 2016, Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
 * @copyright 2017 - 2020, Christof Metzger-Kraus
 * @copyright 2025 - present, Ryan Ammann, Paul Scherrer Institut, Villigen PSI, Switzerland
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

#ifndef OPALX_ParallelTracker_HH
#define OPALX_ParallelTracker_HH

#include "Algorithms/BunchReference.h"
#include "Algorithms/StepSizeConfig.h"
#include "Algorithms/TrackingSetup.h"

#include "Moments/BeamMoments.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeSolveContext.h"

#include "Algorithms/BeamStop.h"

#include "AbsBeamline/ElementBase.h"
#include "Elements/OpalBeamline.h"

#include <functional>
#include <memory>
#include <span>
#include <vector>

class DataSink;

namespace opalx::spacecharge {
    class SpaceChargeSolver;
}  // namespace opalx::spacecharge

/**
 * @brief Implements the main time-based simulation loop for parallel tracking.
 *
 * @note RUN in the input file constructs a ParallelTracker and runs ParallelTracker::execute().
 */
class ParallelTracker {
public:
    /**
     * @brief Construct the tracker for one RUN.
     * @param setup Everything buildTracking() built for this RUN; setup.bunch must be set.
     */
    explicit ParallelTracker(TrackingSetup setup);
    ParallelTracker(const ParallelTracker&)            = delete;
    ParallelTracker& operator=(const ParallelTracker&) = delete;

    /// @brief Run the main tracking loop until all step-size segments complete.
    void execute();

private:
    // --- Order of the run and of one step ---

    /**
     * @brief Everything before the first step that needs the tracker: the start of the
     *        reference particles, the reference pass with autophase and the initial output.
     * @return False if a restart is already at the end of the schedule.
     */
    bool startTracking();

    /// The step loop and the final output; startTracking() must have returned true.
    void trackSteps();

    /**
     * @brief Start a DT segment: its time step, the container activity, the emitting containers
     *        and the beams already stopped.
     * @param stepsCompletedInSegment Nonzero only where a restart resumes inside this segment;
     *        the container activity of the checkpoint is then kept.
     */
    void startSegment(unsigned long long stepsCompletedInSegment);

    /**
     * @brief The particle part of one step: space charge (PRESTEP or MIDPOINT), the
     *        drift-kick-drift push, emission, external fields, spin, and the global processes
     *        with their deletions.
     *
     * @par Frame of reference
     * Particle R and mesh are in the reference frame for the whole step except inside
     * computeSpaceChargeFields() (the solve frame) and computeExternalFields() (each element's
     * local frame, through a FrameScope).
     */
    void advanceParticles();

    /// @brief Advance the bunch time by one step, then the reference particles and the
    ///        reference-to-lab transforms.
    void advanceReference();

    /**
     * @brief Mark particles outside an element aperture or the BOUNDPDESTROY box and delete them
     *        in one flush. Collective.
     * @return Per container the moments computed for the box, holding what the outputs due now
     *         (@p psDump, @p statDump) read; all empty if particles were deleted.
     */
    std::vector<opalx::moments::BeamMoments> removeLostParticles(bool psDump, bool statDump);

    /**
     * @brief Write a checkpoint every CHECKPOINTFREQ global steps, with the schedule position the
     *        next step resumes at: the next segment if @p segmentComplete, else this one.
     */
    void writeCheckpointIfDue(
            bool segmentComplete, std::size_t segmentIndex,
            unsigned long long stepsCompletedInSegment);

    /// @brief The last output, elements off, the solver counters and the final banner.
    void finishTracking();

    // --- Start of tracking ---

    /// @brief Advance the step-size schedule to the checkpointed position.
    /// @throws OpalException if the input changes the time step there.
    void resumeRestartSchedule();

    // --- Time step ---

    /// @brief Set global bunch dt to dtCurrentTrack_m.
    void selectDT();

    /// @brief Set bunch dt from StepSizeConfig and copy to all container dt views.
    void changeDT();

    // --- Particle push and fields ---

    /// @brief First half of the leapfrog step: push all active containers.
    void timeIntegration1();

    /// @brief Second half: kick then push all active containers.
    void timeIntegration2();

    /// @brief The space-charge solve of this step (SpaceChargeSolver::solveStep()).
    /// @throws OpalException if the run has no space-charge solver.
    void computeSpaceChargeFields();

    /// @brief Apply external fields from elements intersecting each active container.
    void computeExternalFields();

    // --- Output ---

    /// @brief Whether this step prints the per-container step line (STEPINFOFREQ).
    [[nodiscard]] bool stepLogDue() const;

    /// @brief Moment groups the outputs due now (and the step line if @p stepLog) read for @p pc.
    [[nodiscard]] opalx::moments::Request outputMomentsDue(
            const ParticleContainer_t& pc, bool psDump, bool statDump, bool stepLog) const;

    /**
     * @brief Moments of every non-empty container holding what outputMomentsDue() lists; empty
     *        containers get an empty value. Entries of @p computed that already hold it are
     *        reused. Collective.
     */
    [[nodiscard]] std::vector<opalx::moments::BeamMoments> momentsOfContainers(
            bool psDump, bool statDump, bool stepLog,
            std::span<const opalx::moments::BeamMoments> computed) const;

    /**
     * @brief Log per-container stats, trigger dumps according to dump flags and print the rank
     *        load when PRINTRANKDISTRFREQ is due.
     * @param psDump   Forwarded to writePhaseSpace when logging occurs.
     * @param statDump Forwarded to writePhaseSpace when logging occurs.
     * @param computed See writePhaseSpace().
     */
    void dumpStats(
            bool psDump, bool statDump, std::span<const opalx::moments::BeamMoments> computed = {});

    /**
     * @brief Write phase space and/or statistics when flags request it.
     * @param psDump   If true, write phase-space snapshot when applicable.
     * @param statDump If true, write statistics (e.g. SDDS).
     * @param computed Moments already computed at this particle state, reused where they hold
     *        the groups the outputs read.
     */
    void writePhaseSpace(
            bool psDump, bool statDump, std::span<const opalx::moments::BeamMoments> computed = {});

    // --- Members ---

    /// The particle bunch to track. Borrowed; opalx::run::RunObjects owns it.
    PartBunch_t* itsBunch_m;
    /// OPTION values of this RUN; the tracker reads no Options:: globals.
    opalx::run::OptionValues options_m;
    SpaceChargeFieldUpdate spaceChargeFieldUpdate_m =
            SpaceChargeFieldUpdate::MIDPOINT;  ///< Self-field time centering for bunch tracking.
    DataSink* itsDataSink_m;  ///< Borrowed beam statistics and phase-space output sink.
    opalx::spacecharge::SpaceChargeSolver*
            spaceChargeSolver_m;  ///< Borrowed run-lifetime space-charge solver.
    opalx::spacecharge::DirichletPlaneConfig dirichletPlane_m;
    OpalBeamline itsOpalBeamline_m;  ///< Tracked field elements and coordinate transforms.
    BeamStops beamStops_m;           ///< Per container: ZSTOP plane and UNTILMAXSTEPS.
    double sStart_m;                 ///< Path-length start position for the track (m).
    double globalPhaseShift_m;       ///< A fresh run starts this long before t = 0 (s).

    /** Step-size segments: s-stop, dt, and steps per segment. */
    StepSizeConfig stepSizes_m;

    double dtCurrentTrack_m;              ///< Global @f$\Delta t@f$ for the current track segment.
    EmittingSamplers emittingSamplers_m;  ///< Per container, the samplers that emit.
    bool restarting_m;                    ///< Preserve state loaded from a checkpoint at startup.
    unsigned long long restartGlobalStep_m;  ///< Completed global steps restored from checkpoint.
    double restartDt_m;                      ///< Time step stored in the checkpoint.
    StepSizeConfig::ResumePosition restartPosition_m;  ///< Saved schedule segment and offset.
    std::string checkpointPath_m;  ///< File the CHECKPOINTFREQ checkpoints are written to.

    // --- Timers ---
    IpplTimings::TimerRef timeIntegrationTimer1_m;
    IpplTimings::TimerRef timeIntegrationTimer2_m;
    IpplTimings::TimerRef fieldEvaluationTimer_m;

    /// Reference particles; borrows the bunch and the beamline. Constructed after the timers
    /// above, so the timing report lists its ReferencePass timer after them.
    BunchReference references_m;
};

#endif  // OPALX_ParallelTracker_HH
