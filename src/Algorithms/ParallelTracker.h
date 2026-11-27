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

#include "Algorithms/BeamStop.h"
#include "Algorithms/BunchReference.h"
#include "Algorithms/StepSizeConfig.h"
#include "Algorithms/TrackingSetup.h"
#include "Distribution/EmittingSamplers.h"
#include "Elements/OpalBeamline.h"
#include "Moments/BeamMoments.h"
#include "Steppers/NumericalIntegrator.h"
#include "Structure/DataSink.h"

#include <memory>
#include <string>
#include <vector>

namespace opalx::spacecharge {
    class SpaceChargeSolver;
}  // namespace opalx::spacecharge

/**
 * @brief The time-step loop of one RUN: the order of operations before the first step and within
 *        a step, and the hand-over of results between the components that own the work.
 *
 * The tracker owns the step-size schedule, the restart position and the checkpoints. The
 * integrator (NumericalIntegrator) owns the order of the particle part of a step, the bunch
 * (PartBunch) the particles and their push, the beamline (OpalBeamline) the elements, frames and
 * apertures, BunchReference the reference particles, EmittingSamplers the emission,
 * SpaceChargeSolver the self fields, BeamStops the stop rules; DataSink writes the output.
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

    /// @brief Put the reference particles at the start: a fresh run sets their poses and moves
    ///        them to sStart_m, a restart keeps the checkpoint's; then the schedule position.
    void placeReferences();

    /// @brief The track-start banner at bunch time @p time [s], the beam stops and the schedule.
    void printTrackStart(double time);

    /// The step loop and the final output; startTracking() must have returned true.
    void trackSteps();

    /**
     * @brief Start a DT segment: its time step, the container activity, the emitting containers
     *        and the beams already stopped.
     * @param stepsCompletedInSegment Nonzero only where a restart resumes inside this segment;
     *        the container activity of the checkpoint is then kept.
     */
    void startSegment(unsigned long long stepsCompletedInSegment);

    /// @brief The particle part of one step: NumericalIntegrator::advance(), the reset of the
    ///        particle dt to the bunch time step, and the global processes with their deletions.
    void advanceParticles();

    /// @brief Advance the bunch time by one step, then the reference particles and the
    ///        reference-to-lab transforms.
    void advanceReference();

    /**
     * @brief Mark particles outside an element aperture or the BOUNDPDESTROY box and delete them
     *        in one flush. Collective.
     * @return Per container the moments computed for the box, holding what the outputs @p due
     *         read; all empty if particles were deleted.
     */
    std::vector<opalx::moments::BeamMoments> removeLostParticles(const DataSink::StepOutputs& due);

    /**
     * @brief Write a checkpoint every CHECKPOINTFREQ global steps, with the schedule position the
     *        next step resumes at: the next segment if @p segmentComplete, else this one.
     */
    void writeCheckpointIfDue(
            bool segmentComplete, std::size_t segmentIndex,
            unsigned long long stepsCompletedInSegment);

    /// @brief The last output, elements off, the solver counters and the final banner.
    void finishTracking();

    /// @brief The time step for the DT of the current segment, dtCurrentTrack_m: the smallest
    ///        emission time step while a source still emits.
    [[nodiscard]] double selectDT() const;

    /// @brief Set the time step of the bunch and of all its containers to selectDT().
    void changeDT();

    // --- Output ---

    /// The external fields at the reference particles for the dumps (OpalBeamline).
    [[nodiscard]] std::function<DataSink::ReferenceFields()> referenceFields();

    // --- Members ---

    /// The particle bunch to track. Borrowed; opalx::run::RunObjects owns it.
    PartBunch_t* itsBunch_m;
    /// OPTION values of this RUN.
    opalx::run::OptionValues options_m;
    DataSink* itsDataSink_m;  ///< Borrowed beam statistics and phase-space output sink.
    opalx::spacecharge::SpaceChargeSolver*
            spaceChargeSolver_m;     ///< Borrowed run-lifetime space-charge solver; never null.
    OpalBeamline itsOpalBeamline_m;  ///< Tracked field elements and coordinate transforms.
    BeamStops beamStops_m;           ///< Per container: ZSTOP plane and UNTILMAXSTEPS.
    EmittingSamplers emittingSamplers_m;  ///< Per container, the samplers that emit.

    StepSizeConfig stepSizes_m;  ///< Step-size segments: s-stop, dt and steps per segment.
    double dtCurrentTrack_m;     ///< DT of the current segment (s).
    double sStart_m;             ///< Path-length start position for the track (m).
    double globalPhaseShift_m;   ///< A fresh run starts this long before t = 0 (s).
    bool restarting_m;           ///< Preserve state loaded from a checkpoint at startup.
    unsigned long long restartGlobalStep_m;  ///< Completed global steps restored from checkpoint.
    double restartDt_m;                      ///< Time step stored in the checkpoint.
    StepSizeConfig::ResumePosition restartPosition_m;  ///< Saved schedule segment and offset.
    std::string checkpointPath_m;  ///< File the CHECKPOINTFREQ checkpoints are written to.
    std::vector<opalx::run::EMFieldDumpDescription> emFieldDumps_m;  ///< DUMPEMFIELDS requests.

    /// Timers of execute() and its parts. They do not overlap and, with the timers of main() and
    /// TrackRun::execute(), add up to mainTimer. Timers started inside them (for example
    /// computeMoments and Write Stat) are detail and not part of that sum.
    IpplTimings::TimerRef trackStartTimer_m;
    IpplTimings::TimerRef stepTimeStepTimer_m;
    /// The particle part of each step, scheme TrackingSetup::timeIntegrator. It borrows the
    /// beamline and the samplers above and creates its timers here, between those around it.
    std::unique_ptr<NumericalIntegrator> integrator_m;
    IpplTimings::TimerRef stepProcessesTimer_m;
    IpplTimings::TimerRef stepReferenceTimer_m;
    IpplTimings::TimerRef stepRemoveLostTimer_m;
    IpplTimings::TimerRef stepOutputTimer_m;
    IpplTimings::TimerRef stepStopRulesTimer_m;
    IpplTimings::TimerRef stepCheckpointTimer_m;
    IpplTimings::TimerRef trackFinishTimer_m;

    /// Reference particles; borrows the bunch, the beamline and integrator_m, which is declared
    /// before it. Constructed after the timers above, so the timing report lists its
    /// ReferencePass timer after them.
    BunchReference references_m;
};

#endif  // OPALX_ParallelTracker_HH
