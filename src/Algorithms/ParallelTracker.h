/**
 * @file ParallelTracker.h
 * @brief Visitor-based parallel tracker with time as the independent variable.
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

#include <optional>
#include "Algorithms/StepSizeConfig.h"
#include "Algorithms/TrackingSetup.h"
#include "Steppers/BorisPusher.h"
#include "Steppers/SpinTBMTPusher.h"
#include "Structure/DataSink.h"

#include "Moments/BeamMoments.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeSolveContext.h"

#include "BasicActions/Option.h"
#include "Utilities/Options.h"

#include "Physics/Physics.h"

#include "Algorithms/BeamStop.h"
#include "Algorithms/GeneralElementMapping.h"

#include "AbsBeamline/Box.h"
#include "AbsBeamline/Collimator.h"
#include "AbsBeamline/ConstantEFieldCavity.h"
#include "AbsBeamline/ConstantFocusing.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/FieldmapElement.h"
#include "AbsBeamline/Laser.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/MultipoleT.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/VariableRFCavity.h"
#include "Beamlines/Beamline.h"
#include "Distribution/SamplingBase.hpp"
#include "Elements/OpalBeamline.h"

#include <functional>
#include <list>
#include <memory>
#include <set>
#include <tuple>
#include <vector>

namespace opalx::spacecharge {
    class SpaceChargeSolver;
}  // namespace opalx::spacecharge

/**
 * @brief Implements the main time-based simulation loop for parallel tracking.
 *
 * @note TRACK and RUN in the input file construct a ParallelTracker and run
 *       ParallelTracker::execute().
 */
class ParallelTracker {
private:
    /// The particle bunch to track. Borrowed; opalx::run::RunObjects owns it.
    PartBunch_t* itsBunch_m;
    /// OPTION values of this RUN; the tracker reads no Options:: globals.
    opalx::run::OptionValues options_m;
    /// Per-container reference momentum [beta*gamma] set on BEAM; empty uses the bunch mean.
    std::vector<std::optional<Vector_t<double, 3>>> beamReferenceMomenta_m;
    /**
     * @brief Everything before the first step that needs the tracker: the start of the
     *        reference particles, the reference pass with autophase and the initial output.
     * @return False if a restart is already at the end of the schedule.
     */
    bool startTracking();
    /// Bunch time, time step and step number of the current step, passed to element apply().
    StepTime stepTime() const {
        return {itsBunch_m->getT(), itsBunch_m->getdT(),
                static_cast<long long>(itsBunch_m->getGlobalTrackStep())};
    }
    /// The step loop and the final output; startTracking() must have returned true.
    void trackSteps();
    SpaceChargeFieldUpdate spaceChargeFieldUpdate_m =
            SpaceChargeFieldUpdate::MIDPOINT;  ///< Self-field time centering for bunch tracking.
    DataSink* itsDataSink_m;  ///< Borrowed beam statistics and phase-space output sink.
    opalx::spacecharge::SpaceChargeSolver*
            spaceChargeSolver_m;  ///< Borrowed run-lifetime space-charge solver.
    opalx::spacecharge::DirichletPlaneConfig dirichletPlane_m;
    std::vector<std::uint8_t> spaceChargeContainerActivity_m;
    OpalBeamline itsOpalBeamline_m;     ///< Cloned field elements and coordinate transforms.
    bool referenceHitMaterial_m;        ///< A reference particle hit an element: the run ends.
    std::vector<BeamStop> beamStops_m;  ///< Per container: ZSTOP plane and UNTILMAXSTEPS.
    double sStart_m;                    ///< Path-length start position for the track (m).
    double globalPhaseShift_m;          ///< A fresh run starts this long before t = 0 (s).

    /** Step-size segments: s-stop, dt, and steps per segment. */
    StepSizeConfig stepSizes_m;

    double dtCurrentTrack_m;  ///< Global @f$\Delta t@f$ for the current track segment.
    std::vector<std::vector<std::shared_ptr<SamplingBase>>>
            emittingSamplers_m;  ///< Per-container emitters.
    bool restarting_m;           ///< Preserve state loaded from a checkpoint at startup.
    unsigned long long restartGlobalStep_m;  ///< Completed global steps restored from checkpoint.
    double restartDt_m;                      ///< Time step stored in the checkpoint.
    StepSizeConfig::ResumePosition restartPosition_m;  ///< Saved schedule segment and offset.
    CavityPhases cavityPhases_m;  ///< Phases autophase set (and a restart restored).

    // --- Timers ---
    IpplTimings::TimerRef timeIntegrationTimer1_m;
    IpplTimings::TimerRef timeIntegrationTimer2_m;
    IpplTimings::TimerRef fieldEvaluationTimer_m;
    IpplTimings::TimerRef referencePassTimer_m;

public:
    /**
     * @brief Construct the tracker for one RUN.
     * @param setup Everything buildTracking() built for this RUN; setup.bunch must be set.
     */
    explicit ParallelTracker(TrackingSetup setup);

    /// @brief Destructor; releases tracker resources.
    virtual ~ParallelTracker();

    /// @brief Run the main tracking loop until all step-size segments complete.
    void execute();

    /**
     * @brief Boris half-kick using E, B and per-particle dt on one container.
     * @param pusher Boris pusher instance.
     * @param pc     Non-null particle container.
     */
    void kickParticles(const BorisPusher& pusher, PartBunch_t::ParticleContainer_t& pc);

    /**
     * @brief Boris position push (unitless positions) on one container.
     * @param pusher Boris pusher instance.
     * @param pc     Non-null particle container.
     */
    void pushParticles(const BorisPusher& pusher, PartBunch_t::ParticleContainer_t& pc);

    /// @brief First half of the leapfrog step: push all active containers.
    void timeIntegration1(BorisPusher& pusher);

    /// @brief Second half: kick then push all active containers.
    void timeIntegration2(BorisPusher& pusher);

    /// @brief First half of the Thomas-BMT spin rotation for all active containers that store
    /// Pol, with the momentum before the kick. timeIntegration2 does the second half with the
    /// momentum after the kick, so the spin update is second order in dt. Call after the
    /// external and space-charge fields at the half-step position have been accumulated.
    void evolveSpinTBMT();

    /// @brief Rotate Pol of every particle in @p pc by half its own dt (Thomas-BMT), with the
    /// current P, E and B. P, E, B and Pol are all along the bunch-frame axes.
    void evolveSpinHalfStep(PartBunch_t::ParticleContainer_t& pc);

    /** @brief Build the current tracker-frame context and dispatch the configured space-charge
     * solve. */
    void computeSpaceChargeFields();

    /// @brief Apply external fields from elements intersecting each active container.
    void computeExternalFields();

    /**
     * @brief Call @p func for each selected element intersecting each active container, with
     * the container in that element's local frame.
     *
     * The container goes from one selected element's frame directly to the next and returns to
     * its reference frame after its last element, so no frame outlives the call. Unselected
     * elements cost no transform.
     *
     * @param selects Whether @p func runs for an element.
     */
    void forEachElementInBunchFrame(
            const std::function<bool(const ElementBase&)>& selects,
            const std::function<
                    void(const std::shared_ptr<ElementBase>&,
                         const std::shared_ptr<ParticleContainer_t>&)>& func);

    /// Elements near container @p ci in the order of the element grid; empty for an inactive
    /// container or one without its own or its root's reference particle.
    std::vector<std::shared_ptr<ElementBase>> candidateElements(size_t ci);

    /// Beamline (emission source) frame for the particles of @p pc.
    ParticleFrame beamlineFrame(const ParticleContainer_t& pc) const;

    /// Local frame of @p element, misalignment included, for the particles of @p pc.
    ParticleFrame elementFrame(
            const std::shared_ptr<ElementBase>& element, const ParticleContainer_t& pc) const;

    /// @brief Mark particles outside the transverse aperture of each nearby element.
    /// @return global number of newly marked particles (allreduced) - collective call.
    size_t applyElementApertures();

    /// @brief Emit macroparticles from configured samplers per container.
    /// @param t  Bunch time (s).
    /// @param dt Global time step (s).
    void emitFromEmissionSources(double t, double dt);

    /// @brief Mark particles moving backward behind an active source/cathode plane.
    size_t markBackwardParticlesAtSourcePlane();

    /// @brief Apply global processes and return the global number of particles marked invalid.
    size_t applyGlobalProcesses(double dt);

    /// @brief Zero E and B on all active particle containers.
    void resetFields();

    /// @brief Set bunch dt from StepSizeConfig and copy to all container dt views.
    void changeDT();

private:
    struct SpaceChargeEmissionProgress {
        bool active     = false;
        double fraction = 1.0;
    };

    /// @brief Build stable native particle and field identities once.
    void initializeSpaceChargeContainerActivity();

    /// @brief The requested solve frame, from the primary container's moments. Collective
    /// unless @p request is the reference frame.
    [[nodiscard]] ParticleFrame makeSpaceChargeFrame(opalx::spacecharge::SolveFrame request) const;
    [[nodiscard]] SpaceChargeEmissionProgress spaceChargeEmissionProgress() const;

    /// @brief Update reference trajectories and lab/reference coordinate transforms.
    void updateReference(const BorisPusher& pusher);

    /// @brief Advance reference positions/momenta through the beamline for one step.
    void updateReferenceParticles(const BorisPusher& pusher);

    /// @brief Refresh each container's reference-to-lab transform from current state.
    void updateRefToLabCSTrafo();

    /**
     * @brief Copy each decay daughter container's root reference particle, frame and path
     *   length into it, without moving particles.
     * @note Used after the reference particles are set up; during tracking the daughter
     *   containers follow their root in updateRefToLabCSTrafo().
     */
    void copyFollowerReferences();

    /**
     * @brief Write phase space and/or statistics when flags request it.
     * @param step     Step index (reserved for diagnostics; may be unused in body).
     * @param psDump   If true, write phase-space snapshot when applicable.
     * @param statDump If true, write statistics (e.g. SDDS).
     * @param computed Moments already computed at this particle state, reused where they hold
     *        the groups the outputs read.
     */
    void writePhaseSpace(
            const long long step, bool psDump, bool statDump,
            std::span<const opalx::moments::BeamMoments> computed = {});

    /**
     * @brief Log per-container stats and trigger dumps according to dump flags.
     * @param step     Current integration step index.
     * @param psDump   Forwarded to writePhaseSpace when logging occurs.
     * @param statDump Forwarded to writePhaseSpace when logging occurs.
     * @param computed See writePhaseSpace().
     */
    void dumpStats(
            long long step, bool psDump, bool statDump,
            std::span<const opalx::moments::BeamMoments> computed = {});

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

    /// @brief Set global bunch dt to dtCurrentTrack_m.
    void selectDT();

    /// @brief Mark a beam done once every particle has passed its stop plane (collective call).
    void updateBeamStops();

    /// @brief Delete particles marked invalid by the central per-container mask.
    size_t deleteInvalidParticles(bool activeOnly, Inform& m, const std::string& reason);

    /// @brief Force-activate containers whose emitting samplers have not yet finished.
    void activateEmittingContainers(double t);

    /**
     * @brief Union of per-container spatial bounds over MPI.
     * @param[out] rmin Corner of the axis-aligned bounding box (min).
     * @param[out] rmax Corner of the axis-aligned bounding box (max).
     */
    void computeInitialBounds(Vector_t<double, 3>& rmin, Vector_t<double, 3>& rmax);

    /**
     * @brief Log reference state for each container at track start.
     * @param m Inform stream for log output.
     */
    void printInitialContainerRefs(Inform& m) const;

    /// @brief Integrate references in time until path length reaches sStart_m.
    void findStartPositions(const BorisPusher& pusher);

    /// @brief Autophase TRAVELINGWAVE and RFCAVITY elements along the reference orbit.
    void autophaseCavities(const BorisPusher& pusher);
};

#endif  // OPALX_ParallelTracker_HH
