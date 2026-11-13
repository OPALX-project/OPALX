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
#include "Algorithms/BorisStepControl.h"
#include "Algorithms/ClosedOrbitInitialState.h"
#include "Algorithms/DeviceExternalField.h"
#include "Algorithms/SpectralTunes.h"
#include "Algorithms/StepSizeConfig.h"
#include "Algorithms/Tracker.h"
#include "Algorithms/TrackingSetup.h"
#include "Steppers/BorisPusher.h"
#include "Steppers/SpinTBMTPusher.h"
#include "Structure/DataSink.h"

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
#include "AbsBeamline/CyclotronSector.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/FieldmapElement.h"
#include "AbsBeamline/Laser.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/MultipoleT.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/ScalingFFAMagnet.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/VariableRFCavity.h"
#include "Algorithms/CyclotronProbeRecorder.h"
#include "Beamlines/Beamline.h"
#include "Distribution/SamplingBase.hpp"
#include "Elements/OpalBeamline.h"

#include <functional>
#include <list>
#include <memory>
#include <set>
#include <tuple>
#include <vector>

class PluginElement;

namespace opalx::spacecharge {
    class SpaceChargeSolver;
}  // namespace opalx::spacecharge

/**
 * @brief Implements the main time-based simulation loop for parallel tracking.
 *
 * @note TRACK and RUN in the input file construct a ParallelTracker and run
 *       ParallelTracker::execute().
 */
class ParallelTracker : public Tracker {
public:
    // Localized TURNS (TrackingSetup::requestedTurns): stop at the Nth forward reference return
    // (single container, supported magnetic/RF elements, no ongoing emission). The complete bunch
    // advances to that reference event time; individual particles need not close.
    // Turn counting is independent of the ordinary device Boris/PIC integration.
    // Compatible bare analytic rings share boundary-controlled substeps so all
    // particles reach each midpoint together. COF/reference remain host-side.
    // This changes rounding and edge truncation errors from the legacy kernels,
    // while retaining the energy-conserving magnetic Boris rotation. Boundary
    // refinement excludes spin. Immutable field descriptors are uploaded once; particle data
    // stay in device views. Supports drifts, passive elements, analytic bends,
    // dipoles, quadrupoles and ideal VariableRFCavity elements. Shared boundary
    // retries are restricted to bare (NONE solver) tracking; PIC retains the
    // ordinary step sequence independently of any diagnostic probe selection.
    // Checkpoint/restart and solver configurations that cannot repeat midpoint
    // trials retain fixed Boris/PIC steps with spatial field selection.
    void visitCyclotronSector(const CyclotronSector& sector) override {
        itsOpalBeamline_m.visit(sector, *this, *itsBunch_m);
    }

    /// Register global diagnostics separately from field/placement elements.
    void visitProbe(const Probe& probe) override {
        probeConfigs_m.push_back(CyclotronProbeRecorder::configuration(probe));
    }

private:
    std::vector<CyclotronProbeRecorder::Config> probeConfigs_m;
    bool bareTracking_m =
            false;  ///< NONE backend; diagnostic eligibility is independent of retries.
    bool allowBoundaryControl_m =
            true;  ///< Restrict retries to bare, undiagnosed Cartesian solves.
    std::optional<ClosedOrbitInitialState> initialOrbit_m;
    /// Per-container reference momentum [beta*gamma] set on BEAM; empty uses the bunch mean.
    std::vector<std::optional<Vector_t<double, 3>>> beamReferenceMomenta_m;
    std::vector<double> tuneInitial_m;
    SpectralTunes::Settings tuneSettings_m;
    bool hasCyclotronGaps();
    std::optional<device_external::CyclotronLattice> cyclotronBunchLattice_m;
    /// Fixed-size first-failure payload, allocated once when device RF transport is prepared.
    Kokkos::View<rf_event::ParticleFailure> cyclotronFailure_m;
    /// Validate fixed bunch requirements once; leave physical stepping order unchanged.
    void prepareCyclotronTransport();
    /// Upload immutable magnetic/gap descriptors only for the device path.
    void prepareCyclotronBunchLattice();
    bool usesHostCyclotronTransport() const;
    void advanceCyclotronParticle();
    /** Host-only event integration for the reference and a bare single proton.
     * R/P are in the laboratory frame, t/dt in seconds. Splits Boris steps at
     * directed gap-plane roots and applies full kicks in time order. Field
     * queries use spatial support, not a nominal closed orbit.
     * Assumes distinct gap planes and a single unpolarized median-plane proton.
     * The caller must choose dt small enough to bracket each crossing without
     * an intervening recrossing of the same plane. Particle E/B diagnostic views
     * are not populated by this event path; only positions/momenta are advanced.
     * @return Elapsed time [s], equal to dt unless EKINSTOP is met at a gap.
     * On target completion the outgoing state is immediately after that complete
     * kick; the remainder magnetic drift is deliberately omitted. With report=true
     * the reference marks energyTargetReached_m and emits the terminal diagnostic.
     * Physical particles (report=false) finish the entire supplied duration even
     * if their own energy exceeds the reference target at an earlier gap.
     * terminalEndpoint applies the shared roundoff-inclusive final endpoint policy
     * to physical particles after that reference stop, never to reference trials.
     * frameSteps counts accepted frame updates including this step; it budgets
     * accumulated coordinate arithmetic in the terminal-only rf_event allowance.
     */
    double advanceCyclotronGaps(
            Vector_t<double, 3>& r, Vector_t<double, 3>& p, double t, double dt, double mass,
            bool report, bool terminalEndpoint = false, unsigned long long frameSteps = 1);
    double kineticEnergyStop_m = 0;      ///< Optional reference kinetic-energy target [eV].
    bool energyTargetReached_m = false;  ///< Latched only by the reference's complete kick.
    /// Energy-mode reference is precomputed so the final bunch clock uses its substep.
    bool pendingEnergyReference_m = false;
    Vector_t<double, 3> pendingReferenceR_m, pendingReferenceP_m;
    unsigned long long requestedTurns_m = 0;
    device_external::Lattice
            deviceRingFields_m;  ///< Geometry for device field selection and trial support checks.
    bool spatialRing_m =
            false;  ///< Select analytic ring fields by physical position, independent of retries.
    bool boundaryControlled_m = false;
    /// Internal subset experiment; one selects every particle (ordinary default).
    unsigned long long boundaryControlStride_m = 1;
    double boundaryStepDt_m                    = 0;  ///< Accepted/trial collective substep cap [s].
    unsigned long long boundaryTrials_m = 0, boundaryRejected_m = 0;
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
    double ringPeriod_m;                ///< One-turn path length for RING, or zero for LINE.

    /** Step-size segments: s-stop, dt, and steps per segment. */
    StepSizeConfig stepSizes_m;

    double dtCurrentTrack_m;      ///< Global @f$\Delta t@f$ for the current track segment.
    double terminalStepDt_m = 0;  ///< Positive final-turn cap [s]; zero means no cap.
    std::vector<std::vector<std::shared_ptr<SamplingBase>>>
            emittingSamplers_m;  ///< Per-container emitters.
    bool restarting_m;           ///< Preserve state loaded from a checkpoint at startup.
    unsigned long long restartGlobalStep_m;  ///< Completed global steps restored from checkpoint.
    double restartDt_m;                      ///< Time step stored in the checkpoint.
    StepSizeConfig::ResumePosition restartPosition_m;  ///< Saved schedule segment and offset.

    // --- Timers ---
    IpplTimings::TimerRef timeIntegrationTimer1_m;
    IpplTimings::TimerRef timeIntegrationTimer2_m;
    IpplTimings::TimerRef fieldEvaluationTimer_m;
    IpplTimings::TimerRef WakeFieldTimer_m;
    IpplTimings::TimerRef PluginElemTimer_m;
    IpplTimings::TimerRef referencePassTimer_m;

public:
    /**
     * @brief Construct tracker with beamline only (no bunch attached here).
     * @param bl       Beamline definition.
     * @param revBeam  If true, reversed beam direction (s = C to 0); OPAL-T parallel
     *                 tracking is forward-only; revTrack is not used.
     */
    explicit ParallelTracker(const Beamline& bl, bool revBeam);

    /**
     * @brief Construct the tracker for one RUN.
     * @param bl    Beamline definition.
     * @param setup Everything TrackRun built for this RUN; setup.bunch must be set.
     */
    ParallelTracker(const Beamline& bl, TrackingSetup setup);

    /// @brief Destructor; releases tracker resources.
    virtual ~ParallelTracker();

    /// @brief Visit the full beamline (iterates elements into OpalBeamline). Overrides
    /// DefaultVisitor.
    void visitBeamline(const Beamline&) override;

    /// @brief Visit a generic element using the base tracker behavior.
    void visitElementBase(const ElementBase&) override;

    /// @brief Apply the algorithm to a constant E-field cavity.
    void visitConstantEFieldCavity(const ConstantEFieldCavity&) override;

    /// @brief Apply the algorithm to a constant linear focusing element.
    void visitConstantFocusing(const ConstantFocusing&) override;
    /// @brief Apply the algorithm to a box absorber.
    void visitBox(const Box&) override;

    /// @brief Apply the algorithm to a collimator.
    void visitCollimator(const Collimator&) override;

    /// @brief Apply the algorithm to a drift.
    void visitDrift(const Drift&) override;

    /// @brief Apply the algorithm to a field-map-driven element.
    void visitFieldmapElement(const FieldmapElement&) override;

    /// @brief Reject laser tracking until dedicated laser tracking is implemented.
    void visitLaser(const Laser&) override;

    /// @brief Apply the algorithm to a monitor.
    void visitMonitor(const Monitor&) override;

    /// @brief Apply the algorithm to a marker.
    void visitMarker(const Marker&) override;

    /// @brief Apply the algorithm to a multipole.
    void visitMultipole(const Multipole&) override;

    /// @brief Apply the algorithm to a multipole (templated type).
    void visitMultipoleT(const MultipoleT&) override;

    /// @brief Apply the algorithm to a rectangular bend.
    void visitRBend(const RBend&) override;

    /// @brief Apply the algorithm to an RF cavity.
    void visitRFCavity(const RFCavity&) override;

    /// @brief Register and initialise an analytic time-dependent RF cavity.
    void visitVariableRFCavity(const VariableRFCavity&) override;

    /// @brief Apply the algorithm to a sector bend.
    void visitSBend(const SBend&) override;

    /// @brief Null op, as everywhere. Why does this exist?
    void visitScalingFFAMagnet(const ScalingFFAMagnet&) override;

    /// @brief Apply the algorithm to a traveling wave cavity.
    void visitTravelingWave(const TravelingWave&) override;

    /// @brief Apply the algorithm to a solenoid.
    void visitSolenoid(const Solenoid&) override;

    void visitEnge(const endfieldmodel::Enge&) override;

    void visitAsymmetricEnge(const endfieldmodel::AsymmetricEnge&) override;

    void visitTanh(const endfieldmodel::Tanh&) override;

    /// @brief Run the main tracking loop until all step-size segments complete.
    void execute() override;

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

    /// @brief Call func for each element intersecting each active container, with the
    ///        container transformed into the element-local frame.
    /// @param func Called per (element, container) pair in the element-local frame.
    void forEachElementInBunchFrame(
            const std::function<
                    void(const std::shared_ptr<ElementBase>&,
                         const std::shared_ptr<ParticleContainer_t>&)>& func);

    /**
     * @brief Internal boundary-controlled device tracking helpers.
     *
     * These entry points remain implementation details, but are public because
     * CUDA extended lambdas require a public enclosing member function. Keep
     * their state private and revisit this interface when the kernels are moved
     * to namespace-scope functors.
     */
    /** Complete the first drift and field gathering at a common physical midpoint.
     * A register-only endpoint trial decides collective subdivision before any
     * momentum kick, reference update, emission or loss is committed. Rejected
     * first drifts are reversed on the solver's current particle ownership.
     */
    void prepareBoundaryStep(BorisPusher&, boris_step::Control&);
    bool boundaryCrossed(double dt);
    void reverseTrialDrift(double dt);

    /// Internal candidate-selection policy; the public two-argument API is unchanged.
    void forEachElementInBunchFrame(
            const std::function<
                    void(const std::shared_ptr<ElementBase>&,
                         const std::shared_ptr<ParticleContainer_t>&)>& func,
            bool spatialCandidates);

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

    /// @brief Reset per-particle dt views to the current global bunch dt.
    void setTime();

private:
    struct SpaceChargeEmissionProgress {
        bool active     = false;
        double fraction = 1.0;
    };

    /// @brief Build stable native particle and field identities once.
    void initializeSpaceChargeContainerActivity();

    [[nodiscard]] opalx::spacecharge::CoordinateFrameTransforms makeSpaceChargeFrameTransforms()
            const;
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
     */
    void writePhaseSpace(const long long step, bool psDump, bool statDump);

    /**
     * @brief Log per-container stats and trigger dumps according to dump flags.
     * @param step     Current integration step index.
     * @param psDump   Forwarded to writePhaseSpace when logging occurs.
     * @param statDump Forwarded to writePhaseSpace when logging occurs.
     */
    void dumpStats(long long step, bool psDump, bool statDump);

    /// @brief Accept beamline visitor, prepare sections, compute and save 3D lattice.
    void prepareSections();

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

    /**
     * @brief Set stored RF phase on the named cavity or traveling-wave element.
     * @param elName  Element name to match.
     * @param maxPhi  RF phase to apply (rad).
     */
    void updateRFElement(std::string elName, double maxPhi);

    /// @brief Print RF phases (debug/diagnostic hook).
    void printRFPhases();

    /// @brief Persist cavity phases to the data sink.
    void saveCavityPhases();

    /// @brief Restore cavity phases from a prior track or restart.
    void restoreCavityPhases();
};

inline void ParallelTracker::visitConstantEFieldCavity(const ConstantEFieldCavity& cav) {
    itsOpalBeamline_m.visit(cav, *this, *itsBunch_m);
}

inline void ParallelTracker::visitConstantFocusing(const ConstantFocusing& focusing) {
    itsOpalBeamline_m.visit(focusing, *this, *itsBunch_m);
}

inline void ParallelTracker::visitBox(const Box& box) {
    itsOpalBeamline_m.visit(box, *this, *itsBunch_m);
}

inline void ParallelTracker::visitCollimator(const Collimator& coll) {
    itsOpalBeamline_m.visit(coll, *this, *itsBunch_m);
}

inline void ParallelTracker::visitDrift(const Drift& drift) {
    itsOpalBeamline_m.visit(drift, *this, *itsBunch_m);
}

inline void ParallelTracker::visitFieldmapElement(const FieldmapElement& fm) {
    itsOpalBeamline_m.visit(fm, *this, *itsBunch_m);
}

inline void ParallelTracker::visitMonitor(const Monitor& monitor) {
    itsOpalBeamline_m.visit(monitor, *this, *itsBunch_m);
}

inline void ParallelTracker::visitMarker(const Marker& marker) {
    itsOpalBeamline_m.visit(marker, *this, *itsBunch_m);
}

inline void ParallelTracker::visitMultipole(const Multipole& mult) {
    itsOpalBeamline_m.visit(mult, *this, *itsBunch_m);
}

inline void ParallelTracker::visitMultipoleT(const MultipoleT& mult) {
    itsOpalBeamline_m.visit(mult, *this, *itsBunch_m);
}

inline void ParallelTracker::visitRBend(const RBend& bend) {
    itsOpalBeamline_m.visit(bend, *this, *itsBunch_m);
}

inline void ParallelTracker::visitRFCavity(const RFCavity& as) {
    itsOpalBeamline_m.visit(as, *this, *itsBunch_m);
}

inline void ParallelTracker::visitVariableRFCavity(const VariableRFCavity& cavity) {
    itsOpalBeamline_m.visit(cavity, *this, *itsBunch_m);
}

inline void ParallelTracker::visitSBend(const SBend& bend) {
    itsOpalBeamline_m.visit(bend, *this, *itsBunch_m);
}

inline void ParallelTracker::visitScalingFFAMagnet(const ScalingFFAMagnet& ffa) {
    itsOpalBeamline_m.visit(ffa, *this, *itsBunch_m);
}

inline void ParallelTracker::visitTravelingWave(const TravelingWave& tw) {
    itsOpalBeamline_m.visit(tw, *this, *itsBunch_m);
}

inline void ParallelTracker::visitSolenoid(const Solenoid& so) {
    itsOpalBeamline_m.visit(so, *this, *itsBunch_m);
}

inline void ParallelTracker::visitEnge(const endfieldmodel::Enge&) {}

inline void ParallelTracker::visitAsymmetricEnge(const endfieldmodel::AsymmetricEnge&) {}

inline void ParallelTracker::visitTanh(const endfieldmodel::Tanh&) {}

#endif  // OPALX_ParallelTracker_HH
