/**
 * @file ParallelTracker.h
 * @brief Visitor-based parallel tracker with time as the independent variable.
 *
 * @copyright Copyright (c) 200x - 2014, Christof Kraus, Paul Scherrer Institut, Villigen PSI, Switzerland
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

#include "Algorithms/StepSizeConfig.h"
#include "Algorithms/Tracker.h"
#include "Steppers/BorisPusher.h"
#include "Structure/DataSink.h"

#include "BasicActions/Option.h"
#include "Utilities/Options.h"

#include "Physics/Physics.h"

#include "Algorithms/IndexMap.h"
#include "Algorithms/OrbitThreader.h"

#include "AbsBeamline/ConstantEFieldCavity.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/MultipoleT.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/ScalingFFAMagnet.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/TravelingWave.h"
#include "Beamlines/Beamline.h"
#include "Distribution/SamplingBase.hpp"
#include "Elements/OpalBeamline.h"

#include <list>
#include <memory>
#include <tuple>
#include <vector>

class ParticleMatterInteractionHandler;
class PluginElement;

/**
 * @brief Implements the main time-based simulation loop for parallel tracking.
 *
 * @note TRACK and RUN in the input file construct a ParallelTracker and run
 *       ParallelTracker::execute().
 */
class ParallelTracker : public Tracker {
private:

    std::shared_ptr<DataSink> itsDataSink_m; ///< Beam statistics and phase-space output.
    OpalBeamline itsOpalBeamline_m;          ///< Cloned field elements and coordinate transforms.
    bool globalEOL_m;                       ///< End-of-line flag (e.g. orbit threader out of bounds).
    double zstart_m;                        ///< Path-length start position for the track (m).

    /** Step-size segments: z-stop, dt, and steps per segment. */
    StepSizeConfig stepSizes_m;

    double dtCurrentTrack_m;                ///< Global @f$\Delta t@f$ for the current track segment.
    unsigned long long repartFreq_m;       ///< Space-charge repartition period (steps); off on one rank.
    std::vector<std::vector<std::shared_ptr<SamplingBase>>> emittingSamplers_m; ///< Per-container emitters.

    // --- Timers ---
    IpplTimings::TimerRef timeIntegrationTimer1_m;
    IpplTimings::TimerRef timeIntegrationTimer2_m;
    IpplTimings::TimerRef fieldEvaluationTimer_m;
    IpplTimings::TimerRef WakeFieldTimer_m;
    IpplTimings::TimerRef PluginElemTimer_m;
    IpplTimings::TimerRef BinRepartTimer_m;
    IpplTimings::TimerRef OrbThreader_m;

public:
    /**
     * @brief Construct tracker with beamline only (no bunch attached here).
     * @param bl       Beamline definition.
     * @param revBeam  If true, reversed beam direction (s = C to 0); OPAL-T parallel
     *                 tracking is forward-only; revTrack is not used.
     */
    explicit ParallelTracker(const Beamline& bl, bool revBeam);

    /**
     * @brief Construct tracker with bunch, output sink, and step-size schedule.
     * @param bl                Beamline definition.
     * @param bunch             Particle bunch (multi-container).
     * @param ds                Data sink for statistics and dumps.
     * @param revBeam           Reversed beam flag (see single-argument constructor).
     * @param maxSTEPS          Max integration steps per z-segment (parallel to zstop/dt).
     * @param zstart            Starting path length (m).
     * @param zstop             Stop path length per segment (m).
     * @param dt                Time step per segment (s).
     * @param emittingSamplers  Optional per-container samplers for emitParticles(t, dt).
     */
    explicit ParallelTracker(const Beamline& bl, std::shared_ptr<PartBunch_t> bunch,
        const std::shared_ptr<DataSink>& ds, bool revBeam,
        const std::vector<unsigned long long>& maxSTEPS, 
        double zstart, const std::vector<double>& zstop, 
        const std::vector<double>& dt,
        const std::vector<std::vector<std::shared_ptr<SamplingBase>>>& emittingSamplers = {});

    /// @brief Destructor; releases tracker resources.
    virtual ~ParallelTracker();

    /// @brief Visit the full beamline (iterates elements into OpalBeamline). Overrides DefaultVisitor.
    virtual void visitBeamline(const Beamline&);

    /// @brief Visit a generic component; rejects LASER (not implemented).
    virtual void visitComponent(const Component&);

    /// @brief Apply the algorithm to a constant E-field cavity.
    virtual void visitConstantEFieldCavity(const ConstantEFieldCavity&);

    /// @brief Apply the algorithm to a drift.
    virtual void visitDrift(const Drift&);

    /// @brief Apply the algorithm to a marker.
    virtual void visitMarker(const Marker&);

    /// @brief Apply the algorithm to a multipole.
    virtual void visitMultipole(const Multipole&);

    /// @brief Apply the algorithm to a multipole (templated type).
    virtual void visitMultipoleT(const MultipoleT&);

    /// @brief Apply the algorithm to an RF cavity.
    virtual void visitRFCavity(const RFCavity&);

    /// @brief Apply the algorithm to a traveling wave cavity.
    virtual void visitTravelingWave(const TravelingWave&);

    /// @brief Apply the algorithm to a solenoid.
    virtual void visitSolenoid(const Solenoid&);

    /// @brief Run the main tracking loop until all step-size segments complete.
    virtual void execute();

    /**
     * @brief Boris half-kick using E, B and per-particle dt on one container.
     * @param pusher Boris pusher instance.
     * @param pc     Particle container.
     */
    void kickParticles(const BorisPusher& pusher,
                       const std::shared_ptr<PartBunch_t::ParticleContainer_t>& pc);
                       
    /**
     * @brief Boris position push (unitless positions) on one container.
     * @param pusher Boris pusher instance.
     * @param pc     Particle container.
     */
    void pushParticles(const BorisPusher& pusher,
                       const std::shared_ptr<PartBunch_t::ParticleContainer_t>& pc);

    /// @brief First half of the leapfrog step: push all active containers.
    void timeIntegration1(BorisPusher& pusher);

    /// @brief Second half: kick then push all active containers.
    void timeIntegration2(BorisPusher& pusher);

    /**
     * @brief Self-fields in beam frame (primary container); optional binary repartition.
     * @param step Global step index (used for repartition cadence).
     */
    void computeSpaceChargeFields(unsigned long long step);

    /// @brief Apply external fields from elements intersecting each active container.
    /// @param oth Orbit threader for element queries.
    void computeExternalFields(OrbitThreader& oth);

    /// @brief Emit macroparticles from configured samplers per container.
    /// @param t  Bunch time (s).
    /// @param dt Global time step (s).
    void emitFromEmissionSources(double t, double dt);

    /// @brief Zero E and B on all active particle containers.
    void resetFields();

    /// @brief Set bunch dt from StepSizeConfig and copy to all container dt views.
    void changeDT();

    /// @brief Reset per-particle dt views to the current global bunch dt.
    void setTime();

private:

    /// @brief Update reference trajectories and lab/reference coordinate transforms.
    void updateReference(const BorisPusher& pusher);

    /// @brief Advance reference positions/momenta through the beamline for one step.
    void updateReferenceParticles(const BorisPusher& pusher);

    /// @brief Refresh each container's reference-to-lab transform from current state.
    void updateRefToLabCSTrafo();

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

    /// @brief Load REPARTFREQ and related options from input.
    void setOptionalVariables();

    /**
     * @brief Whether tracking should stop at end-of-line (global reduction).
     * @param globalBoundingBox Spatial bounds from the orbit threader.
     */
    bool hasEndOfLineReached(const BoundingBox& globalBoundingBox);

    /// @brief Trigger binary repartition for the field solver if configured.
    void doBinaryRepartition();

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

    /// @brief Integrate references in time until path length reaches zstart_m.
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
    itsOpalBeamline_m.visit(cav, *this, itsBunch_m);
}

inline void ParallelTracker::visitDrift(const Drift& drift) {
    itsOpalBeamline_m.visit(drift, *this, itsBunch_m);
}

inline void ParallelTracker::visitMarker(const Marker& marker) {
    itsOpalBeamline_m.visit(marker, *this, itsBunch_m);
}

inline void ParallelTracker::visitMultipole(const Multipole& mult) {
    itsOpalBeamline_m.visit(mult, *this, itsBunch_m);
}

inline void ParallelTracker::visitMultipoleT(const MultipoleT& mult) {
    itsOpalBeamline_m.visit(mult, *this, itsBunch_m);
}

inline void ParallelTracker::visitRFCavity(const RFCavity& as) {
    itsOpalBeamline_m.visit(as, *this, itsBunch_m);
}

inline void ParallelTracker::visitTravelingWave(const TravelingWave& tw) {
    itsOpalBeamline_m.visit(tw, *this, itsBunch_m);
}

inline void ParallelTracker::visitSolenoid(const Solenoid& so) {
    itsOpalBeamline_m.visit(so, *this, itsBunch_m);
}

#endif  // OPALX_ParallelTracker_HH
