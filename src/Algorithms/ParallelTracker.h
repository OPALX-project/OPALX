//
// Class ParallelTracker
//   The visitor class for tracking particles with time as independent
//   variable.
//
// Copyright (c) 200x - 2014, Christof Kraus, Paul Scherrer Institut, Villigen PSI, Switzerland
//               2015 - 2016, Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
//               2017 - 2020, Christof Metzger-Kraus
//               2025 - present, Ryan Ammann, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#ifndef OPAL_ParallelTracker_HH
#define OPAL_ParallelTracker_HH

#include "Algorithms/StepSizeConfig.h"
#include "Algorithms/Tracker.h"
#include "Algorithms/OrbitThreader.h"

#include "Steppers/BorisPusher.h"

#include "Structure/DataSink.h"

#include "AbsBeamline/ConstantEFieldCavity.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/MultipoleT.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/Solenoid.h"

#include "Beamlines/Beamline.h"

#include "Distribution/SamplingBase.hpp"

#include "Elements/OpalBeamline.h"

#include <memory>
#include <vector>

class ParticleMatterInteractionHandler;
class PluginElement;

/**
 * @brief Implements the simulation loop
 * 
 * @note ParallelTracker implements the main simulation loop. 
 * TRACK and RUN commands in the inputfile invoke the creation of the 
 * ParalleTracker object and the execution of ParallelTracker::execute()
 */
class ParallelTracker : public Tracker {
private:
    /* ============================= Variables ============================= */
    
    // Responsible for writing beam statistics
    std::shared_ptr<DataSink> itsDataSink_m;

    // Beamline Object which holds a list of pointers to beamline components
    OpalBeamline itsOpalBeamline_m;

    // Flag indicating if End Of Line has been reached
    bool globalEOL_m;

    // Starting position
    double zstart_m;

    /**   
     * stores informations where to change the time step and
     * where to stop the simulation,
     * the time step sizes and
     * the number of time steps with each configuration 
     */
    StepSizeConfig stepSizes_m;

    // The current time stepsize dt (controlled by StepSizeConfig)
    double dtCurrentTrack_m;

    // Controls the frequency of load balancing 
    unsigned long long repartFreq_m;

    /// Per-container time-dependent (emitting) sources; emitParticles(t, dt) each step.
    std::vector<std::vector<std::shared_ptr<SamplingBase>>> emittingSamplers_m;

    /* ============================== Timers =============================== */
    IpplTimings::TimerRef timeIntegrationTimer1_m;
    IpplTimings::TimerRef timeIntegrationTimer2_m;
    IpplTimings::TimerRef fieldEvaluationTimer_m;
    IpplTimings::TimerRef WakeFieldTimer_m;
    IpplTimings::TimerRef PluginElemTimer_m;
    IpplTimings::TimerRef BinRepartTimer_m;
    IpplTimings::TimerRef OrbThreader_m;

public:
    /* ============================ Constructors =========================== */

    /** 
     * Constructor for single track
     * The Beamline object "bl"
     * If "revBeam" is true, the beam runs from s = C to s = 0.
     * OPAL-T parallel tracking is forward (beam) direction only (revTrack is not used).
    */ 
    explicit ParallelTracker(const Beamline& bl, bool revBeam);

    /** 
     * Constructor for single track
     * The Beamline object "bl", the ParticleBunch "bunch", the DataSink "ds".
     * If "revBeam" is true, the beam runs from s = C to s = 0.
     * Vector of maxSteps per track; starting position "zstart"; zstop and dt per segment.
     * Optional per-container emitting samplers (emitParticles(t, dt) called each step)
    */
    explicit ParallelTracker(const Beamline& bl, std::shared_ptr<PartBunch_t> bunch,
        const std::shared_ptr<DataSink>& ds, bool revBeam,
        const std::vector<unsigned long long>& maxSTEPS, 
        double zstart, const std::vector<double>& zstop, 
        const std::vector<double>& dt,
        const std::vector<std::vector<std::shared_ptr<SamplingBase>>>& emittingSamplers = {});

    // Destructor
    virtual ~ParallelTracker();

    /* ========================= Visit Functions =========================== */
    /// Apply the algorithm to a beam line.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void visitBeamline(const Beamline&);

    /// Apply the algorithm to an arbitrary component.
    virtual void visitComponent(const Component&);

    /// Apply the algorithm to a constant E-field cavity element.
    virtual void visitConstantEFieldCavity(const ConstantEFieldCavity&);

    /// Apply the algorithm to a drift space.
    virtual void visitDrift(const Drift&);

    /// Apply the algorithm to a marker.
    virtual void visitMarker(const Marker&);

    /// Apply the algorithm to a multipole.
    virtual void visitMultipole(const Multipole&);

    /// Apply the algorithm to an arbitrary multipole.
    virtual void visitMultipoleT(const MultipoleT&);

    /// Apply the algorithm to a RF cavity.
    virtual void visitRFCavity(const RFCavity&);

    /// Apply the algorithm to a RF cavity.
    virtual void visitSolenoid(const Solenoid&);

    /* ========================= Start Simulation ========================== */
    virtual void execute();

    /* ========================= PIC Functions ============================= */
    void kickParticles(const BorisPusher& pusher,
                       const std::shared_ptr<PartBunch_t::ParticleContainer_t>& pc);
    void pushParticles(const BorisPusher& pusher,
                       const std::shared_ptr<PartBunch_t::ParticleContainer_t>& pc);
    void timeIntegration1(BorisPusher& pusher);
    void timeIntegration2(BorisPusher& pusher);
    void computeSpaceChargeFields(unsigned long long step);    
    void computeExternalFields(OrbitThreader& oth);
    void emitFromEmissionSources(double t, double dt);
    void resetFields();
    /* ===================================================================== */ 
    /* =========================== Functions =============================== */
    // Control the dtview of the bunch
    void changeDT();
    void setTime();

private:
    /// @brief Update ref particles and trafos 
    void updateReference(const BorisPusher& pusher);

    /// @brief Update ref particles positions
    void updateReferenceParticles(const BorisPusher& pusher);

    /// @brief Update ref particle trafos
    void updateRefToLabCSTrafo();

    /// @brief Write phase space
    void writePhaseSpace(const long long step, bool psDump, bool statDump);

    /// @brief dump stats
    void dumpStats(long long step, bool psDump, bool statDump);

    /// @brief Prepare all beamline elements
    void prepareSections();

    /// @brief Sets itsBunch_m::dt to dtCurrentTrack_m 
    void selectDT();

    /// @brief void prepareOpalBeamlineSections();
    void setOptionalVariables();

    /// @brief Checks if the reference particle has reached the end of the beamline
    bool hasEndOfLineReached(const BoundingBox& globalBoundingBox);

    /// @brief Load balancing
    void doBinaryRepartition();
    void computeInitialBounds(Vector_t<double, 3>& rmin, Vector_t<double, 3>& rmax);
    void printInitialContainerRefs(Inform& m) const;

    /// @brief Finds start for reference particle
    void findStartPositions(const BorisPusher& pusher);

    /* ========================== Autophasing ============================== */
    /// @brief Setup for TRAVERLINGWAVE and RFCAVITY
    void autophaseCavities(const BorisPusher& pusher);
    void updateRFElement(std::string elName, double maxPhi);
    void printRFPhases();
    void saveCavityPhases();
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

inline void ParallelTracker::visitSolenoid(const Solenoid& so) {
    itsOpalBeamline_m.visit(so, *this, itsBunch_m);
}

#endif  // OPAL_ParallelTracker_HH
