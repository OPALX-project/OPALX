//
// Class ParallelTracker
//   OPAL-T tracker.
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
#include "Algorithms/CoordinateSystemTrafo.h"
#include "Steppers/BorisPusher.h"
#include "Structure/DataSink.h"

#include "BasicActions/Option.h"
#include "Utilities/Options.h"

#include "Physics/Physics.h"

#include "Algorithms/IndexMap.h"
#include "Algorithms/OrbitThreader.h"

#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/MultipoleT.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/Ring.h"
#include "AbsBeamline/ScalingFFAMagnet.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/TravelingWave.h"
#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"

#include <list>
#include <memory>
#include <tuple>
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
    DataSink* itsDataSink_m;

    // Beamline Object which holds a list of pointers to beamline components
    OpalBeamline itsOpalBeamline_m;

    // Pointer to OpalRing
    Ring* opalRing_m;

    // Flag indicating if End Of Line has been reached
    bool globalEOL_m;

    // Flag indicating if particles are lost
    bool deletedParticles_m;

    // Path along beamline
    double pathLength_m;

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

    // This variable controls the minimal number of steps of 
    // emission (using bins) before we can merge the bins
    int minStepforReBin_m;

    // Controls the frequency of load balancing 
    unsigned long long repartFreq_m;

    // Total number of particles in the whole simulation
    size_t numParticlesInSimulation_m;
    /* ===================================================================== */
    /* ============================== Timers =============================== */
    IpplTimings::TimerRef timeIntegrationTimer1_m;
    IpplTimings::TimerRef timeIntegrationTimer2_m;
    IpplTimings::TimerRef fieldEvaluationTimer_m;
    IpplTimings::TimerRef WakeFieldTimer_m;
    IpplTimings::TimerRef PluginElemTimer_m;
    IpplTimings::TimerRef BinRepartTimer_m;
    IpplTimings::TimerRef OrbThreader_m;
    /* ===================================================================== */
    /* ========================== Ring Variables =========================== */
    typedef std::pair<double[8], Component*> element_pair;
    typedef std::pair<ElementType, element_pair> type_pair;
    typedef std::list<type_pair*> beamline_list;
    
    unsigned turnnumber_m;

    std::list<Component*> myElements;
    beamline_list FieldDimensions;

    double bega;
    double referenceR;
    double referenceTheta;
    double referenceZ = 0.0;

    double referencePr;
    double referencePt;
    double referencePz = 0.0;
    double referencePtot;

    double referencePsi;
    double referencePhi;

    double sinRefTheta_m;
    double cosRefTheta_m; 
    
    std::vector<PluginElement*> pluginElements_m;
    /* ===================================================================== */ 
    /* ========================== NOT IMPLEMENTED ========================== */
    // Particle - Matter interaction
    std::set<ParticleMatterInteractionHandler*> 
        activeParticleMatterInteractionHandlers_m;
    bool particleMatterStatus_m;
   
    // Does nothing ...
    unsigned int emissionSteps_m;
    
    // Wakefield stuff - Does nothing... 
    bool wakeStatus_m;
    WakeFunction* wakeFunction_m;
    /* ===================================================================== */ 
public:
    /* ============================ Constructors =========================== */
    /*
     * Constructor for single track
     * The Beamline object "bl"
     * The reference data object "data" 
     * If "revBeam" is true, the beam runs from s = C to s = 0.
     * If "revTrack" is true, we track against the beam.
    */ 
    explicit ParallelTracker(const Beamline& bl, const PartData& data, 
        bool revBeam, bool revTrack);

    /*
     * Constructor for single track
     * The Beamline object "bl"
     * The ParticleBunch "bunch"
     * The DataSink "ds"
     * The reference data object "data" 
     * If "revBeam" is true, the beam runs from s = C to s = 0.
     * If "revTrack" is true, we track against the beam.
     * Vector of maxSteps per track
     * Starting position "zstart"
     * Vector of ends of the individual tracks
     * Vector of different timesteps for individual tracks
    */
    explicit ParallelTracker(const Beamline& bl, PartBunch_t* bunch, 
        DataSink& ds, const PartData& data, bool revBeam,
        bool revTrack, const std::vector<unsigned long long>& maxSTEPS, 
        double zstart, const std::vector<double>& zstop, 
        const std::vector<double>& dt);

    // Destructor
    virtual ~ParallelTracker();
    /* ===================================================================== */
    /* ========================= Visit Functions =========================== */
    /// Apply the algorithm to a beam line.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void visitBeamline(const Beamline&);

    /// Apply the algorithm to a drift space.
    virtual void visitDrift(const Drift&);

    /// Apply the algorithm to a ring
    virtual void visitRing(const Ring& ring);

    /// Apply the algorithm to a marker.
    virtual void visitMarker(const Marker&);

    /// Apply the algorithm to a multipole.
    virtual void visitMultipole(const Multipole&);

    /// Apply the algorithm to an arbitrary multipole.
    virtual void visitMultipoleT(const MultipoleT&);

    /// Apply the algorithm to a offset (placement).
    virtual void visitOffset(const Offset&);

    /// Apply the algorithm to a RF cavity.
    virtual void visitRFCavity(const RFCavity&);

    /// Apply the algorithm to a RF cavity.
    virtual void visitSolenoid(const Solenoid&);

    /// Apply the algorithm to a traveling wave.
    virtual void visitTravelingWave(const TravelingWave&);

    /// Apply the algorithm to a scaling FFA.
    virtual void visitScalingFFAMagnet(const ScalingFFAMagnet& bend);

    /// Apply the algorithm to a vertical FFA magnet.
    virtual void visitVerticalFFAMagnet(const VerticalFFAMagnet& bend);
    /* ===================================================================== */
    /* ========================= Start Simulation ========================== */
    virtual void execute();
    /* ===================================================================== */
    /* ========================= PIC Functions ============================= */
    void kickParticles(const BorisPusher& pusher);
    void pushParticles(const BorisPusher& pusher);
    void timeIntegration1(BorisPusher& pusher);
    void timeIntegration2(BorisPusher& pusher);
    void computeSpaceChargeFields(unsigned long long step);    
    void computeExternalFields(OrbitThreader& oth);
    void resetFields();
    /* ===================================================================== */ 
    /* =========================== Functions =============================== */
    // Control the dtview of the bunch
    void changeDT(bool backTrack = false);
    void setTime();

private:
    // Reference Particle 
    void updateReferenceParticle(const BorisPusher& pusher);
    void updateReference(const BorisPusher& pusher);
    void updateRefToLabCSTrafo();

    // File writing
    void writePhaseSpace(const long long step, bool psDump, bool statDump);
    void dumpStats(long long step, bool psDump, bool statDump);

    // Visits all elements inside itsBeamline_m
    // Sorts elements in itsOpalBeamline_m
    // Sets itsOpalBeamline_m::prepared_m = true
    void prepareSections();

    // Sets itsBunch_m::dt to dtCurrentTrack_m 
    void selectDT(bool backTrack = false);

    // void prepareOpalBeamlineSections();
    void setOptionalVariables();

    // Checks if the reference particle has reached the end of the beamline
    bool hasEndOfLineReached(const BoundingBox& globalBoundingBox);

    // Load balancing
    void doBinaryRepartition();

    // Applies a (fractional) step tau
    void applyFractionalStep(const BorisPusher& pusher, double tau);

    // Finds start for reference particle
    void findStartPosition(const BorisPusher& pusher);
    /* ===================================================================== */
    /* ========================== Autophasing ============================== */
    // Setup for TRAVERLINGWAVE and RFCAVITY
    void autophaseCavities(const BorisPusher& pusher);
    void updateRFElement(std::string elName, double maxPhi);
    void printRFPhases();
    void saveCavityPhases();
    void restoreCavityPhases();
    /* ===================================================================== */
    /* ========================== RING FUNCTIONS =========================== */
    void buildupFieldList(double BcParameter[], ElementType elementType, 
        Component* elptr);
    bool applyPluginElements(const double dt);
    /* ===================================================================== */
    /* ========================== NOT IMPLEMENTED ========================== */
    ParallelTracker();
    ParallelTracker(const ParallelTracker&);
    void operator=(const ParallelTracker&);
    void emitParticles(long long step);
    void computeWakefield(IndexMap::value_t& elements);
    void computeParticleMatterInteraction(IndexMap::value_t elements, 
        OrbitThreader& oth);
    void handleRestartRun();
    void transformBunch(const CoordinateSystemTrafo& trafo);
    /* ===================================================================== */
};

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

inline void ParallelTracker::visitTravelingWave(const TravelingWave& as) {
    itsOpalBeamline_m.visit(as, *this, itsBunch_m);
}

#endif  // OPAL_ParallelTracker_HH
