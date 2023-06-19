//
// Class ParallelTTracker
//   OPAL-T tracker.
//   The visitor class for tracking particles with time as independent
//   variable.
//
// Copyright (c) 200x - 2014, Christof Kraus, Paul Scherrer Institut, Villigen PSI, Switzerland
//               2015 - 2016, Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
//               2017 - 2020, Christof Metzger-Kraus
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
#ifndef OPAL_ParallelTTracker_HH
#define OPAL_ParallelTTracker_HH

#include "Algorithms/Tracker.h"
#include "Steppers/BorisPusher.h"
#include "Structure/DataSink.h"
#include "Algorithms/StepSizeConfig.h"

#include "BasicActions/Option.h"
#include "Utilities/Options.h"

#include "Physics/Physics.h"

#include "Algorithms/OrbitThreader.h"
#include "Algorithms/IndexMap.h"

#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/MultipoleT.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/TravelingWave.h"

#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"

#include <list>
#include <vector>

class ParticleMatterInteractionHandler;

class ParallelTTracker: public Tracker {

public:
    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is initially empty.
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ParallelTTracker(const Beamline &bl,
                              const PartData &data,
                              bool revBeam,
                              bool revTrack);

    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is taken from [b]bunch[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ParallelTTracker(const Beamline &bl,
                              PartBunchBase<double, 3> *bunch,
                              DataSink &ds,
                              const PartData &data,
                              bool revBeam,
                              bool revTrack,
                              const std::vector<unsigned long long> &maxSTEPS,
                              double zstart,
                              const std::vector<double> &zstop,
                              const std::vector<double> &dt);


    virtual ~ParallelTTracker();

    /// Apply the algorithm to the top-level beamline.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void execute();

    /// Apply the algorithm to a beam line.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void visitBeamline(const Beamline &);

    /// Apply the algorithm to a drift space.
    virtual void visitDrift(const Drift &);

    /// Apply the algorithm to a marker.
    virtual void visitMarker(const Marker &);

    /// Apply the algorithm to a multipole.
    virtual void visitMultipole(const Multipole &);

    /// Apply the algorithm to an arbitrary multipole.
    virtual void visitMultipoleT(const MultipoleT &);

    /// Apply the algorithm to a RF cavity.
    virtual void visitRFCavity(const RFCavity &);

    /// Apply the algorithm to a traveling wave.
    virtual void visitTravelingWave(const TravelingWave &);

private:

    // Not implemented.
    ParallelTTracker();
    ParallelTTracker(const ParallelTTracker &);
    void operator=(const ParallelTTracker &);

    /******************** STATE VARIABLES ***********************************/

    DataSink *itsDataSink_m;

    OpalBeamline itsOpalBeamline_m;

    bool globalEOL_m;

    bool wakeStatus_m;

    bool deletedParticles_m;

    WakeFunction* wakeFunction_m;

    double pathLength_m;

    /// where to start
    double zstart_m;

    /// stores informations where to change the time step and
    /// where to stop the simulation,
    /// the time step sizes and
    /// the number of time steps with each configuration
    StepSizeConfig stepSizes_m;

    double dtCurrentTrack_m;

    // This variable controls the minimal number of steps of emission (using bins)
    // before we can merge the bins
    int minStepforReBin_m;

    // this variable controls the minimal number of steps until we repartition the particles
    unsigned int repartFreq_m;

    unsigned int emissionSteps_m;

    size_t numParticlesInSimulation_m;

    IpplTimings::TimerRef timeIntegrationTimer1_m;
    IpplTimings::TimerRef timeIntegrationTimer2_m;
    IpplTimings::TimerRef fieldEvaluationTimer_m ;
    IpplTimings::TimerRef BinRepartTimer_m;
    IpplTimings::TimerRef WakeFieldTimer_m;

    std::set<ParticleMatterInteractionHandler*> activeParticleMatterInteractionHandlers_m;
    bool particleMatterStatus_m;

    /********************** END VARIABLES ***********************************/

    void kickParticles(const BorisPusher &pusher);
    void pushParticles(const BorisPusher &pusher);
    void updateReferenceParticle(const BorisPusher &pusher);

    void writePhaseSpace(const long long step, bool psDump, bool statDump);

    /********** BEGIN AUTOPHSING STUFF **********/
    void updateRFElement(std::string elName, double maxPhi);
    void printRFPhases();
    void saveCavityPhases();
    void restoreCavityPhases();
    /************ END AUTOPHSING STUFF **********/

    void prepareSections();

    void timeIntegration1(BorisPusher & pusher);
    void timeIntegration2(BorisPusher & pusher);
    void selectDT(bool backTrack = false);
    void changeDT(bool backTrack = false);
    void emitParticles(long long step);
    void computeExternalFields(OrbitThreader &oth);
    void computeWakefield(IndexMap::value_t &elements);
    void computeParticleMatterInteraction(IndexMap::value_t elements, OrbitThreader &oth);

    void computeSpaceChargeFields(unsigned long long step);
    // void prepareOpalBeamlineSections();
    void dumpStats(long long step, bool psDump, bool statDump);
    void setOptionalVariables();
    bool hasEndOfLineReached(const BoundingBox& globalBoundingBox);
    void handleRestartRun();

    void setTime();
    void doBinaryRepartition();

    void transformBunch(const CoordinateSystemTrafo &trafo);

    void updateReference(const BorisPusher &pusher);
    void updateRefToLabCSTrafo();
    void applyFractionalStep(const BorisPusher &pusher, double tau);
    void findStartPosition(const BorisPusher &pusher);
    void autophaseCavities(const BorisPusher &pusher);

};


inline void ParallelTTracker::visitDrift(const Drift &drift) {
    itsOpalBeamline_m.visit(drift, *this, itsBunch_m);
}

inline void ParallelTTracker::visitMarker(const Marker &marker) {
    itsOpalBeamline_m.visit(marker, *this, itsBunch_m);
}

inline void ParallelTTracker::visitMultipole(const Multipole &mult) {
    itsOpalBeamline_m.visit(mult, *this, itsBunch_m);
}

inline void ParallelTTracker::visitMultipoleT(const MultipoleT &mult) {
    itsOpalBeamline_m.visit(mult, *this, itsBunch_m);
}


inline void ParallelTTracker::visitRFCavity(const RFCavity &as) {
    itsOpalBeamline_m.visit(as, *this, itsBunch_m);
}

inline void ParallelTTracker::visitTravelingWave(const TravelingWave &as) {
    itsOpalBeamline_m.visit(as, *this, itsBunch_m);
}

inline void ParallelTTracker::kickParticles(const BorisPusher &pusher) {
    int localNum = itsBunch_m->getLocalNum();
    for (int i = 0; i < localNum; ++i)
        pusher.kick(itsBunch_m->R[i], itsBunch_m->P[i], itsBunch_m->Ef[i], itsBunch_m->Bf[i], itsBunch_m->dt[i]);
}

inline void ParallelTTracker::pushParticles(const BorisPusher &pusher) {
    itsBunch_m->switchToUnitlessPositions(true);

    for (unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {
        pusher.push(itsBunch_m->R[i], itsBunch_m->P[i], itsBunch_m->dt[i]);
    }
    itsBunch_m->switchOffUnitlessPositions(true);
}

#endif // OPAL_ParallelTTracker_HH