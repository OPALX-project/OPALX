//
// Class ParallelTracker
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
#ifndef OPAL_ParallelTracker_HH
#define OPAL_ParallelTracker_HH

#include "Algorithms/StepSizeConfig.h"
#include "Algorithms/Tracker.h"
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

class ParallelTracker : public Tracker {

    DataSink* itsDataSink_m;

    OpalBeamline itsOpalBeamline_m;

    /*
      Ring specifics
    */

    // we store a pointer explicitly to the Ring
    Ring* opalRing_m;

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
    unsigned long long repartFreq_m;

    unsigned int emissionSteps_m;

    size_t numParticlesInSimulation_m;

    IpplTimings::TimerRef timeIntegrationTimer1_m;
    IpplTimings::TimerRef timeIntegrationTimer2_m;
    IpplTimings::TimerRef fieldEvaluationTimer_m;
    IpplTimings::TimerRef WakeFieldTimer_m;
    IpplTimings::TimerRef PluginElemTimer_m;
    IpplTimings::TimerRef BinRepartTimer_m;
    IpplTimings::TimerRef OrbThreader_m;
    
    std::set<ParticleMatterInteractionHandler*> activeParticleMatterInteractionHandlers_m;
    bool particleMatterStatus_m;

public:
    typedef std::vector<double> dvector_t;
    typedef std::vector<int> ivector_t;
    typedef std::pair<double[8], Component*> element_pair;
    typedef std::pair<ElementType, element_pair> type_pair;
    typedef std::list<type_pair*> beamline_list;

    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is initially empty.
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ParallelTracker(const Beamline& bl, const PartData& data, bool revBeam, bool revTrack);

    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is taken from [b]bunch[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ParallelTracker(
        const Beamline& bl, PartBunch_t* bunch, DataSink& ds, const PartData& data, bool revBeam,
        bool revTrack, const std::vector<unsigned long long>& maxSTEPS, double zstart,
        const std::vector<double>& zstop, const std::vector<double>& dt);

    virtual ~ParallelTracker();

    /// Apply the algorithm to the top-level beamline.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void execute();

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

    // made following public: __host__ __device__ lambda cannot have private or protected access within its class
    void kickParticles(const BorisPusher& pusher);

    void pushParticles(const BorisPusher& pusher);

    void timeIntegration2(BorisPusher& pusher);

    void changeDT(bool backTrack = false);

    void computeSpaceChargeFields(unsigned long long step);

    void setTime();
private:
    // Not implemented.
    ParallelTracker();
    ParallelTracker(const ParallelTracker&);
    void operator=(const ParallelTracker&);

    /*
      Ring specifics
    */

    unsigned turnnumber_m;

    std::list<Component*> myElements;
    beamline_list FieldDimensions;
    /// The reference variables
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

    void buildupFieldList(double BcParameter[], ElementType elementType, Component* elptr);

    std::vector<PluginElement*> pluginElements_m;





    /********************** END VARIABLES ***********************************/

    void updateReferenceParticle(const BorisPusher& pusher);

    void writePhaseSpace(const long long step, bool psDump, bool statDump);

    /********** BEGIN AUTOPHSING STUFF **********/
    void updateRFElement(std::string elName, double maxPhi);
    void printRFPhases();
    void saveCavityPhases();
    void restoreCavityPhases();
    /************ END AUTOPHSING STUFF **********/

    void prepareSections();

    void timeIntegration1(BorisPusher& pusher);
    void selectDT(bool backTrack = false);
    void emitParticles(long long step);
    void computeExternalFields(OrbitThreader& oth);
    void computeWakefield(IndexMap::value_t& elements);
    void computeParticleMatterInteraction(IndexMap::value_t elements, OrbitThreader& oth);

    // void prepareOpalBeamlineSections();
    void dumpStats(long long step, bool psDump, bool statDump);
    void setOptionalVariables();
    bool hasEndOfLineReached(const BoundingBox& globalBoundingBox);
    void handleRestartRun();

    void doBinaryRepartition();

    void transformBunch(const CoordinateSystemTrafo& trafo);

    void updateReference(const BorisPusher& pusher);
    void updateRefToLabCSTrafo();
    void applyFractionalStep(const BorisPusher& pusher, double tau);
    void findStartPosition(const BorisPusher& pusher);
    void autophaseCavities(const BorisPusher& pusher);

    bool applyPluginElements(const double dt);
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

inline void ParallelTracker::kickParticles(const BorisPusher& /*pusher*/) {
    
    // auto Rview  = itsBunch_m->getParticleContainer()->R.getView();
    auto Pview  = itsBunch_m->getParticleContainer()->P.getView();
    // auto dtview = itsBunch_m->getParticleContainer()->dt.getView(); 
    auto Efview = itsBunch_m->getParticleContainer()->E.getView();
    auto Bfview = itsBunch_m->getParticleContainer()->B.getView();

    /// \todo Apparently, we want mass in eV and charge in elementary charges here to match OPAL's BorisPusher
    double mass = itsBunch_m->getMassPerParticle(); // itsReference.getM();
    double charge = itsBunch_m->getChargePerParticle(); // itsReference.getQ();

    Kokkos::parallel_for(
        "kickParticles", ippl::getRangePolicy(Pview),
        KOKKOS_LAMBDA(const size_t i) {
            /**
             *
             * Implementation follows chapter 4-4, pp. 61–63, from:
             * Birdsall, C. K. and Langdon, A. B. (1985). Plasma Physics via Computer Simulation.
             *
             * Up to finite-precision effects, the new implementation is equivalent to the old one,
             * but uses fewer floating-point operations.
             *
             * The relativistic variant implemented below is described in chapter 15-4, pp. 356–357.
             * However, since different units are used here, small modifications are required.
             * The relativistic variant can be derived from the nonrelativistic one by replacing:
             *   mass
             * with:
             *   gamma * rest mass
             * and transforming the units accordingly.
             *
             * Parameters:
             *   R      Scaled position, R = x / (c * dt). Not used here.
             *   P      Scaled velocity, P = (v / c) * gamma.
             *   Ef     Electric field.
             *   Bf     Magnetic field.
             *   dt     Timestep.
             *   mass   Rest energy, i.e., rest mass * c^2.
             *   charge Particle charge.
             *
             */
           
            Vector_t<double, 3> p = Pview(i); 
            pusher.kick(0, p, Efview(i), Bfview(i), 0, mass, charge); /// \todo might want to remove dt and R altogether from the kick!
            Pview(i) = p; 
        });
        
    /// \todo unnecessary update? kick does not modify positions
    itsBunch_m->getParticleContainer()->update();
    ippl::Comm->barrier();
    *gmsg << "passed BorisPusher argument not used in ParallelTracker::kickParticles" << endl;
}

inline void ParallelTracker::pushParticles(const BorisPusher& /*pusher*/) {

    /// \todo use false for now, since I am not sure how well integrated "dt_per_particle" is (needs to be consistent with particle emission later!).
    itsBunch_m->switchToUnitlessPositions(false);

    auto Rview  = itsBunch_m->getParticleContainer()->R.getView();
    auto Pview  = itsBunch_m->getParticleContainer()->P.getView();
    //auto dtview = itsBunch_m->getParticleContainer()->dt.getView();

    Kokkos::parallel_for(
        "pushParticles", ippl::getRangePolicy(Rview),
        KOKKOS_LAMBDA(const size_t i) {
            /** 
             * \f[ \vec{x}_{n+1/2} = \vec{x}_{n} + \frac{1}{2}\vec{v}_{n-1/2}\quad (= \vec{x}_{n} +
             * \frac{\Delta t}{2} \frac{\vec{\beta}_{n-1/2}\gamma_{n-1/2}}{\gamma_{n-1/2}}) \f]
             *
             * \code
             * R[i] += 0.5 * P[i] * recpgamma;
             * \endcode
             */
            /// \todo check +-
            
            Vector_t<double, 3> x = Rview(i); 
            pusher.push(x, Pview(i), 0); // this 0 is "dt" that is not used with unitless positions!
            Rview(i) = x;
        });

    itsBunch_m->switchOffUnitlessPositions(false);
    itsBunch_m->getParticleContainer()->update();
    ippl::Comm->barrier();
    *gmsg << "pusher not used" << endl;
}

#endif  // OPAL_ParallelTracker_HH
