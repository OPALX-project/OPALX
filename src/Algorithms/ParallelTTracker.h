#ifndef OPAL_ParallelTTracker_HH
#define OPAL_ParallelTTracker_HH

// ------------------------------------------------------------------------
// $RCSfile: ParallelTTracker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelTTracker
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/Tracker.h"
#include "Algorithms/PartPusher.h"
#include "Structure/DataSink.h"

#include "BasicActions/Option.h"
#include "Utilities/OpalOptions.h"
#include "Utilities/Options.h"

#include "Physics/Physics.h"
#ifdef HAVE_AMR_SOLVER
#include <Amr.H>
#endif

class BMultipoleField;
class PartBunch;
class AlignWrapper;
class BeamBeam;
#include "AbsBeamline/Collimator.h"
class Corrector;
class Degrader;
class Diagnostic;
class Drift;
class ElementBase;
class Lambertson;
class Marker;
class Monitor;
class Multipole;
class Probe;
class RBend;
class RFCavity;
class TravelingWave;
class RFQuadrupole;
class SBend;
class Separator;
class Septum;
class Solenoid;
class ParallelPlate;
class CyclotronValley;

#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"
#include "Solvers/WakeFunction.hh"

#include <list>
#include <vector>
#include <queue>

class BorisPusher;



// Class ParallelTTracker
// ------------------------------------------------------------------------
/// Track using a time integration scheme
// [p]
// Phase space coordinates numbering:
// [tab 3 b]
// [row]number [&]name          [&]unit  [/row]
// [row]0      [&]$x$           [&]metres [/row]
// [row]1      [&]$p_x/p_r$     [&]1      [/row]
// [row]2      [&]$y$           [&]metres [/row]
// [row]3      [&]$p_y/p_r$     [&]1      [/row]
// [row]4      [&]$v*delta_t$   [&]metres [/row]
// [row]5      [&]$delta_p/p_r$ [&]1      [/row]
// [/tab][p]
// Where $p_r$ is the constant reference momentum defining the reference
// frame velocity, $m$ is the rest mass of the particles, and $v$ is the
// instantaneous velocity of the particle.
// [p]
// Other units used:
// [tab 2 b]
// [row]quantity             [&]unit           [/row]
// [row]reference momentum   [&]electron-volts [/row]
// [row]velocity             [&]metres/second  [/row]
// [row]accelerating voltage [&]volts          [/row]
// [row]separator voltage    [&]volts          [/row]
// [row]frequencies          [&]hertz          [/row]
// [row]phase lags           [&]$2*pi$         [/row]
// [/tab][p]
// Approximations used:
// [ul]
// [li] blah
// [li] blah
// [li] blah
// [/ul]
//
// On going through an element, we use the following steps:
// To complete the map, we propagate the closed orbit and add that to the map.

typedef std::pair<double, Vector_t > PhaseEnT;
// typedef std::pair<std::string, double > MaxPhasesT;

class ParallelTTracker: public Tracker {

public:
    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is initially empty.
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ParallelTTracker(const Beamline &bl, const PartData &data,
                              bool revBeam, bool revTrack, size_t N);

    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is taken from [b]bunch[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ParallelTTracker(const Beamline &bl, PartBunch &bunch, DataSink &ds,
                              const PartData &data, bool revBeam, bool revTrack,
                              const std::vector<unsigned long long> &maxSTEPS, const std::vector<double> &zstop,
                              int timeIntegrator, const std::vector<double> &dt, size_t N);

    /// Constructor
    //  Amr pointer is taken
#ifdef HAVE_AMR_SOLVER
    explicit ParallelTTracker(const Beamline &bl, PartBunch &bunch, DataSink &ds,
                              const PartData &data, bool revBeam, bool revTrack,
                              const std::vector<unsigned long long> &maxSTEPS, const std::vector<double> &zstop,
                              int timeIntegrator, const std::vector<double> &dt, size_t N, Amr* amrptr_in);
#endif

    virtual ~ParallelTTracker();

    virtual void visitAlignWrapper(const AlignWrapper &);

    /// Apply the algorithm to a BeamBeam.
    virtual void visitBeamBeam(const BeamBeam &);

    /// Apply the algorithm to a collimator.
    virtual void visitCollimator(const Collimator &);


    /// Apply the algorithm to a Corrector.
    virtual void visitCorrector(const Corrector &);

    /// Apply the algorithm to a Degrader.
    virtual void visitDegrader(const Degrader &);

    /// Apply the algorithm to a Diagnostic.
    virtual void visitDiagnostic(const Diagnostic &);

    /// Apply the algorithm to a Drift.
    virtual void visitDrift(const Drift &);

    /// Apply the algorithm to a Lambertson.
    virtual void visitLambertson(const Lambertson &);

    /// Apply the algorithm to a Marker.
    virtual void visitMarker(const Marker &);

    /// Apply the algorithm to a Monitor.
    virtual void visitMonitor(const Monitor &);

    /// Apply the algorithm to a Multipole.
    virtual void visitMultipole(const Multipole &);

    /// Apply the algorithm to a Probe.
    virtual void visitProbe(const Probe &);

    /// Apply the algorithm to a RBend.
    virtual void visitRBend(const RBend &);

    /// Apply the algorithm to a RFCavity.
    virtual void visitRFCavity(const RFCavity &);

    /// Apply the algorithm to a RFCavity.
    virtual void visitTravelingWave(const TravelingWave &);

    /// Apply the algorithm to a RFQuadrupole.
    virtual void visitRFQuadrupole(const RFQuadrupole &);

    /// Apply the algorithm to a SBend.
    virtual void visitSBend(const SBend &);

    /// Apply the algorithm to a Separator.
    virtual void visitSeparator(const Separator &);

    /// Apply the algorithm to a Septum.
    virtual void visitSeptum(const Septum &);

    /// Apply the algorithm to a Solenoid.
    virtual void visitSolenoid(const Solenoid &);

    /// Apply the algorithm to a ParallelPlate.
    virtual void visitParallelPlate(const ParallelPlate &);

    /// Apply the algorithm to a CyclotronValley.
    virtual void visitCyclotronValley(const CyclotronValley &);


    //
    // Apply the emission considered schottky effect
    double schottkyLoop(double coeff);

    // Apply the correction to each particle for schottky effect
    void applySchottkyCorrection(PartBunch &itsBunch, int ne, double t, double rescale_coeff);


    /// Apply the algorithm to the top-level beamline.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void execute();

    /// Apply the algorithm to a beam line.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void visitBeamline(const Beamline &);

    /// set Multipacting flag
    inline void setMpacflg(bool mpacflg) {

        mpacflg_m = mpacflg;
    }

private:

    // Not implemented.
    ParallelTTracker();
    ParallelTTracker(const ParallelTTracker &);
    void operator=(const ParallelTTracker &);
#ifdef HAVE_AMR_SOLVER
    Amr* amrptr;
#endif

    /******************** STATE VARIABLES ***********************************/

    PartBunch        *itsBunch;
    DataSink         *itsDataSink_m;
    BoundaryGeometry *bgf_m;

    OpalBeamline itsOpalBeamline_m;
    LineDensity  lineDensity_m;


    Vector_t RefPartR_zxy_m;
    Vector_t RefPartP_zxy_m;
    Vector_t RefPartR_suv_m;
    Vector_t RefPartP_suv_m;


#ifdef OPAL_DKS
  DKSBase dksbase;

  void *r_ptr;
  void *p_ptr;
  void *x_ptr;

  void *lastSec_ptr;
  void *orient_ptr;
  void *dt_ptr;

  int ierr;

  int stream1;
  int stream2;

  unsigned int numDeviceElements;
#endif

    bool globalEOL_m;

    bool wakeStatus_m;

    /*--------- Added by Xiaoying Pang 04/22/2014 ---------------
     * This WakeFunction pointer is used to store a dipole's wake function and
     * to be used in the following drift if CSR calculation is requested in
     * the drift. */
    WakeFunction* wakeFunction_m;

    bool surfaceStatus_m;

    int secondaryFlg_m;

    /// multipacting flag
    bool mpacflg_m;

    bool nEmissionMode_m;

    /// where to stop
    std::queue<double> zStop_m;

    /// The scale factor for dimensionless variables (FIXME: move to PartBunch)
    double scaleFactor_m;

    // Vector of the scale factor for dimensionless variables (FIXME: move to PartBunch)
    Vector_t vscaleFactor_m;

    double recpGamma_m;

    double rescale_coeff_m;

    double dtCurrentTrack_m;
    std::queue<double> dtAllTracks_m;

    double surfaceEmissionStop_m;


    size_t specifiedNPart_m;

    // This variable controls the minimal number of steps of emission (using bins)
    // before we can merge the bins
    int minStepforReBin_m;

    // The space charge solver crashes if we use less than ~10 particles.
    // This variable controls the number of particles to be emitted before we use
    // the space charge solver.
    size_t minBinEmitted_m;

    // this variable controls the minimal number of steps until we repartition the particles
    int repartFreq_m;

    int lastVisited_m;

    /// The number of refinements of the search range of the phase
    int numRefs_m;

    /// ??
    int gunSubTimeSteps_m;

    unsigned int emissionSteps_m;

    /// The maximal number of steps the system is integrated per TRACK
    std::queue<unsigned long long> localTrackSteps_m;

    size_t maxNparts_m;
    size_t numberOfFieldEmittedParticles_m;

    // flag which indicates whether any particle is within the influence of bending element.
    // if this is the case we track the reference particle as if it were a real particle,
    // otherwise the reference particle is defined as the centroid particle of the bunch
    unsigned long bends_m;

    size_t numParticlesInSimulation_m;

    Tenzor<double, 3> space_orientation_m;

    FieldList::iterator currently_ap_cavity_m;

    // Vector of the scale factor for dimensionless variables (FIXME: move to PartBunch)

    IpplTimings::TimerRef timeIntegrationTimer1_m;
    IpplTimings::TimerRef timeIntegrationTimer2_m;
    IpplTimings::TimerRef timeFieldEvaluation_m ;
    IpplTimings::TimerRef BinRepartTimer_m;
    IpplTimings::TimerRef WakeFieldTimer_m;

  IpplTimings::TimerRef timeIntegrationTimer1Loop1_m;
  IpplTimings::TimerRef timeIntegrationTimer1Loop2_m;
  IpplTimings::TimerRef timeIntegrationTimer2Loop1_m;
  IpplTimings::TimerRef timeIntegrationTimer2Loop2_m;

    // 1 --- LF-2 (Boris-Buneman)
    // 3 --- AMTS (Adaptive Boris-Buneman with multiple time stepping)
    int timeIntegrator_m;


    size_t Nimpact_m;
    double SeyNum_m;


    std::map<long, SurfacePhysicsHandler*> sphys_m;


    /********************** END VARIABLES ***********************************/

    int LastVisited;

    // Fringe fields for entrance and exit of magnetic elements.
    void applyEntranceFringe(double edge, double curve,
                             const BMultipoleField &field, double scale);
    void applyExitFringe(double edge, double curve,
                         const BMultipoleField &field, double scale);

    inline double getEnergyMeV(Vector_t p) {
        return (sqrt(dot(p, p) + 1.0) - 1.0) * itsBunch->getM() * 1e-6;
    }

    void kickParticles(const BorisPusher &pusher);
    void kickParticles(const BorisPusher &pusher, const int &flg);
    void updateReferenceParticle();

    void updateSpaceOrientation(const bool &move = false);
    Vector_t TransformTo(const Vector_t &vec, const Vector_t &ori) const;
    Vector_t TransformBack(const Vector_t &vec, const Vector_t &ori) const;

    void kickReferenceParticle(const Vector_t &externalE, const Vector_t &externalB);
    void writePhaseSpace(const long long step, const double &sposRef, bool psDump, bool statDump);

    /********** BEGIN AUTOPHSING STUFF **********/
    void updateRFElement(std::string elName, double maxPhi);
    void printRFPhases();
    void handleAutoPhasing();
    /************ END AUTOPHSING STUFF **********/

    void handleOverlappingMonitors();
    void prepareSections();

    double ptoEMeV(Vector_t p);

    void bgf_main_collision_test();
    void timeIntegration1(BorisPusher & pusher);
    void timeIntegration1_bgf(BorisPusher & pusher);
    void timeIntegration2(BorisPusher & pusher);
    void timeIntegration2_bgf(BorisPusher & pusher);
    void selectDT();
    void changeDT();
    void emitParticles(long long step);
    void computeExternalFields();
    void handleBends();
    void switchElements(double scaleMargin = 3.0);
    void computeSpaceChargeFields();
    void prepareOpalBeamlineSections();
    void dumpStats(long long step, bool psDump, bool statDump);
    void setOptionalVariables();
    bool hasEndOfLineReached();
    void doSchottyRenormalization();
    void setupSUV(bool updateReference = true);
    void handleRestartRun();
    void prepareEmission();
    void setTime();
    void initializeBoundaryGeometry();
    void doBinaryRepartition();
    void executeDefaultTracker();
#ifdef HAVE_AMR_SOLVER
    void executeAMRTracker();
#endif
    void executeAMTSTracker();
    void push(double h);
    void kick(double h, bool avoidGammaCalc = false);
    void computeExternalFields_AMTS();
    void borisExternalFields(double h, bool isFirstSubstep, bool isLastSubstep);
    double calcG(); // Time step chooser for adaptive variant
    Vector_t calcMeanR() const;
    Vector_t calcMeanP() const;
    double pSqr(Vector_t p1, Vector_t p2);
};

inline double ParallelTTracker::pSqr(Vector_t p1, Vector_t p2) {
  return p1[0]*p2[0] + p1[0]*p2[1] + p1[0]*p2[2] + p1[1]*p2[0] + p1[1]*p2[1] + p1[1]*p2[2] + p1[2]*p2[0] + p1[2]*p2[1] + p1[2]*p2[2];
}

inline double ParallelTTracker::ptoEMeV(Vector_t p) {
    return (sqrt(dot(p, p) + 1.0) - 1.0) * itsBunch->getM() * 1e-6;
}


inline void ParallelTTracker::visitBeamline(const Beamline &bl) {
    itsBeamline_m.iterate(*dynamic_cast<BeamlineVisitor *>(this), false);
}

inline void ParallelTTracker::visitAlignWrapper(const AlignWrapper &wrap) {
    itsOpalBeamline_m.visit(wrap, *this, itsBunch);
}

inline void ParallelTTracker::visitBeamBeam(const BeamBeam &bb) {
    itsOpalBeamline_m.visit(bb, *this, itsBunch);
}


inline void ParallelTTracker::visitCollimator(const Collimator &coll) {
    itsOpalBeamline_m.visit(coll, *this, itsBunch);
}


inline void ParallelTTracker::visitCorrector(const Corrector &corr) {
    itsOpalBeamline_m.visit(corr, *this, itsBunch);
}


inline void ParallelTTracker::visitDegrader(const Degrader &deg) {
    itsOpalBeamline_m.visit(deg, *this, itsBunch);
}


inline void ParallelTTracker::visitDiagnostic(const Diagnostic &diag) {
    itsOpalBeamline_m.visit(diag, *this, itsBunch);
}


inline void ParallelTTracker::visitDrift(const Drift &drift) {
    itsOpalBeamline_m.visit(drift, *this, itsBunch);
}


inline void ParallelTTracker::visitLambertson(const Lambertson &lamb) {
    itsOpalBeamline_m.visit(lamb, *this, itsBunch);
}


inline void ParallelTTracker::visitMarker(const Marker &marker) {
    itsOpalBeamline_m.visit(marker, *this, itsBunch);
}


inline void ParallelTTracker::visitMonitor(const Monitor &mon) {
    itsOpalBeamline_m.visit(mon, *this, itsBunch);
}


inline void ParallelTTracker::visitMultipole(const Multipole &mult) {
    itsOpalBeamline_m.visit(mult, *this, itsBunch);
}

inline void ParallelTTracker::visitProbe(const Probe &prob) {
    itsOpalBeamline_m.visit(prob, *this, itsBunch);
}


inline void ParallelTTracker::visitRBend(const RBend &bend) {
    itsOpalBeamline_m.visit(bend, *this, itsBunch);
}


inline void ParallelTTracker::visitRFCavity(const RFCavity &as) {
    itsOpalBeamline_m.visit(as, *this, itsBunch);
}

inline void ParallelTTracker::visitTravelingWave(const TravelingWave &as) {
    itsOpalBeamline_m.visit(as, *this, itsBunch);
}


inline void ParallelTTracker::visitRFQuadrupole(const RFQuadrupole &rfq) {
    itsOpalBeamline_m.visit(rfq, *this, itsBunch);
}

inline void ParallelTTracker::visitSBend(const SBend &bend) {
    itsOpalBeamline_m.visit(bend, *this, itsBunch);
}


inline void ParallelTTracker::visitSeparator(const Separator &sep) {
    itsOpalBeamline_m.visit(sep, *this, itsBunch);
}


inline void ParallelTTracker::visitSeptum(const Septum &sept) {
    itsOpalBeamline_m.visit(sept, *this, itsBunch);
}


inline void ParallelTTracker::visitSolenoid(const Solenoid &solenoid) {
    itsOpalBeamline_m.visit(solenoid, *this, itsBunch);
}


inline void ParallelTTracker::visitParallelPlate(const ParallelPlate &pplate) {
    itsOpalBeamline_m.visit(pplate, *this, itsBunch);
}

inline void ParallelTTracker::visitCyclotronValley(const CyclotronValley &cv) {
    itsOpalBeamline_m.visit(cv, *this, itsBunch);
}

inline void ParallelTTracker::kickParticles(const BorisPusher &pusher) {
    int localNum = itsBunch->getLocalNum();
    for(int i = 0; i < localNum; ++i)
        pusher.kick(itsBunch->R[i], itsBunch->P[i], itsBunch->Ef[i], itsBunch->Bf[i], itsBunch->dt[i]);
    itsBunch->calcBeamParameters();
}

// BoundaryGeometry version of kickParticles function
inline void ParallelTTracker::kickParticles(const BorisPusher &pusher, const int &flg) {
    int localNum = itsBunch->getLocalNum();
    for(int i = 0; i < localNum; ++i) {
        if(itsBunch->TriID[i] == 0) { //TriID[i] will be 0 if particles don't collide the boundary before kick;
            pusher.kick(itsBunch->R[i], itsBunch->P[i], itsBunch->Ef[i], itsBunch->Bf[i], itsBunch->dt[i]);
        }
    }
    itsBunch->calcBeamParameters();
}

inline void ParallelTTracker::updateReferenceParticle() {
    RefPartR_suv_m = itsBunch->get_rmean() * scaleFactor_m;
    RefPartP_suv_m = itsBunch->get_pmean();

    /* Update the position of the reference particle in ZXY-coordinates. The angle between the ZXY- and the SUV-coordinate
     *  system is determined by the momentum of the reference particle. We calculate the momentum of the reference
     *  particle by rotating the centroid momentum (= momentum of the reference particle in the SUV-coordinate system).
     *  Then we push the reference particle with this momentum for half a time step.
     */

    double gamma = sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));

    /* First update the momentum of the reference particle in zxy coordinate system, then update its position     */
    RefPartP_zxy_m = dot(space_orientation_m, RefPartP_suv_m);

    RefPartR_zxy_m += RefPartP_zxy_m * scaleFactor_m / (2. * gamma);
    RefPartR_suv_m += RefPartP_suv_m * scaleFactor_m / (2. * gamma);

}

inline void ParallelTTracker::updateSpaceOrientation(const bool &move) {
    /*
       Update the position of the reference particle in
       ZXY-coordinates. The angle between the ZXY- and the
       SUV-coordinate system is determined by the momentum of the
       reference particle. We calculate the momentum of the reference
       particle by rotating the centroid momentum (= momentum of the
       reference particle in the SUV-coordinate system). Then we push
       the reference particle with this momentum for half a time step.
     */
    itsBunch->calcBeamParameters();
    RefPartR_suv_m = itsBunch->get_rmean();
    RefPartP_suv_m = itsBunch->get_pmean();

    double AbsMomentum = sqrt(dot(RefPartP_suv_m, RefPartP_suv_m));
    double AbsMomentumProj = sqrt(RefPartP_suv_m(0) * RefPartP_suv_m(0) + RefPartP_suv_m(2) * RefPartP_suv_m(2));
    RefPartP_zxy_m *= AbsMomentum / sqrt(dot(RefPartP_zxy_m, RefPartP_zxy_m));

    space_orientation_m(0, 0) = RefPartP_zxy_m(2) / AbsMomentumProj;
    space_orientation_m(0, 1) = -RefPartP_zxy_m(0) * RefPartP_zxy_m(1) / (AbsMomentum * AbsMomentumProj);
    space_orientation_m(0, 2) = RefPartP_zxy_m(0) / AbsMomentum;
    space_orientation_m(1, 0) = 0.;
    space_orientation_m(1, 1) = AbsMomentumProj / AbsMomentum;
    space_orientation_m(1, 2) = RefPartP_zxy_m(1) / AbsMomentum;
    space_orientation_m(2, 0) = -RefPartP_zxy_m(0) / AbsMomentumProj;
    space_orientation_m(2, 1) = -RefPartP_zxy_m(2) * RefPartP_zxy_m(1) / (AbsMomentum * AbsMomentumProj);
    space_orientation_m(2, 2) = RefPartP_zxy_m(2) / AbsMomentum;

    Vector_t EulerAngles;

    if(RefPartP_suv_m(2) < 1.0e-12) {
        EulerAngles = Vector_t(0.0);
    } else {
        EulerAngles = Vector_t(-atan(RefPartP_suv_m(0) / RefPartP_suv_m(2)),          \
                               -asin(RefPartP_suv_m(1) / AbsMomentum),    \
                               0.0);
    }

    // Rotate the local coordinate system of all sections which are online
    Vector_t smin(0.0), smax(0.0);
    itsBunch->get_bounds(smin, smax);
    itsOpalBeamline_m.updateOrientation(EulerAngles, RefPartR_suv_m, smin(2), smax(2));

    itsBunch->rotateAbout(RefPartR_suv_m, RefPartP_suv_m);
    if(move)       // move the bunch such that the new centroid location is at (0,0,z)
        itsBunch->moveBy(Vector_t(-RefPartR_suv_m(0), -RefPartR_suv_m(1), 0.0));
}

inline Vector_t ParallelTTracker::TransformTo(const Vector_t &vec, const Vector_t &ori) const {

    // Rotate to the the element's local coordinate system.
    //
    // 1) Rotate about the z axis by angle negative ori(2).
    // 2) Rotate about the y axis by angle negative ori(0).
    // 3) Rotate about the x axis by angle ori(1).

    const double sina = sin(ori(0));
    const double cosa = cos(ori(0));
    const double sinb = sin(ori(1));
    const double cosb = cos(ori(1));
    const double sinc = sin(ori(2));
    const double cosc = cos(ori(2));

    Vector_t temp(0.0, 0.0, 0.0);

    temp(0) = (cosa * cosc) * vec(0) + (cosa * sinc) * vec(1) - sina *        vec(2);
    temp(1) = (-cosb * sinc - sina * sinb * cosc) * vec(0) + (cosb * cosc - sina * sinb * sinc) * vec(1) - cosa * sinb * vec(2);
    temp(2) = (-sinb * sinc + sina * cosb * cosc) * vec(0) + (sinb * cosc + sina * cosb * sinc) * vec(1) + cosa * cosb * vec(2);

    return temp;
}

inline Vector_t ParallelTTracker::TransformBack(const Vector_t &vec, const Vector_t &ori) const {

    // Rotate out of the element's local coordinate system.
    //
    // 1) Rotate about the x axis by angle negative ori(1).
    // 2) Rotate about the y axis by angle ori(0).
    // 3) Rotate about the z axis by angle ori(3).

    const double sina = sin(ori(0));
    const double cosa = cos(ori(0));
    const double sinb = sin(ori(1));
    const double cosb = cos(ori(1));
    const double sinc = sin(ori(2));
    const double cosc = cos(ori(2));

    Vector_t temp(0.0, 0.0, 0.0);

    temp(0) =  cosa * cosc * vec(0) + (-sina * sinb * cosc - cosb * sinc) * vec(1) + (sina * cosb * cosc - sinb * sinc) * vec(2);
    temp(1) =  cosa * sinc * vec(0) + (-sina * sinb * sinc + cosb * cosc) * vec(1) + (sina * cosb * sinc + sinb * cosc) * vec(2);
    temp(2) = -sina        * vec(0) + (-cosa * sinb) * vec(1) + (cosa * cosb) * vec(2);

    return temp;
}

inline void ParallelTTracker::kickReferenceParticle(const Vector_t &externalE, const Vector_t &externalB) {
    using Physics::c;

    // track reference particle
    Vector_t um = RefPartP_suv_m + 0.5 * itsReference.getQ() * itsBunch->getdT() / itsReference.getM() * c * externalE;
    double gamma = sqrt(1.0 + dot(um, um));

    double tmp = 0.5 * itsReference.getQ() * c * c * itsBunch->getdT() / (itsReference.getM() * gamma);
    Vector_t a = tmp * externalB;

    Vector_t s = um + tmp * cross(um, externalB);

    tmp = 1.0 + dot(a, a);

    um(0) = ((1.0 + a(0) * a(0))    * s(0) + (a(0) * a(1) + a(2)) * s(1) + (a(0) * a(2) - a(1)) * s(2)) / tmp;
    um(1) = ((a(0) * a(1) - a(2)) * s(0) + (1.0 + a(1) * a(1)) * s(1) + (a(1) * a(2) + a(0)) * s(2)) / tmp;
    um(2) = ((a(0) * a(2) + a(1)) * s(0) + (a(1) * a(2) - a(0)) * s(1) + (1.0 + a(2) * a(2)) * s(2)) / tmp;

    RefPartP_suv_m = um + 0.5 * itsReference.getQ()  * itsBunch->getdT() / itsReference.getM() * c * externalE;

}

#endif // OPAL_ParallelTTracker_HH