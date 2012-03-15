#ifndef OPAL_AutoPhase_HH
#define OPAL_AutoPhase_HH

// ------------------------------------------------------------------------
// $RCSfile: AutoPhase.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AutoPhase
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include <vector>
#include <list>
#include <utility>
#include "Algorithms/ParallelTTracker.h"
#include "Algorithms/PartPusher.h"
#include "Structure/DataSink.h"
#include "Utilities/Options.h"

#include "Physics/Physics.h"
#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"

class BMultipoleField;
class PartBunch;
class PlanarArcGeometry;


// Class AutoPhase
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

typedef pair<double, Vector_t > PhaseEnT;

class AutoPhase: public ParallelTTracker {

public:
    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is initially empty.
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit AutoPhase(const Beamline &bl, const PartData &data,
                       bool revBeam, bool revTrack);

    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is taken from [b]bunch[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit AutoPhase(const Beamline &bl,
                       PartBunch &bunch,
                       DataSink &ds,
                       const PartData &data,
                       bool revBeam,
                       bool revTrack,
                       int maxSTEPS,
                       double zstop,
                       long long actStep,
                       double actT,
                       int numRefinements);

    virtual ~AutoPhase();

    virtual void visitAlignWrapper(const AlignWrapper &);

    /// Apply the algorithm to a BeamBeam.
    virtual void visitBeamBeam(const BeamBeam &);

    /// Apply the algorithm to a collimator.
    virtual void visitCollimator(const Collimator &);

    /// Apply the algorithm to a Corrector.
    virtual void visitCorrector(const Corrector &);

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

    /// Apply the algorithm to the top-level beamline.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void execute();

    /// Apply the algorithm to a beam line.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void visitBeamline(const Beamline &);

private:

    // Not implemented.
    AutoPhase();
    AutoPhase(const AutoPhase &);
    void operator=(const AutoPhase &);


    pair<FieldList::iterator, bool> checkCavity(double s);

    pair<FieldList::iterator, bool> doOneStep(BorisPusher pusher);

    FieldList cavities_m;
    FieldList travelingwaves_m;


    int LastVisited;
    const Beamline &itsBeamline_m;
    OpalBeamline itsOpalBeamline_m;

    PartBunch *itsBunch;
    Vector_t RefPartR_zxy_m;
    Vector_t RefPartP_zxy_m;
    Vector_t RefPartR_suv_m;
    Vector_t RefPartP_suv_m;

    DataSink *itsDataSink;

    /// The number of refinements of the search range of the phase
    int numRefs_m;

    /// The maximal number of steps the system is integrated
    long long maxSteps_m;

    /// where to stop
    double zstop_m;

    long long step_m;

    double actT_m;

    /// The scale factor for dimensionless variables
    double scaleFactor_m;
    double SpaceOrientation_m[9];

    void kickParticles(const BorisPusher &pusher);

    void updateReferenceParticle();
    void updateSpaceOrientation(const bool &move = false);
    Vector_t TransformTo(const Vector_t &vec, const Vector_t &ori) const;
    Vector_t TransformBack(const Vector_t &vec, const Vector_t &ori) const;

    void kickReferenceParticle(const Vector_t &externalE, const Vector_t &externalB);
    void writePhaseSpace(const long long step, const double &sposRef);

    double ptoEMeV(Vector_t p) ;

    IpplTimings::TimerRef timeIntegrationTimer1_m;
    IpplTimings::TimerRef timeIntegrationTimer2_m;
    IpplTimings::TimerRef timeFieldEvaluation_m;

};


inline double AutoPhase::ptoEMeV(Vector_t p) {
    return (sqrt(dot(p, p) + 1.0) - 1.0) * itsBunch->getM() * 1e-6;
}

inline void AutoPhase::visitBeamline(const Beamline &bl) {
    itsBeamline_m.iterate(*dynamic_cast<BeamlineVisitor *>(this), false);
}

inline void AutoPhase::visitAlignWrapper(const AlignWrapper &wrap) {
    itsOpalBeamline_m.visit(wrap, *this, itsBunch);
}

inline void AutoPhase::visitBeamBeam(const BeamBeam &bb) {
    itsOpalBeamline_m.visit(bb, *this, itsBunch);
}


inline void AutoPhase::visitCollimator(const Collimator &coll) {
    itsOpalBeamline_m.visit(coll, *this, itsBunch);
}


inline void AutoPhase::visitCorrector(const Corrector &corr) {
    itsOpalBeamline_m.visit(corr, *this, itsBunch);
}


inline void AutoPhase::visitDiagnostic(const Diagnostic &diag) {
    itsOpalBeamline_m.visit(diag, *this, itsBunch);
}


inline void AutoPhase::visitDrift(const Drift &drift) {
    itsOpalBeamline_m.visit(drift, *this, itsBunch);
}


inline void AutoPhase::visitLambertson(const Lambertson &lamb) {
    itsOpalBeamline_m.visit(lamb, *this, itsBunch);
}


inline void AutoPhase::visitMarker(const Marker &marker) {
    itsOpalBeamline_m.visit(marker, *this, itsBunch);
}


inline void AutoPhase::visitMonitor(const Monitor &mon) {
    itsOpalBeamline_m.visit(mon, *this, itsBunch);
}


inline void AutoPhase::visitMultipole(const Multipole &mult) {
    itsOpalBeamline_m.visit(mult, *this, itsBunch);
}

inline void AutoPhase::visitProbe(const Probe &prob) {
    itsOpalBeamline_m.visit(prob, *this, itsBunch);
}


inline void AutoPhase::visitRBend(const RBend &bend) {
    itsOpalBeamline_m.visit(bend, *this, itsBunch);
}


inline void AutoPhase::visitRFCavity(const RFCavity &as) {
    itsOpalBeamline_m.visit(as, *this, itsBunch);
}

inline void AutoPhase::visitTravelingWave(const TravelingWave &as) {
    itsOpalBeamline_m.visit(as, *this, itsBunch);
}


inline void AutoPhase::visitRFQuadrupole(const RFQuadrupole &rfq) {
    itsOpalBeamline_m.visit(rfq, *this, itsBunch);
}

inline void AutoPhase::visitSBend(const SBend &bend) {
    itsOpalBeamline_m.visit(bend, *this, itsBunch);
}


inline void AutoPhase::visitSeparator(const Separator &sep) {
    itsOpalBeamline_m.visit(sep, *this, itsBunch);
}


inline void AutoPhase::visitSeptum(const Septum &sept) {
    itsOpalBeamline_m.visit(sept, *this, itsBunch);
}


inline void AutoPhase::visitSolenoid(const Solenoid &solenoid) {
    itsOpalBeamline_m.visit(solenoid, *this, itsBunch);
}

inline void AutoPhase::kickParticles(const BorisPusher &pusher) {
    using Physics::c;

    int localNum = itsBunch->getLocalNum();
    for(int i = 0; i < localNum; ++i)
        pusher.kick(itsBunch->R[i], itsBunch->P[i], itsBunch->Ef[i], itsBunch->Bf[i], itsBunch->dt[i]);
    itsBunch->calcBeamParameters();
}

inline void AutoPhase::updateReferenceParticle() {
    RefPartR_suv_m = itsBunch->get_rmean() * scaleFactor_m;
    RefPartP_suv_m = itsBunch->get_pmean();

    /* Update the position of the reference particle in ZXY-coordinates. The angle between the ZXY- and the SUV-coordinate
     *  system is determined by the momentum of the reference particle. We calculate the momentum of the reference
     *  particle by rotating the centroid momentum (= momentum of the reference particle in the SUV-coordinate system).
     *  Then we push the reference particle with this momentum for half a time step.
     */

    double gamma = sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));

    /* First update the momentum of the reference particle in zxy coordinate system, then update its position     */
    RefPartP_zxy_m(0) = SpaceOrientation_m[0] * RefPartP_suv_m(0) + SpaceOrientation_m[1] * RefPartP_suv_m(1) + SpaceOrientation_m[2] * RefPartP_suv_m(2);
    RefPartP_zxy_m(1) = SpaceOrientation_m[3] * RefPartP_suv_m(0) + SpaceOrientation_m[4] * RefPartP_suv_m(1) + SpaceOrientation_m[5] * RefPartP_suv_m(2);
    RefPartP_zxy_m(2) = SpaceOrientation_m[6] * RefPartP_suv_m(0) + SpaceOrientation_m[7] * RefPartP_suv_m(1) + SpaceOrientation_m[8] * RefPartP_suv_m(2);

    RefPartR_zxy_m += RefPartP_zxy_m * scaleFactor_m / (2. * gamma);
    RefPartR_suv_m += RefPartP_suv_m * scaleFactor_m / (2. * gamma);

}

inline void AutoPhase::updateSpaceOrientation(const bool &move) {
    /* Update the position of the reference particle in ZXY-coordinates. The angle between the ZXY- and the SUV-coordinate
     *  system is determined by the momentum of the reference particle. We calculate the momentum of the reference
     *  particle by rotating the centroid momentum (= momentum of the reference particle in the SUV-coordinate system).
     *  Then we push the reference particle with this momentum for half a time step.
     */
    double AbsMomentum = sqrt(dot(RefPartP_suv_m, RefPartP_suv_m));
    double AbsMomentumProj = sqrt(RefPartP_suv_m(0) * RefPartP_suv_m(0) + RefPartP_suv_m(2) * RefPartP_suv_m(2));

    SpaceOrientation_m[0] = RefPartP_zxy_m(2) / AbsMomentumProj;
    SpaceOrientation_m[1] = -RefPartP_zxy_m(0) * RefPartP_zxy_m(1) / (AbsMomentum * AbsMomentumProj);
    SpaceOrientation_m[2] = RefPartP_zxy_m(0) / AbsMomentum;
    SpaceOrientation_m[3] = 0.;
    SpaceOrientation_m[4] = AbsMomentumProj / AbsMomentum;
    SpaceOrientation_m[5] = RefPartP_zxy_m(1) / AbsMomentum;
    SpaceOrientation_m[6] = -RefPartP_zxy_m(0) / AbsMomentumProj;
    SpaceOrientation_m[7] = -RefPartP_zxy_m(2) * RefPartP_zxy_m(1) / (AbsMomentum * AbsMomentumProj);
    SpaceOrientation_m[8] = RefPartP_zxy_m(2) / AbsMomentum;

    Vector_t EulerAngles;
    if(RefPartP_suv_m(2) < 1.0e-12) {
        EulerAngles = Vector_t(0.0);
    } else {
        EulerAngles = Vector_t(-atan(RefPartP_suv_m(0) / RefPartP_suv_m(2)),          \
                               -asin(RefPartP_suv_m(1) / AbsMomentum),    \
                               0.0);
    }
    // rotate the local coordinate system of all sections which are online
    Vector_t smin(0.0), smax(0.0);
    itsBunch->get_bounds(smin, smax);
    itsOpalBeamline_m.updateOrientation(EulerAngles, RefPartR_suv_m * scaleFactor_m, smin(2)*scaleFactor_m, smax(2)*scaleFactor_m);
    itsBunch->rotateAbout(RefPartR_suv_m, RefPartP_suv_m);
    if(move)       // move the bunch such that the new centroid location is at (0,0,z)
        itsBunch->moveBy(Vector_t(-RefPartR_suv_m(0), -RefPartR_suv_m(1), 0.0));
    itsBunch->calcBeamParameters();
}

inline Vector_t AutoPhase::TransformTo(const Vector_t &vec, const Vector_t &ori) const {
    const double sina = sin(ori(0));
    const double cosa = cos(ori(0));
    const double sinb = sin(ori(1));
    const double cosb = cos(ori(1));

    Vector_t temp(0.0, 0.0, 0.0);

    temp(0) =  cosa *        vec(0)                 - sina *        vec(2);
    temp(1) = -sina * sinb * vec(0) + cosb * vec(1) - cosa * sinb * vec(2);
    temp(2) =  sina * cosb * vec(0) + sinb * vec(1) + cosa * cosb * vec(2);

    return temp;
}

inline Vector_t AutoPhase::TransformBack(const Vector_t &vec, const Vector_t &ori) const {
    const double sina = sin(ori(0));
    const double cosa = cos(ori(0));
    const double sinb = sin(ori(1));
    const double cosb = cos(ori(1));

    Vector_t temp(0.0, 0.0, 0.0);

    temp(0) =  cosa * vec(0) - sina * sinb * vec(1) + sina * cosb * vec(2);
    temp(1) =                  cosb *        vec(1) + sinb *        vec(2);
    temp(2) = -sina * vec(0) + cosa * sinb * vec(1) + cosa * cosb * vec(2);

    return temp;
}

inline void AutoPhase::kickReferenceParticle(const Vector_t &externalE, const Vector_t &externalB) {
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


inline void AutoPhase::writePhaseSpace(const long long step, const double &sposRef) {
    Vector_t externalE, externalB;
    Vector_t FDext[6];  // FDext = {BHead, EHead, BRef, ERef, BTail, ETail}.

    // Sample fields at (xmin, ymin, zmin), (xmax, ymax, zmax) and the centroid location. We
    // are sampling the electric and magnetic fields at the back, front and
    // center of the beam.

    Vector_t pos[3];
    pos[0] = Vector_t(0.0, 0.0, sposRef);
    pos[1] = Vector_t(0.0, 0.0, sposRef);
    pos[2] = Vector_t(0.0, 0.0, sposRef);

    for(int k = 0; k < 3; ++k) {
        externalB = Vector_t(0.0);
        externalE = Vector_t(0.0);
        itsOpalBeamline_m.getFieldAt(pos[k], itsBunch->get_rmean(), itsBunch->getT() - 0.5 * itsBunch->getdT(), externalE, externalB);

        FDext[2*k]   = externalB;
        FDext[2*k+1] = externalE;
    }

    if(step % Options::statDumpFreq == 0) {
        itsDataSink->writeStatData(*itsBunch, FDext, sposRef, sposRef, sposRef);
    }
}


/*
namespace std{
    std::ostream& operator<<(std::ostream& stream, const PhaseEnT& sample)
    {
        return stream << "Phase: " << sample.first << ", Energy: " << sample.second;
    }
}
*/


#endif // OPAL_AutoPhase_HH
