#ifndef OPAL_ParallelSliceTracker_HH
#define OPAL_ParallelSliceTracker_HH

// ------------------------------------------------------------------------
// $RCSfile: ParallelSliceTracker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelSliceTracker
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include <vector>
#include <list>

#include "Algorithms/Tracker.h"
#include "Algorithms/PartPusher.h"
#include "Structure/DataSink.h"
#include "Utilities/Options.h"

#include "Ippl.h"

#include "Physics/Physics.h"

#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/BeamBeam.h"
#include "AbsBeamline/Collimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/ParallelPlate.h"
#include "AbsBeamline/CyclotronValley.h"

#include "Algorithms/ParallelTTracker.h"

#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"

class BMultipoleField;
class EnvelopeBunch;
class PlanarArcGeometry;


// Class ParallelSliceTracker
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

class ParallelSliceTracker: public Tracker {

public:

    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is initially empty.
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ParallelSliceTracker(const Beamline &bl, const PartData &data,
                                  bool revBeam, bool revTrack);

    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is taken from [b]bunch[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ParallelSliceTracker(const Beamline &bl, EnvelopeBunch &bunch, DataSink &ds,
                                  const PartData &data, bool revBeam, bool revTrack, int maxSTEPS, double zstop,
                                  ParallelTTracker &mySlApTracker);

    virtual ~ParallelSliceTracker();

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
    virtual void visitMultipole(const Multipole &
                               );
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

    /// Apply the algorithm to a ParallelPlate, it is empty for slice tracker.
    virtual void visitParallelPlate(const ParallelPlate &);

    /// Apply the algorithm to a CyclotronValley, it is empty for slice tracker.
    virtual void visitCyclotronValley(const CyclotronValley &);

    /// Apply the algorithm to the top-level beamline.
    virtual void execute();

    /// Apply the algorithm to a beam line.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void visitBeamline(const Beamline &);

private:

    // Not implemented.
    ParallelSliceTracker();
    ParallelSliceTracker(const ParallelSliceTracker &);
    void operator=(const ParallelSliceTracker &);

    void updateRFElement(string elName, double maxPhi);
    void updateAllRFElements();
    double getCavityPhase(FieldList cav, string name);

    int LastVisited;

    OpalBeamline *itsOpalBeamline_m;

    EnvelopeBunch *itsBunch;

    ParallelTTracker *mySlApTracker_m;

    Vector_t RefPartR_zxy_m;
    Vector_t RefPartP_zxy_m;
    Vector_t RefPartR_suv_m;
    Vector_t RefPartP_suv_m;

    DataSink *itsDataSink;

    /// The maximal number of steps the system is integrated
    long long maxSteps_m;

    double zstop_m;

    /// The scale factor for dimensionless variables
    double scaleFactor_m;
    double SpaceOrientation_m[9];

    // Fringe fields for entrance and exit of magnetic elements.
    void applyEntranceFringe(double edge, double curve,
                             const BMultipoleField &field, double scale);
    void applyExitFringe(double edge, double curve,
                         const BMultipoleField &field, double scale);

    void FieldsOutput(double z, double Ex, double Ey, double Ez,  double Bx, double By, double Bz); //Function to write field output

    void kickParticles();

    void updateReferenceParticle();
    void updateSpaceOrientation(const bool &move = false);
    void kickReferenceParticle(const Vector_t &externalE, const Vector_t &externalB);
    void writePhaseSpace(const long long step, const double &sposRef);
    void writeLastStepPhaseSpace(const long long step, const double &sposRef);

    IpplTimings::TimerRef timeIntegrationTimer1_m;
    IpplTimings::TimerRef timeIntegrationTimer2_m;
    IpplTimings::TimerRef timeFieldEvaluation_m ;
    IpplTimings::TimerRef BinRepartTimer_m;
    IpplTimings::TimerRef WakeFieldTimer_m;

};

inline void ParallelSliceTracker::visitBeamline(const Beamline &bl) {
    itsBeamline_m.iterate(*dynamic_cast<BeamlineVisitor *>(this), false);
}

inline void ParallelSliceTracker::visitAlignWrapper(const AlignWrapper &wrap) {
    itsOpalBeamline_m->visit(wrap, *this, itsBunch);
}

inline void ParallelSliceTracker::visitBeamBeam(const BeamBeam &bb) {
    itsOpalBeamline_m->visit(bb, *this, itsBunch);
}

inline void ParallelSliceTracker::visitCollimator(const Collimator &coll) {
    itsOpalBeamline_m->visit(coll, *this, itsBunch);
}


inline void ParallelSliceTracker::visitCorrector(const Corrector &corr) {
    itsOpalBeamline_m->visit(corr, *this, itsBunch);
}


inline void ParallelSliceTracker::visitDiagnostic(const Diagnostic &diag) {
    itsOpalBeamline_m->visit(diag, *this, itsBunch);
}


inline void ParallelSliceTracker::visitDrift(const Drift &drift) {
    itsOpalBeamline_m->visit(drift, *this, itsBunch);
}


inline void ParallelSliceTracker::visitLambertson(const Lambertson &lamb) {
    itsOpalBeamline_m->visit(lamb, *this, itsBunch);
}


inline void ParallelSliceTracker::visitMarker(const Marker &marker) {
    itsOpalBeamline_m->visit(marker, *this, itsBunch);
}


inline void ParallelSliceTracker::visitMonitor(const Monitor &mon) {
    itsOpalBeamline_m->visit(mon, *this, itsBunch);
}


inline void ParallelSliceTracker::visitMultipole(const Multipole &mult) {
    itsOpalBeamline_m->visit(mult, *this, itsBunch);
}


inline void ParallelSliceTracker::visitProbe(const Probe &prob) {
    itsOpalBeamline_m->visit(prob, *this, itsBunch);
}


inline void ParallelSliceTracker::visitRBend(const RBend &bend) {
    itsOpalBeamline_m->visit(bend, *this, itsBunch);
}


inline void ParallelSliceTracker::visitRFCavity(const RFCavity &as) {
    itsOpalBeamline_m->visit(as, *this, itsBunch);
}

inline void ParallelSliceTracker::visitTravelingWave(const TravelingWave &as) {
    itsOpalBeamline_m->visit(as, *this, itsBunch);
}


inline void ParallelSliceTracker::visitRFQuadrupole(const RFQuadrupole &rfq) {
    itsOpalBeamline_m->visit(rfq, *this, itsBunch);
}

inline void ParallelSliceTracker::visitSBend(const SBend &bend) {
    itsOpalBeamline_m->visit(bend, *this, itsBunch);
}


inline void ParallelSliceTracker::visitSeparator(const Separator &sep) {
    itsOpalBeamline_m->visit(sep, *this, itsBunch);
}


inline void ParallelSliceTracker::visitSeptum(const Septum &sept) {
    itsOpalBeamline_m->visit(sept, *this, itsBunch);
}


inline void ParallelSliceTracker::visitSolenoid(const Solenoid &solenoid) {
    itsOpalBeamline_m->visit(solenoid, *this, itsBunch);
}

inline void ParallelSliceTracker::visitParallelPlate(const ParallelPlate &pplate) {
    //do nothing.
}

inline void ParallelSliceTracker::visitCyclotronValley(const CyclotronValley &cv) {
    // Do nothing here.
}

inline void ParallelSliceTracker::kickParticles() {
    using Physics::c;
}

inline void ParallelSliceTracker::updateReferenceParticle() {
    itsBunch->calcBeamParameters();

}

inline void ParallelSliceTracker::updateSpaceOrientation(const bool &move) {
    itsBunch->calcBeamParameters();
}

inline void ParallelSliceTracker::kickReferenceParticle(const Vector_t &externalE, const Vector_t &externalB) {
    using Physics::c;
}

inline void ParallelSliceTracker::writePhaseSpace(const long long step, const double &sposRef) {
    Inform msg("ParallelSliceTracker");
    Vector_t externalE, externalB;
    Vector_t FDext[6];  // FDext = {BHead, EHead, BRef, ERef, BTail, ETail}.

    // Sample fields at (xmin, ymin, zmin), (xmax, ymax, zmax) and the centroid location. We
    // are sampling the electric and magnetic fields at the back, front and
    // center of the beam.
    Vector_t rmin, rmax;
    itsBunch->get_bounds(rmin, rmax);

    Vector_t pos[3];
    pos[0] = Vector_t(rmax(0), rmax(1), rmax(2));
    pos[1] = Vector_t(0.0, 0.0, sposRef);
    pos[2] = Vector_t(rmin(0), rmin(1), rmin(2));

    for(int k = 0; k < 3; ++k) {
        externalB = Vector_t(0.0);
        externalE = Vector_t(0.0);
        itsOpalBeamline_m->getFieldAt(pos[k], itsBunch->get_rmean(), itsBunch->getT() - 0.5 * itsBunch->getdT(), externalE, externalB);

        FDext[2*k]   = externalB;
        FDext[2*k+1] = externalE * 1e-6;
    }

    if(step % Options::psDumpFreq == 0) {
        // Write fields to .h5 file.
        //itsDataSink->stashPhaseSpaceEnvelope(*itsBunch, FDext, rmax(2), sposRef, rmin(2));
        itsDataSink->writePhaseSpaceEnvelope(*itsBunch, FDext, rmax(2), sposRef, rmin(2));
        msg << "* Wrote beam phase space." << endl;
        msg << *itsBunch << endl;
    }

    if(step % Options::statDumpFreq == 0) {
        // Write statistical data.
        itsDataSink->writeStatData(*itsBunch, FDext, rmax(2), sposRef, rmin(2));
        INFOMSG("* Wrote beam statistics." << endl);
    }
}

inline void ParallelSliceTracker::writeLastStepPhaseSpace(const long long step, const double &sposRef) {
    Inform msg("ParallelSliceTracker");
    if(itsBunch->isValid_m) {
        Vector_t externalE, externalB;
        Vector_t FDext[6];
        Vector_t rmin, rmax;
        itsBunch->get_bounds(rmin, rmax);

        Vector_t pos[3];
        pos[0] = Vector_t(rmax(0), rmax(1), rmax(2));
        pos[1] = Vector_t(0.0, 0.0, sposRef);
        pos[2] = Vector_t(rmin(0), rmin(1), rmin(2));

        for(int k = 0; k < 3; ++k) {
            externalB = Vector_t(0.0);
            externalE = Vector_t(0.0);
            FDext[2*k]   = externalB;
            FDext[2*k+1] = externalE * 1e-6;
        }

        // Write statistical data.
        itsDataSink->writeStatData(*itsBunch, FDext, rmax(2), sposRef, rmin(2));
        INFOMSG("* Wrote beam statistics." << endl);
    } else {
        INFOMSG("* Invalid bunch! No statistics dumped." << endl);
    }
}

#endif // OPAL_ParallelSliceTracker_HH
