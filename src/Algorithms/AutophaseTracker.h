#ifndef OPAL_AutophaseTracker_H
#define OPAL_AutophaseTracker_H

//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Algorithms/PartBunch.h"
#include "Utilities/ClassicField.h"
#include "Elements/OpalBeamline.h"
#include "Algorithms/PartPusher.h"

#include "Beamlines/Beamline.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/BeamBeam.h"
#include "AbsBeamline/Collimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Degrader.h"
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

#include "Algorithms/DefaultVisitor.h"

#include <queue>

class BMultipoleField;
class PartData;

class AutophaseTracker: public DefaultVisitor {
public:
    AutophaseTracker(const Beamline &beamline, const PartData &ref, const double &T0);
    virtual ~AutophaseTracker();

    FieldList executeAutoPhaseForSliceTracker();
    void execute(const std::queue<double> &dtAllTracks,
                 const std::queue<double> &maxZ,
                 const std::queue<unsigned long long> &maxTrackSteps);


    virtual void visitBeamline(const Beamline &bl);
    virtual void visitAlignWrapper(const AlignWrapper &);
    virtual void visitBeamBeam(const BeamBeam &);
    virtual void visitCollimator(const Collimator &);
    virtual void visitCorrector(const Corrector &);
    virtual void visitDegrader(const Degrader &);
    virtual void visitDiagnostic(const Diagnostic &);
    virtual void visitDrift(const Drift &);
    virtual void visitLambertson(const Lambertson &);
    virtual void visitMarker(const Marker &);
    virtual void visitMonitor(const Monitor &);
    virtual void visitMultipole(const Multipole &);
    virtual void visitProbe(const Probe &);
    virtual void visitRBend(const RBend &);
    virtual void visitRFCavity(const RFCavity &);
    virtual void visitTravelingWave(const TravelingWave &);
    virtual void visitRFQuadrupole(const RFQuadrupole &);
    virtual void visitSBend(const SBend &);
    virtual void visitSeparator(const Separator &);
    virtual void visitSeptum(const Septum &);
    virtual void visitSolenoid(const Solenoid &);
    virtual void visitParallelPlate(const ParallelPlate &);
    virtual void visitCyclotronValley(const CyclotronValley &);

private:
    ClassicField* checkCavity(double s);
    void updateRFElement(std::string elName, double maxPhi);
    void updateAllRFElements(double phiShift);
    void adjustCavityPhase();
    double APtrack(Component *cavity, double cavity_start, const double &phi) const;
    double getGlobalPhaseShift();
    void handleOverlappingMonitors();
    void prepareSections();
    double getEnergyMeV(const Vector_t &p);

    OpalBeamline itsOpalBeamline_m;
    FieldList cavities_m;
    FieldList::iterator currentAPCavity_m;
    PartBunch itsBunch_m;
    BorisPusher itsPusher_m;
};

inline void AutophaseTracker::visitBeamline(const Beamline &bl) {
    bl.iterate(*dynamic_cast<BeamlineVisitor *>(this), false);
}

inline void AutophaseTracker::visitAlignWrapper(const AlignWrapper &wrap) {
    itsOpalBeamline_m.visit(wrap, *this, &itsBunch_m);
}

inline void AutophaseTracker::visitBeamBeam(const BeamBeam &bb) {
    itsOpalBeamline_m.visit(bb, *this, &itsBunch_m);
}


inline void AutophaseTracker::visitCollimator(const Collimator &coll) {
    itsOpalBeamline_m.visit(coll, *this, &itsBunch_m);
}


inline void AutophaseTracker::visitCorrector(const Corrector &corr) {
    itsOpalBeamline_m.visit(corr, *this, &itsBunch_m);
}


inline void AutophaseTracker::visitDegrader(const Degrader &deg) {
    itsOpalBeamline_m.visit(deg, *this, &itsBunch_m);
}


inline void AutophaseTracker::visitDiagnostic(const Diagnostic &diag) {
    itsOpalBeamline_m.visit(diag, *this, &itsBunch_m);
}


inline void AutophaseTracker::visitDrift(const Drift &drift) {
    itsOpalBeamline_m.visit(drift, *this, &itsBunch_m);
}


inline void AutophaseTracker::visitLambertson(const Lambertson &lamb) {
    itsOpalBeamline_m.visit(lamb, *this, &itsBunch_m);
}


inline void AutophaseTracker::visitMarker(const Marker &marker) {
    itsOpalBeamline_m.visit(marker, *this, &itsBunch_m);
}


inline void AutophaseTracker::visitMonitor(const Monitor &mon) {
    itsOpalBeamline_m.visit(mon, *this, &itsBunch_m);
}


inline void AutophaseTracker::visitMultipole(const Multipole &mult) {
    itsOpalBeamline_m.visit(mult, *this, &itsBunch_m);
}

inline void AutophaseTracker::visitProbe(const Probe &prob) {
    itsOpalBeamline_m.visit(prob, *this, &itsBunch_m);
}


inline void AutophaseTracker::visitRBend(const RBend &bend) {
    itsOpalBeamline_m.visit(bend, *this, &itsBunch_m);
}


inline void AutophaseTracker::visitRFCavity(const RFCavity &as) {
    itsOpalBeamline_m.visit(as, *this, &itsBunch_m);
}

inline void AutophaseTracker::visitTravelingWave(const TravelingWave &as) {
    itsOpalBeamline_m.visit(as, *this, &itsBunch_m);
}


inline void AutophaseTracker::visitRFQuadrupole(const RFQuadrupole &rfq) {
    itsOpalBeamline_m.visit(rfq, *this, &itsBunch_m);
}

inline void AutophaseTracker::visitSBend(const SBend &bend) {
    itsOpalBeamline_m.visit(bend, *this, &itsBunch_m);
}


inline void AutophaseTracker::visitSeparator(const Separator &sep) {
    itsOpalBeamline_m.visit(sep, *this, &itsBunch_m);
}


inline void AutophaseTracker::visitSeptum(const Septum &sept) {
    itsOpalBeamline_m.visit(sept, *this, &itsBunch_m);
}


inline void AutophaseTracker::visitSolenoid(const Solenoid &solenoid) {
    itsOpalBeamline_m.visit(solenoid, *this, &itsBunch_m);
}


inline void AutophaseTracker::visitParallelPlate(const ParallelPlate &pplate) {
    itsOpalBeamline_m.visit(pplate, *this, &itsBunch_m);
}

inline void AutophaseTracker::visitCyclotronValley(const CyclotronValley &cv) {
    itsOpalBeamline_m.visit(cv, *this, &itsBunch_m);
}

inline double AutophaseTracker::getEnergyMeV(const Vector_t &p) {
    return (sqrt(dot(p, p) + 1.0) - 1.0) * itsBunch_m.getM() * 1e-6;
}

#endif