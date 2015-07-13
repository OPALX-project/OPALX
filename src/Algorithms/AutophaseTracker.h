#ifndef OPAL_AutophaseTracker_H
#define OPAL_AutophaseTracker_H

//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Algorithms/PartBunch.h"
#include "Algorithms/PartPusher.h"
#include "Algorithms/PartData.h"
#include "Algorithms/PBunchDefs.h"

class ClassicField;
#include "Elements/OpalBeamline.h"

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

#define AP_VISITELEMENT(elem) virtual void visit##elem(const elem &el) \
    { itsOpalBeamline_m.visit(el, *this, &itsBunch_m); }

#define AP_IGNOREELEMENT(elem) virtual void visit##elem(const elem &) { }

class AutophaseTracker: public DefaultVisitor {
public:
    AutophaseTracker(const Beamline &beamline,
                     const PartData &ref,
                     const double &T0,
                     double initialR,
                     double initialP);
    virtual ~AutophaseTracker();

    void execute(const std::queue<double> &dtAllTracks,
                 const std::queue<double> &maxZ,
                 const std::queue<unsigned long long> &maxTrackSteps);

    void save(const std::string &fname);


    virtual void visitBeamline(const Beamline &bl);
    AP_VISITELEMENT(AlignWrapper)
    AP_IGNOREELEMENT(BeamBeam)
    AP_IGNOREELEMENT(Collimator)
    AP_IGNOREELEMENT(Corrector)
    AP_IGNOREELEMENT(CyclotronValley)
    AP_IGNOREELEMENT(Degrader)
    AP_IGNOREELEMENT(Diagnostic)
    AP_IGNOREELEMENT(Drift)
    AP_IGNOREELEMENT(Lambertson)
    AP_IGNOREELEMENT(Marker)
    AP_IGNOREELEMENT(Monitor)
    AP_IGNOREELEMENT(Multipole)
    AP_IGNOREELEMENT(ParallelPlate)
    AP_IGNOREELEMENT(Probe)
    AP_IGNOREELEMENT(RBend)
    AP_VISITELEMENT(RFCavity)
    AP_IGNOREELEMENT(RFQuadrupole)
    AP_IGNOREELEMENT(SBend)
    AP_IGNOREELEMENT(Separator)
    AP_IGNOREELEMENT(Septum)
    AP_IGNOREELEMENT(Solenoid)
    AP_VISITELEMENT(TravelingWave)

private:
    double APtrack(const std::shared_ptr<Component> &cavity, double cavity_start, const double &phi) const;
    void track(double uptoSPos,
               size_t &step,
               std::queue<double> &dtAllTracks,
               std::queue<double> &maxZ,
               std::queue<unsigned long long> &maxTrackSteps);

    double getEnergyMeV(const Vector_t &p);
    void evaluateField();
    std::shared_ptr<Component> getNextCavity(const std::shared_ptr<Component> &current);
    void advanceTime(size_t & step, const size_t maxSteps, const double beginNextCavity);
    double guessCavityPhase(const std::shared_ptr<Component> &);
    double optimizeCavityPhase(const std::shared_ptr<Component> &,
                               double initialPhase,
                               size_t currentStep,
                               double dt);
    double getBeginCavity(const std::shared_ptr<Component> &);
    double getEndCavity(const std::shared_ptr<Component> &);
    void sendCavityPhases();
    void receiveCavityPhases();
    void printCavityPhases();

    OpalBeamline itsOpalBeamline_m;
    PartBunch itsBunch_m;
    BorisPusher itsPusher_m;
    Layout_t *particleLayout_m;

    IpplTimings::TimerRef timeIntegrationTimer_m;
    IpplTimings::TimerRef fieldEvaluationTimer_m;

};

inline
void AutophaseTracker::visitBeamline(const Beamline &bl) {
    bl.iterate(*dynamic_cast<BeamlineVisitor *>(this), false);
}

inline
double AutophaseTracker::getEnergyMeV(const Vector_t &p) {
    return (sqrt(dot(p, p) + 1.0) - 1.0) * itsBunch_m.getM() * 1e-6;
}

inline
double AutophaseTracker::getBeginCavity(const std::shared_ptr<Component> &comp) {
    if (comp == NULL) return -1e6;

    double startComp = 0.0, endComp = 0.0;
    comp->getDimensions(startComp, endComp);

    return startComp;
}

inline
double AutophaseTracker::getEndCavity(const std::shared_ptr<Component> &comp) {
    if (comp == NULL) return -1e6;

    double startComp = 0.0, endComp = 0.0;
    comp->getDimensions(startComp, endComp);

    return endComp;
}

#endif