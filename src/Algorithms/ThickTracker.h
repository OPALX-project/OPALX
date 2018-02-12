#ifndef OPAL_ThickTracker_HH
#define OPAL_ThickTracker_HH

// ------------------------------------------------------------------------
// $RCSfile: ThickTracker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ThickTracker
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/Tracker.h"
#include "Steppers/BorisPusher.h"
#include "Structure/DataSink.h"                                                                      

#include "Classic/FixedAlgebra/FTps.h"
#include "Classic/FixedAlgebra/FVps.h"


#include "Algorithms/OrbitThreader.h"
#include "Algorithms/IndexMap.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/BeamBeam.h"
#include "AbsBeamline/Collimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RBend3D.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/ParallelPlate.h"
#include "AbsBeamline/CyclotronValley.h"
                       
#include "Elements/OpalBeamline.h"

#define PSdim 6

class BMultipoleField;

template <class T, unsigned Dim>
class PartBunchBase;

class PlanarArcGeometry;


// Class ThickTracker
// ------------------------------------------------------------------------
/// Track using thick-lens algorithm.
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

class ThickTracker: public Tracker {

public:

    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is initially empty.
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ThickTracker(const Beamline &bl, const PartData &data,
                          bool revBeam, bool revTrack);

    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is taken from [b]bunch[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ThickTracker(const Beamline &bl, PartBunchBase<double, 3> *bunch,
			  DataSink &ds,
                          const PartData &data, 
			  bool revBeam, bool revTrack,
			  const std::vector<unsigned long long> &maxSTEPS,
			  double zstart,
			  const std::vector<double> &zstop,
			  const std::vector<double> &dt,
                          const int& truncOrder);
    
    virtual ~ThickTracker();


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



    /// Apply the algorithm to a beam line.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void visitBeamline(const Beamline &);

    /// Apply the algorithm to the top-level beamline.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void execute();

    void updateReferenceParticle(const BorisPusher &pusher);
    void updateRFElement(std::string elName, double maxPhase);
    void prepareSections();  
    void saveCavityPhases();
    void restoreCavityPhases();
    void changeDT();
    void selectDT();
    void autophaseCavities(const BorisPusher &pusher);
    void findStartPosition(const BorisPusher &pusher);


    typedef FTps<double, PSdim> series_t;

    typedef FVps<double, PSdim> map_t;

    series_t x, px, y, py, delta; //z


    struct structMapTracking {
                	std::string elementName;
                	map_t elementMap;
                	std::size_t nSlices;
                	double elementPos;
                	double stepSize;
                };

    void createHamiltonian(series_t& H, std::shared_ptr<Component> element, double& stepSize, std::size_t& nSlices);
    void fillDrift(std::list<structMapTracking>& mapBeamLine,double& elementPos, double& undefSpace);

    void setHamiltonianDrift(series_t& H, double& beta0, double& gamma0, double& q );
    void setHamiltonianSBend(series_t& H, double& beta0, double& gamma0, double& q, double& h, double& K0 );
    void setHamiltonianRBend(series_t& H, double& beta0, double& gamma0, double& q, double& h, double& K0 );
    void setHamiltonianQuadrupole(series_t& H, double& beta0, double& gamma0, double& q, double& K1 );

    void dumpStats(long long step, bool psDump, bool statDump);
    void writePhaseSpace(const long long step, bool psDump, bool statDump);

private:

    // Not implemented.
    ThickTracker();
    ThickTracker(const ThickTracker &);
    void operator=(const ThickTracker &);

    // Fringe fields for entrance and exit of magnetic elements.
    void applyEntranceFringe(double edge, double curve,
                             const BMultipoleField &field, double scale);
    void applyExitFringe(double edge, double curve,
                         const BMultipoleField &field, double scale);

    Vector_t RefPartR_m; 
    Vector_t RefPartP_m; 

    DataSink *itsDataSink_m; 

    OpalBeamline itsOpalBeamline_m;

    double pathLength_m;

    /// where to start
    double zstart_m;


    /// where to stop
    std::queue<double> zStop_m;

    double dtCurrentTrack_m;  
    std::queue<double> dtAllTracks_m;

    /// The maximal number of steps the system is integrated per TRACK 
    std::queue<unsigned long long> localTrackSteps_m;

    CoordinateSystemTrafo referenceToLabCSTrafo_m;

    bool globalEOL_m;
    
    
    int truncOrder_m;
};

inline void ThickTracker::visitAlignWrapper(const AlignWrapper &wrap) {
    itsOpalBeamline_m.visit(wrap, *this, itsBunch_m);
}

inline void ThickTracker::visitBeamBeam(const BeamBeam &bb) {
    itsOpalBeamline_m.visit(bb, *this, itsBunch_m);
}


inline void ThickTracker::visitCollimator(const Collimator &coll) {
    itsOpalBeamline_m.visit(coll, *this, itsBunch_m);
}


inline void ThickTracker::visitCorrector(const Corrector &corr) {
    itsOpalBeamline_m.visit(corr, *this, itsBunch_m);
}


inline void ThickTracker::visitDegrader(const Degrader &deg) {
    itsOpalBeamline_m.visit(deg, *this, itsBunch_m);
}


inline void ThickTracker::visitDiagnostic(const Diagnostic &diag) {
    itsOpalBeamline_m.visit(diag, *this, itsBunch_m);
}


inline void ThickTracker::visitDrift(const Drift &drift) {
    itsOpalBeamline_m.visit(drift, *this, itsBunch_m);
}


inline void ThickTracker::visitLambertson(const Lambertson &lamb) {
    itsOpalBeamline_m.visit(lamb, *this, itsBunch_m);
}


inline void ThickTracker::visitMarker(const Marker &marker) {
    itsOpalBeamline_m.visit(marker, *this, itsBunch_m);
}


inline void ThickTracker::visitMonitor(const Monitor &mon) {
    itsOpalBeamline_m.visit(mon, *this, itsBunch_m);
}


inline void ThickTracker::visitMultipole(const Multipole &mult) {
    itsOpalBeamline_m.visit(mult, *this, itsBunch_m);
}

inline void ThickTracker::visitProbe(const Probe &prob) {
    itsOpalBeamline_m.visit(prob, *this, itsBunch_m);
}


inline void ThickTracker::visitRBend(const RBend &bend) {
    itsOpalBeamline_m.visit(bend, *this, itsBunch_m);
}


inline void ThickTracker::visitRFCavity(const RFCavity &as) {
    itsOpalBeamline_m.visit(as, *this, itsBunch_m);
}

inline void ThickTracker::visitTravelingWave(const TravelingWave &as) {
    itsOpalBeamline_m.visit(as, *this, itsBunch_m);
}


inline void ThickTracker::visitRFQuadrupole(const RFQuadrupole &rfq) {
    itsOpalBeamline_m.visit(rfq, *this, itsBunch_m);
}

inline void ThickTracker::visitSBend(const SBend &bend) {
    itsOpalBeamline_m.visit(bend, *this, itsBunch_m);
}


inline void ThickTracker::visitSeparator(const Separator &sep) {
    itsOpalBeamline_m.visit(sep, *this, itsBunch_m);
}


inline void ThickTracker::visitSeptum(const Septum &sept) {
    itsOpalBeamline_m.visit(sept, *this, itsBunch_m);
}


inline void ThickTracker::visitSolenoid(const Solenoid &solenoid) {
    itsOpalBeamline_m.visit(solenoid, *this, itsBunch_m);
}


inline void ThickTracker::visitParallelPlate(const ParallelPlate &pplate) {
    itsOpalBeamline_m.visit(pplate, *this, itsBunch_m);
}

inline void ThickTracker::visitCyclotronValley(const CyclotronValley &cv) {
    itsOpalBeamline_m.visit(cv, *this, itsBunch_m);
}



#endif // OPAL_ThickTracker_HH
