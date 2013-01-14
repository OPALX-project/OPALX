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

class BMultipoleField;
class PartBunch;
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
    explicit ThickTracker(const Beamline &bl, const PartBunch &bunch,
                          const PartData &data, bool revBeam, bool revTrack);

    virtual ~ThickTracker();

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

    // Apply the algorithm to a CyclotronValley.
    virtual void visitCyclotronValley(const CyclotronValley &);

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

};

#endif // OPAL_ThickTracker_HH
