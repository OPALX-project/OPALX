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
#include "Structure/DataSink.h"
#include "Utilities/Options.h"
#include <vector>
#include <list>

class BMultipoleField;
class PartBunch;
class PlanarArcGeometry;


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

class ParallelTTracker: public Tracker {

public:
struct FieldListType1Entry
  {
    Component* Element;
    double Start;
    double End;
    FieldListType1Entry(Component* aComp, double start, double end)
    {
      Element = aComp;
      Start = start;
      End = end;
    }
    static bool SortAscByStart(FieldListType1Entry e1, FieldListType1Entry e2)
    {
      return (e1.Start < e2.Start);
    }
    
    static bool SortAscByEnd(FieldListType1Entry e1, FieldListType1Entry e2)
    {
      return (e1.End < e2.End);
    }
    static bool SortDescByStart(FieldListType1Entry e1, FieldListType1Entry e2)
    {
      return (e1.Start > e2.Start);
    }
    
    static bool SortDescByEnd(FieldListType1Entry e1, FieldListType1Entry e2)
    {
      return (e1.End > e2.End);
    }
    
  };

  typedef list<FieldListType1Entry> FieldListType1;
  typedef list<FieldListType1Entry>::iterator FieldListType1Iterator;

  struct FieldListType2Entry
  {
    list<Component*> Elements;
    double Start;
    double End;
    FieldListType2Entry(list<Component*> aList, double start, double end)
    {
      Elements.assign(aList.begin(), aList.end());
      Start = start;
      End = end;
    }
    static bool SortAscByStart(FieldListType2Entry e1, FieldListType2Entry e2)
    {
      return (e1.Start < e2.Start);
    }

    static bool SortAscByEnd(FieldListType2Entry e1, FieldListType2Entry e2)
    {
      return (e1.End < e2.End);
    }
    static bool SortDescByStart(FieldListType2Entry e1, FieldListType2Entry e2)
    {
      return (e1.Start > e2.Start);
    }

    static bool SortDescByEnd(FieldListType2Entry e1, FieldListType2Entry e2)
    {
      return (e1.End > e2.End);
    }
  };

  typedef vector<FieldListType2Entry> FieldListType2;
  typedef vector<FieldListType2Entry>::iterator FieldListType2Iterator;
  typedef list<Component*>::iterator ElementIterator;




  /// Constructor.
  //  The beam line to be tracked is "bl".
  //  The particle reference data are taken from "data".
  //  The particle bunch tracked is initially empty.
  //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
  //  If [b]revTrack[/b] is true, we track against the beam.
  explicit ParallelTTracker(const Beamline &bl, const PartData &data,
		       bool revBeam, bool revTrack);

  /// Constructor.
  //  The beam line to be tracked is "bl".
  //  The particle reference data are taken from "data".
  //  The particle bunch tracked is taken from [b]bunch[/b].
  //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
  //  If [b]revTrack[/b] is true, we track against the beam.
  explicit ParallelTTracker(const Beamline &bl, PartBunch &bunch, DataSink &ds,
                        const PartData &data, bool revBeam, bool revTrack, int maxSTEPS);

  virtual ~ParallelTTracker();

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
  ParallelTTracker();
  ParallelTTracker(const ParallelTTracker &);
  void operator=(const ParallelTTracker &);

  FieldListType1 myElements;
  FieldListType2 myFieldList;

  int LastVisited;
  Beamline *itsBeamline;

  PartBunch *itsBunch;

  DataSink *itsDataSink;

  /// The maximal number of steps the system is integrated
  int maxSteps_m;

  /// The scale factor for dimensionless variables
  double scaleFactor_m;

  // Fringe fields for entrance and exit of magnetic elements.
  void applyEntranceFringe(double edge, double curve,
			   const BMultipoleField &field, double scale);
  void applyExitFringe(double edge, double curve,
		       const BMultipoleField &field, double scale);

  /// Build up a 2D map of elements
  //
  void buildupFieldList();

  IpplTimings::TimerRef timeIntegrationTimer1_m;
  IpplTimings::TimerRef timeIntegrationTimer2_m;
  IpplTimings::TimerRef timeFieldEvaluation_m ;
  IpplTimings::TimerRef BinRepartTimer_m;
  IpplTimings::TimerRef WakeFieldTimer_m;

};

#endif // OPAL_ParallelTTracker_HH
