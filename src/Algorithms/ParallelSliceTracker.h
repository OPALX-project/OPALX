#ifdef HAVE_ENVELOPE_SOLVER
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

#include "Algorithms/Tracker.h"
#include "Algorithms/bet/SLDataSink.h"
#include "Utilities/Options.h"
#include <vector>
#include <list>

#include "Physics/Physics.h"

#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/Collimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"

class BMultipoleField;
class SLPartBunch;
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
  struct FieldListEntry
  {
    Component* Element;
    double Start;
    double End;
    bool goneLive;
    bool goneOff;
    FieldListEntry(Component* aComp, double start, double end)
    {
      Element = aComp;
      Start = start;
      End = end;
      goneLive = false;
      goneOff = false;
    }
    static bool SortAscByStart(FieldListEntry e1, FieldListEntry e2)
    {
      return (e1.Start < e2.Start);
    }
    
    static bool SortAscByEnd(FieldListEntry e1, FieldListEntry e2)
    {
      return (e1.End < e2.End);
    }
    static bool SortDescByStart(FieldListEntry e1, FieldListEntry e2)
    {
      return (e1.Start > e2.Start);
    }
    
    static bool SortDescByEnd(FieldListEntry e1, FieldListEntry e2)
    {
      return (e1.End > e2.End);
    }
    
    static bool ZeroLength(FieldListEntry e1)
    {
      return e1.Start == e1.End;
    }
  };

  typedef list<FieldListEntry> FieldList;
  typedef list<FieldListEntry>::iterator FieldListIterator;

  typedef list<Component*>::iterator ElementIterator;

  struct SectionListEntry
  {
    list<Component*> Elements;
    double Start;
    double End;
    bool bends;
    bool hasWake;
    SectionListEntry(list<Component*> aList, double start, double end)
    {
      Elements.assign(aList.begin(), aList.end());
      Start = start;
      End = end;
      bends = false;
      hasWake = false;
      for (ElementIterator el_it = aList.begin(); el_it != aList.end(); ++el_it) {
        bends = bends || (*el_it)->bends();
        if ((*el_it)->bends())
          *gmsg << (*el_it)->getName() << " bends" << endl;
        hasWake = hasWake || (*el_it)->hasWake();
        if ((*el_it)->hasWake())
          *gmsg << (*el_it)->getName() << " has a wake function attached" << endl;
      }
      if (hasWake)
        *gmsg << "in section from " << Start << " to " << End << " an element has a wake attached" << endl;
    }
    static bool SortAscByStart(SectionListEntry e1, SectionListEntry e2)
    {
      return (e1.Start < e2.Start);
    }

    static bool SortAscByEnd(SectionListEntry e1, SectionListEntry e2)
    {
      return (e1.End < e2.End);
    }
    static bool SortDescByStart(SectionListEntry e1, SectionListEntry e2)
    {
      return (e1.Start > e2.Start);
    }

    static bool SortDescByEnd(SectionListEntry e1, SectionListEntry e2)
    {
      return (e1.End > e2.End);
    }
  };

  typedef vector<SectionListEntry> SectionList;
  typedef vector<SectionListEntry>::iterator SectionListIterator;
  




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
  explicit ParallelSliceTracker(const Beamline &bl, SLPartBunch &bunch, SLDataSink &ds,
                            const PartData &data, bool revBeam, bool revTrack, int maxSTEPS);

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
  ParallelSliceTracker();
  ParallelSliceTracker(const ParallelSliceTracker &);
  void operator=(const ParallelSliceTracker &);

  FieldList allElements;
  SectionList allSections;

  int LastVisited;
  Beamline *itsBeamline;

  SLPartBunch *itsBunch;
  Vector_t RefPartR_zxy_m;
  Vector_t RefPartP_zxy_m;
  Vector_t RefPartR_suv_m;
  Vector_t RefPartP_suv_m;

  SLDataSink *itsDataSink;

  /// The maximal number of steps the system is integrated
  long long maxSteps_m;

  /// The scale factor for dimensionless variables
  double scaleFactor_m;
  double SpaceOrientation_m[9];

  // Fringe fields for entrance and exit of magnetic elements.
  void applyEntranceFringe(double edge, double curve,
                           const BMultipoleField &field, double scale);
  void applyExitFringe(double edge, double curve,
                       const BMultipoleField &field, double scale);

  /// Build up a 2D map of elements
  //
  void buildupFieldList();
  
  void kickParticles();
  
  void updateReferenceParticle();
  void updateSpaceOrientation(bool move = false);
  void kickReferenceParticle(const Vector_t &externalE, const Vector_t &externalB);
  bool getExternalField(const Vector_t &R, const double &t, Vector_t &extE, Vector_t &extB);
  void writePhaseSpace(const double &sposRef);

  IpplTimings::TimerRef timeIntegrationTimer1_m;
  IpplTimings::TimerRef timeIntegrationTimer2_m;
  IpplTimings::TimerRef timeFieldEvaluation_m ;
  IpplTimings::TimerRef BinRepartTimer_m;
  IpplTimings::TimerRef WakeFieldTimer_m;

};

inline void ParallelSliceTracker::kickParticles()
{
  using Physics::c;


}

inline void ParallelSliceTracker::updateReferenceParticle()
{
  itsBunch->calcBeamParameters();

}

inline void ParallelSliceTracker::updateSpaceOrientation(bool move)
{
  itsBunch->calcBeamParameters();    
}

inline void ParallelSliceTracker::kickReferenceParticle(const Vector_t &externalE, const Vector_t &externalB)
{
  using Physics::c;

}

inline bool ParallelSliceTracker::getExternalField(const Vector_t &R, const double &t, Vector_t &extE, Vector_t &extB)
{
  bool bends = false;

  extE = Vector_t(0,0,0);
  extB = Vector_t(0,0,0);

  ElementIterator element_it;
  SectionListIterator sections_it;
  SectionListIterator sections_end_it = allSections.end();

  if ( R(2) <= allSections.back().End 
       && 
       R(2) >= allSections.front().Start )
    {
      for (sections_it = allSections.begin(); 
           sections_it != sections_end_it;
           ++sections_it)
        {
          if (R(2) >= (*sections_it).Start 
              && 
              R(2) <  (*sections_it).End)
            {
              for (element_it = (*sections_it).Elements.begin(); 
                   element_it != (*sections_it).Elements.end(); 
                   ++element_it)
                {
                  if (!(*element_it)->Online())
                    {
                      *gmsg << "* ************** W A R N I N G *****************************************************" << endl;
                      *gmsg << "* trying to get field from element " << (*element_it)->getName() << " which is not online yet" << endl;
                      *gmsg << "* **********************************************************************************" << endl;
                    }
                  else
                    {
                      bends = bends || (*element_it)->bends(); // is the reference particle in a bending field?
                      (*element_it)->apply(R, t, extE, extB);
                    }
                }
              break;
            }
        }
    }
  return bends;

}
inline void ParallelSliceTracker::writePhaseSpace(const double &sposRef)
{
  Vector_t externalE, externalB;
  //FDext = {BHead, EHead, BRef, ERef, BTail, ETail}
  Vector_t FDext[6];
  
  //sample fields at (0,0,rmin), (0,0,rmax) and the centroid location
  Vector_t rmin, rmax;
  itsBunch->get_bounds(rmin,rmax);
  
  Vector_t pos[3];
  pos[0] = Vector_t(rmax(0),rmax(1),rmax(2));
  pos[1] = Vector_t(0.00,0.00,sposRef);
  pos[2] = Vector_t(rmin(0),rmin(1),rmin(2));
  
  for(int k=0; k < 3; ++k) 
    {
      getExternalField(pos[k], itsBunch->getT() - 0.5 * itsBunch->getdT(), externalE, externalB);

      FDext[2*k]   = externalB;
      FDext[2*k+1] = externalE;
    }
  
  itsDataSink->writePhaseSpace(*itsBunch);
  *gmsg << *itsBunch << endl;

}

#endif // OPAL_ParallelSliceTracker_HH
#endif
