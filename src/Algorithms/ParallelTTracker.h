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

#include "Solvers/WakeFunction.hh"
#include "Structure/Wake.h"

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
    bool live;
    WakeFunction *Wake;
    Vector_t Orientation;
    
    SectionListEntry(list<Component*> aList, double start, double end)
    {
      Elements.assign(aList.begin(), aList.end());
      Start = start;
      End = end;
      bends = false;
      hasWake = false;
      live = false;
      Orientation = Vector_t(0.0);
      for (ElementIterator el_it = aList.begin(); el_it != aList.end(); ++el_it) {
        if ((*el_it)->bends())
          {
            bends = true;
/*             *gmsg << (*el_it)->getName() << " bends" << endl; */
            Orientation = (*el_it)->getOrientation();
          }

        if (!hasWake && (*el_it)->hasWake())
          {
            hasWake = true;
            Wake = (*el_it)->getWake()->wf_m;
            *gmsg << (*el_it)->getName() << " has a wake function attached " << endl;
          }
      }
    }
    double getStart(double u, double v)
    {
      return Start - (sin(Orientation(0))*u + tan(Orientation(1))*v) / cos(Orientation(0)); // ???? is this the correct formula????
    }

    double getEnd(double u, double v)
    {
      return End - (sin(Orientation(0))*u + tan(Orientation(1))*v) / cos(Orientation(0)); // ???? is this the correct formula????
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
  ParallelTTracker();
  ParallelTTracker(const ParallelTTracker &);
  void operator=(const ParallelTTracker &);

  FieldList allElements;
  SectionList allSections;

  int LastVisited;
  Beamline *itsBeamline;

  PartBunch *itsBunch;
  Vector_t RefPartR_zxy_m;
  Vector_t RefPartP_zxy_m;
  Vector_t RefPartR_suv_m;
  Vector_t RefPartP_suv_m;

  DataSink *itsDataSink;

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
  
  void kickParticles(const BorisPusher &pusher);
  
  void updateReferenceParticle();
  void updateSpaceOrientation(bool move = false);
  Vector_t TransformTo(const Vector_t &vec, const Vector_t &ori) const;
  Vector_t TransformBack(const Vector_t &vec, const Vector_t &ori) const;

  void kickReferenceParticle(const Vector_t &externalE, const Vector_t &externalB);
  bool getExternalField(const Vector_t &R, const double &t, Vector_t &extE, Vector_t &extB);
  void writePhaseSpace(const double &sposRef);

  IpplTimings::TimerRef timeIntegrationTimer1_m;
  IpplTimings::TimerRef timeIntegrationTimer2_m;
  IpplTimings::TimerRef timeFieldEvaluation_m ;
  IpplTimings::TimerRef BinRepartTimer_m;
  IpplTimings::TimerRef WakeFieldTimer_m;

};

inline void ParallelTTracker::kickParticles(const BorisPusher &pusher)
{
  using Physics::c;

  int localNum = itsBunch->getLocalNum();
  for (int i = 0; i < localNum; ++i) 
    pusher.kick(itsBunch->R[i], itsBunch->P[i], itsBunch->Ef[i], itsBunch->Bf[i], itsBunch->dt[i]);
}

inline void ParallelTTracker::updateReferenceParticle()
{
  itsBunch->calcBeamParameters();
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

inline void ParallelTTracker::updateSpaceOrientation(bool move)
{
  /* Update the position of the reference particle in ZXY-coordinates. The angle between the ZXY- and the SUV-coordinate
   *  system is determined by the momentum of the reference particle. We calculate the momentum of the reference
   *  particle by rotating the centroid momentum (= momentum of the reference particle in the SUV-coordinate system).
   *  Then we push the reference particle with this momentum for half a time step.
   */
  double AbsMomentum = sqrt(dot(RefPartP_zxy_m,RefPartP_zxy_m));
  double AbsMomentumProj = sqrt(RefPartP_zxy_m(0) * RefPartP_zxy_m(0) + RefPartP_zxy_m(2) * RefPartP_zxy_m(2));
  
  SpaceOrientation_m[0] = RefPartP_zxy_m(2) / AbsMomentumProj;
  SpaceOrientation_m[1] = -RefPartP_zxy_m(0) * RefPartP_zxy_m(1)/(AbsMomentum * AbsMomentumProj);
  SpaceOrientation_m[2] = RefPartP_zxy_m(0) / AbsMomentum;
  SpaceOrientation_m[3] = 0.;
  SpaceOrientation_m[4] = AbsMomentumProj / AbsMomentum;
  SpaceOrientation_m[5] = RefPartP_zxy_m(1) / AbsMomentum;
  SpaceOrientation_m[6] = -RefPartP_zxy_m(0) / AbsMomentumProj;
  SpaceOrientation_m[7] = -RefPartP_zxy_m(2) * RefPartP_zxy_m(1) / (AbsMomentum * AbsMomentumProj);
  SpaceOrientation_m[8] = RefPartP_zxy_m(2) / AbsMomentum;

  Vector_t EulerAngles  = Vector_t(-atan(RefPartP_suv_m(0)/RefPartP_suv_m(2)), \
                                   -asin(RefPartP_suv_m(1)/AbsMomentum), \
                                   0.0);
  
  // rotate the local coordinate system of all sections which are online
  for (SectionListIterator Sections_it = allSections.begin(); Sections_it != allSections.end(); ++Sections_it)
    if ((*Sections_it).live)
      {
        (*Sections_it).Orientation += EulerAngles;
      }

  itsBunch->rotateAbout(RefPartR_suv_m,EulerAngles);   
  if (move)      // move the bunch such that the new centroid location is at (0,0,z)
      itsBunch->moveBy(Vector_t(-RefPartR_suv_m(0),-RefPartR_suv_m(1),0.0));
  itsBunch->calcBeamParameters();
}

inline Vector_t ParallelTTracker::TransformTo(const Vector_t &vec, const Vector_t &ori) const
{
  double sina = sin(ori(0));
  double cosa = cos(ori(0));
  double sinb = sin(ori(1));
  double cosb = cos(ori(1));

  Vector_t temp(0.0, 0.0, 0.0);
  
/*   temp(0) =  cosa * vec(0)         - sina * vec(2); */
/*   temp(1) =                vec(1); */
/*   temp(2) =  sina * vec(0)         + cosa * vec(2); */
  temp(0) =  cosa *        vec(0)                 - sina *        vec(2);
  temp(1) = -sina * sinb * vec(0) + cosb * vec(1) - cosa * sinb * vec(2);
  temp(2) =  sina * cosb * vec(0) + sinb * vec(1) + cosa * cosb * vec(2);

  return temp;
}

inline Vector_t ParallelTTracker::TransformBack(const Vector_t &vec, const Vector_t &ori) const
{
  double sina = sin(ori(0));
  double cosa = cos(ori(0));
  double sinb = sin(ori(1));
  double cosb = cos(ori(1));

  Vector_t temp(0.0, 0.0, 0.0);
  
/*   temp(0) =  cosa * vec(0)          + sina * vec(2); */
/*   temp(1) =                  vec(1); */
/*   temp(2) = -sina * vec(0)          + cosa * vec(2); */
  temp(0) =  cosa * vec(0) - sina * sinb * vec(1) + sina * cosb * vec(2);
  temp(1) =                  cosb *        vec(1) + sinb *        vec(2);
  temp(2) = -sina * vec(0) + cosa * sinb * vec(1) + cosa * cosb * vec(2);

  return temp;
}

inline void ParallelTTracker::kickReferenceParticle(const Vector_t &externalE, const Vector_t &externalB)
{
  using Physics::c;

  // track reference particle
  Vector_t um = RefPartP_suv_m + 0.5 * itsReference.getQ() * itsBunch->getdT() / itsReference.getM() * c * externalE;
  double gamma = sqrt(1.0 + dot(um,um));

  double tmp = 0.5 * itsReference.getQ() * c*c * itsBunch->getdT() / (itsReference.getM() * gamma);
  Vector_t a = tmp * externalB;

  Vector_t s = um + tmp * cross(um,externalB);

  tmp = 1.0 + dot(a,a);
        
  um(0) = ((1.0 + a(0)*a(0))    * s(0) + (a(0) * a(1) + a(2)) * s(1) + (a(0) * a(2) - a(1)) * s(2)) / tmp;
  um(1) = ((a(0) * a(1) - a(2)) * s(0) +    (1.0 + a(1)*a(1)) * s(1) + (a(1) * a(2) + a(0)) * s(2)) / tmp;
  um(2) = ((a(0) * a(2) + a(1)) * s(0) + (a(1) * a(2) - a(0)) * s(1) +    (1.0 + a(2)*a(2)) * s(2)) / tmp;
        
  RefPartP_suv_m = um + 0.5 * itsReference.getQ()  * itsBunch->getdT() / itsReference.getM() * c * externalE;

}


inline bool ParallelTTracker::getExternalField(const Vector_t &R, const double &t, Vector_t &extE, Vector_t &extB)
{
  bool bends = false;

  extE = Vector_t(0.0);
  extB = Vector_t(0.0);

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
                  bends = bends || (*element_it)->bends(); // Is the reference particle in a bending field?
                  (*element_it)->apply(R, t, extE, extB);
                }
              if (bends)
                {
                  extE = TransformBack(extE, (*sections_it).Orientation);
                  extB = TransformBack(extB, (*sections_it).Orientation);
                }
              break;
            }
        }
    }
  return bends;

}
inline void ParallelTTracker::writePhaseSpace(const double &sposRef)
{
  Vector_t externalE, externalB;
  Vector_t FDext[6];  // FDext = {BHead, EHead, BRef, ERef, BTail, ETail}.

  // Sample fields at (0,0,zmin), (0,0,zmax) and the centroid location. We
  // are sampling the electric and magnetic fields at the back, front and
  // center of the beam.
  Vector_t rmin, rmax;
  itsBunch->get_bounds(rmin,rmax);

  Vector_t pos[3];
  pos[0] = Vector_t(rmax(0), rmax(1), rmax(2));
  pos[1] = Vector_t(0.0, 0.0, sposRef);
  pos[2] = Vector_t(rmin(0), rmin(1), rmin(2));

  for(int k = 0; k < 3; ++k) 
    {
      getExternalField(pos[k], itsBunch->getT() - 0.5 * itsBunch->getdT(), externalE, externalB);

      FDext[2*k]   = externalB;
      FDext[2*k+1] = externalE;
    }

  // Write fields to .h5 file.
  itsDataSink->writePhaseSpace(*itsBunch, FDext, rmax(2), sposRef, rmin(2));
  *gmsg << *itsBunch << endl;
  //                   itsBunch->printBinHist();
}

#endif // OPAL_ParallelTTracker_HH
