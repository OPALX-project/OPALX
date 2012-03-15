#ifndef OPAL_ParallelCyclotronTracker_HH
#define OPAL_ParallelCyclotronTracker_HH

// ------------------------------------------------------------------------
// $RCSfile: ParallelCyclotronTracker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelCyclotron
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

class BMultipoleField;
class PartBunch;
class PlanarArcGeometry;


// Class ParallelCyclotronTracker
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

enum CyclOperationModeT {SINGLEP, MULTIP, TUNECALC};

class ParallelCyclotronTracker: public Tracker {

public:

  typedef pair<double[8],Component*>        element_pair;
  typedef pair<string, element_pair>       type_pair;
  typedef list<type_pair*>                 beamline_list;
  /// Constructor.
  //  The beam line to be tracked is "bl".
  //  The particle reference data are taken from "data".
  //  The particle bunch tracked is initially empty.
  //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
  //  If [b]revTrack[/b] is true, we track against the beam.
  explicit ParallelCyclotronTracker(const Beamline &bl, const PartData &data,
		       bool revBeam, bool revTrack);

  /// Constructor.
  //  The beam line to be tracked is "bl".
  //  The particle reference data are taken from "data".
  //  The particle bunch tracked is taken from [b]bunch[/b].
  //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
  //  If [b]revTrack[/b] is true, we track against the beam.
  explicit ParallelCyclotronTracker(const Beamline &bl, PartBunch &bunch, DataSink &ds,
				    const PartData &data, bool revBeam, bool revTrack, int maxSTEPS, int timeIntegrator);

  virtual ~ParallelCyclotronTracker();

  /// Apply the algorithm to a Cyclotorn
  virtual void visitCyclotron(const Cyclotron &cycl);

  /// Apply the algorithm to a RFCavity.
  virtual void visitRFCavity(const RFCavity &);

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

  /// set total number of tracked bunches
  inline void setNumBunch(int n){ numBunch_m = n; }

  /// get total number of tracked bunches
  inline int  getNumBunch(){ return numBunch_m; }
  
  /// set the working sub-mode for multi-bunch mode: "FORCE" or "AUTO"  
  inline void  setMultiBunchMode(const int flag ){multiBunchMode_m = flag; }

  /// set last dumped step 
  inline int setLastDumpedStep(const int para ){lastDumpedStep_m = para ; }
  
  /// set the control parameter for "AUTO" sub-mode
  inline void  setParaAutoMode(const double para ){CoeffDBunches_m = para; }

  /// set coefficients for track (default value is for proton)
  void setTrackCoeff(const double para);

 private:

  // Not implemented.
  ParallelCyclotronTracker();
  ParallelCyclotronTracker(const ParallelCyclotronTracker &);
  void operator=(const ParallelCyclotronTracker &);

  beamline_list FieldDimensions;
  list<Component*> myElements;
  int LastVisited;
  Beamline *itsBeamline;

  PartBunch *itsBunch;

  DataSink *itsDataSink;

  /// The maximal number of steps the system is integrated
  int maxSteps_m;

  /// The scale factor for dimensionless variables
  double scaleFactor_m;

  double referenceR;
  double referenceTheta;
  double referenceZ;

  double referencePr;
  double referencePt;
  double referencePz;
  double referencePtot;

  double sinRefTheta;
  double cosRefTheta;
  /// The number of bunches specified in TURNS of RUN commond 
  int numBunch_m;

  // 0 for single bunch (default),
  // 1 for FORCE,
  // 2 for AUTO 
  int multiBunchMode_m;
  
  // control parameter for AUTO multi-bunch mode 
  double CoeffDBunches_m;

  int lastDumpedStep_m;
  // Charge-Mass ratio in [proton unit] 
  double ratioCh_M_m; 

  double chtmc_m;
  double chtm_m;
  
  double PathLength_m;

  // the name of time integrator
  // The ID of time integrator 
  // 0 --- RK-4(default) 
  // 1 --- LF-2
	int  timeIntegrator_m;

  void Tracker_LF();
  void Tracker_RK4();

  /*
   Local Variables both used by the integration methods
  */
	
  Vector_t rold_m, pold_m, rnew_m, ptmp_m;

  long long step_m;
  long long restartStep0_m;

  // temporal 6 phase space varibles of particle [x,y,z,px,py,pz]. Unit: mm & dimensionless
  double variable_m[6];
  // temporal 3 real space varibles of particle ID=0 [x,y,z]. for tune with SC.  Unit: mm
  Vector_t variable_tune0_m;
  // temporal 3 real space varibles of particle ID=1 [x,y,z]. for tune with SC.  Unit: mm
  Vector_t variable_tune1_m;

  // vector of [angle, x, y] of SEO read in from external file for tune with SC. Unit : rad, mm
  vector<Vector_t> variable_SEO_m;

  // save initial phase space distribution (in global Cartesian frame ) for multi-bunch simultion. FixMe: not used
  Vector_t *initialR_m, *initialP_m;

  // record how many bunches has already been injected. ONLY FOR MPM
  int BunchCount_m;

  // decide how many energy bins. ONLY FOR MPM
  // For the time being, we set bin number equal to bunch number.
  int BinCount_m;

  // used for automatic injection in multi-bunch mode
  double RLastTurn_m , RThisTurn_m;

  // start time to record tune data
  double StartTime_m;
 
  // external field arrays for dumping
  Vector_t FDext_m[2], extE_m, extB_m;

  // mark the dumpstep to inject new bunch from here for AUTO mode of restart run of multibunch
  int backupDumpStep_m;

  // flag to determine whether the tune of betatron oscillation is calculated or not for many paticles
  // FixMe: read in from input file
  bool flagDoTune_m;

  const int myNode_m; 
  const size_t initialLocalNum_m; 
  const size_t initialTotalNum_m;

  ofstream outfTheta0_m;
  ofstream outfTheta1_m;
  ofstream outfTheta2_m;
  ofstream outfThetaEachTurn_m;

  void openFiles(string fn);
  void closeFiles();



  // Fringe fields for entrance and exit of magnetic elements.
  void applyEntranceFringe(double edge, double curve,
			   const BMultipoleField &field, double scale);
  void applyExitFringe(double edge, double curve,
		       const BMultipoleField &field, double scale);

  void buildupFieldList(double BcParameter[], string ElementType, Component* elptr);

  bool derivate(double *y, double t, double *yp, int Pindex);

  bool rk4(double x[], double t, double tau,int Pindex);

  // angle range [0~2PI) degree
  double calculateAngle(double x, double y);
  // angle range [-PI~PI) degree
  double calculateAngle2(double x, double y);

  void readOneBunch(const int BeamCount, const int step);
  
  bool checkGapCross(Vector_t Rold, Vector_t Rnew, Component *elptr, double &DistOld);

  bool getTunes(vector<double> &t,  vector<double> &r,  vector<double> &z,int lastTurn,double &nur,double &nuz);

  IpplTimings::TimerRef IntegrationTimer_m;
  IpplTimings::TimerRef DumpTimer_m ;
  IpplTimings::TimerRef TransformTimer_m;
  IpplTimings::TimerRef BinRepartTimer_m;
  
};

#endif // OPAL_ParallelCyclotronTracker_HH
