// ------------------------------------------------------------------------
// $RCSfile: ParallelTTracker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelTTracker
//   The visitor class for tracking particles with time as independent 
//   variable.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include <cfloat>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Algorithms/ParallelTTracker.h"

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

#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "BeamlineGeometry/RBendGeometry.h"
#include "Beamlines/Beamline.h"

#include "Fields/BMultipoleField.h"
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FTpsMath.h"
#include "FixedAlgebra/FVps.h"

#include "Physics/Physics.h"

#include "Utilities/NumToStr.h"
#include "Distribution/Distribution.h"

class Beamline;
class PartData;
using Physics::c;

#define PSdim 6
typedef FVector<double,PSdim> Vector;
typedef FMatrix<double,PSdim,PSdim> Matrix;
typedef FTps<double,PSdim> Series;
typedef FVps<double,PSdim> Map, VSeries;
typedef FMatrix<FTps<double,PSdim>,PSdim,PSdim> MxSeries;

static Vector implicitIntStep(const Vector &zin, const VSeries &f, double s, double ds,
                              int nx = 20);
static Vector implicitInt2(const Vector &zin, const VSeries &f, double s, double ds,
                           int nx = 20, int cx = 4);
static Vector implicitInt4(const Vector &zin, const VSeries &f, double s, double ds,
                           int nx = 20, int cx = 4);
static Vector fixedPointInt2(const Vector &zin, const VSeries &f, double ds,
                             int nx = 50);
static Vector fixedPointInt4(const Vector &zin, const VSeries &f, double s, double ds,
                             int nx = 50, int cx = 4);

extern Inform* gmsg2all;

// Class ParallelTTracker
// ------------------------------------------------------------------------

ParallelTTracker::ParallelTTracker(const Beamline &beamline,
                                   const PartData &reference,
		                           bool revBeam, 
                                   bool revTrack):
  Tracker(beamline, reference, revBeam, revTrack),
  myElements(),
  myFieldList()
{
  itsBeamline = dynamic_cast<Beamline*>(beamline.clone());
}


ParallelTTracker::ParallelTTracker(const Beamline &beamline,
                                   PartBunch &bunch,
                                   DataSink &ds,
                                   const PartData &reference,
                                   bool revBeam,
                                   bool revTrack,
                                   int maxSTEPS):
  Tracker(beamline, reference, revBeam, revTrack),
  myElements(),
  myFieldList(),
  maxSteps_m(maxSTEPS)
{ 
  itsBeamline = dynamic_cast<Beamline*>(beamline.clone());
  itsBunch = &bunch;
  itsDataSink = &ds;
  scaleFactor_m = itsBunch->getdT() * c;

  timeIntegrationTimer1_m  = IpplTimings::getTimer("Time integration1");
  timeIntegrationTimer2_m  = IpplTimings::getTimer("Time integration2");
  timeFieldEvaluation_m  = IpplTimings::getTimer("Field evaluation");
  
  BinRepartTimer_m   = IpplTimings::getTimer("Time of Binary repart.");
  WakeFieldTimer_m   = IpplTimings::getTimer("Time of Wake Field calc.");

}


ParallelTTracker::~ParallelTTracker()
{
  for (FieldListType1Iterator compindex = myElements.begin(); compindex != myElements.end(); compindex++)
    {
      delete (*compindex).Element;
    }
  myElements.clear();
  myFieldList.clear();
  delete itsBeamline;
}

void ParallelTTracker::visitAlignWrapper(const AlignWrapper &wrap)
{
  //   *gmsg << "In visitAlignWrapper " << wrap.getName() << " of type " << wrap.getType() <<  endl;
  if (wrap.getType() == "beamline")
    {
      Beamline *bl = dynamic_cast<Beamline*>(wrap.getElement());
      bl->iterate(*this, false);
    }
  else
    {
      wrap.getElement()->accept(*this);
    }
}

void ParallelTTracker::visitBeamBeam(const BeamBeam &)
{
  //  *gmsg << "In BeamBeam; "<< endl;
}


void ParallelTTracker::visitCollimator(const Collimator &coll)
{
  //  *gmsg << "In Collimator; L= " << coll.getElementLength() << endl;
  myElements.push_back(ParallelTTracker::FieldListType1Entry(dynamic_cast<Collimator*>(coll.clone()),0.,0.));
}


void ParallelTTracker::visitCorrector(const Corrector &corr)
{
  //  *gmsg << "In Corrector; L= " << corr.getElementLength() << endl;
  myElements.push_back(FieldListType1Entry(dynamic_cast<Corrector*>(corr.clone()),0.,0.));
}


void ParallelTTracker::visitDiagnostic(const Diagnostic &diag)
{
  //   *gmsg << "In Diagnostic; L= " << diag.getElementLength() << endl;
  myElements.push_back(FieldListType1Entry(dynamic_cast<Diagnostic*>(diag.clone()),0.,0.));
}


void ParallelTTracker::visitDrift(const Drift &drift)
{
  //   *gmsg << "In drift L= " << drift.getElementLength() << endl;
  myElements.push_back(FieldListType1Entry(dynamic_cast<Drift*>(drift.clone()),0.,0.));
  
}


void ParallelTTracker::visitLambertson(const Lambertson &lamb)
{
  //   *gmsg << "In Lambertson; L= " << lamb.getElementLength() << endl;
  myElements.push_back(FieldListType1Entry(dynamic_cast<Lambertson*>(lamb.clone()),0.,0.));
}


void ParallelTTracker::visitMarker(const Marker &marker)
{
  //   *gmsg << "In Marker; L= " << marker.getElementLength() << endl;
  myElements.push_back(FieldListType1Entry(dynamic_cast<Marker*>(marker.clone()),0.,0.));
}


void ParallelTTracker::visitMonitor(const Monitor &mon)
{
  //   *gmsg << "In Monitor; L= " << corr.getElementLength() << endl;
  Monitor *elptr = dynamic_cast<Monitor*>(mon.clone()->removeWrappers());
  if (!elptr->hasAttribute("ELEMEDGE"))
    {
      *gmsg << "Multipole: no position of the element given!" << endl;
      return;
    }
  double startField = elptr->getAttribute("ELEMEDGE");
  double endField;
  elptr->initialise(itsBunch, startField, endField, 1.0);
  myElements.push_back(FieldListType1Entry(elptr,startField, endField));
}


void ParallelTTracker::visitMultipole(const Multipole &mult)
{
  //   *gmsg << "In Multipole; L= " << mult.getElementLength() << endl;
  Multipole *elptr = dynamic_cast<Multipole*>(mult.clone()->removeWrappers());
  if (!elptr->hasAttribute("ELEMEDGE"))
    {
      *gmsg << "Multipole: no position of the element given!" << endl;
      return;
    }
  double startField = elptr->getAttribute("ELEMEDGE");
  double endField;
  elptr->initialise(itsBunch, startField, endField, 1.0);
  myElements.push_back(FieldListType1Entry(elptr,startField, endField));
}


void ParallelTTracker::visitRBend(const RBend &bend)
{
  //   *gmsg << "In RBend; L= " << bend.getElementLength() << endl;
  RBend *elptr = dynamic_cast<RBend*>(bend.clone()->removeWrappers());
  if (!elptr->hasAttribute("ELEMEDGE"))
    {
      *gmsg << "RBend: no position of the element given!" << endl;
      return;
    }
  double startField = elptr->getAttribute("ELEMEDGE");
  double endField;
  elptr->initialise(itsBunch, startField, endField, 1.0);
  myElements.push_back(FieldListType1Entry(elptr,startField, endField));

}


void ParallelTTracker::visitRFCavity(const RFCavity &as)
{
  //   *gmsg << "In RFCavity; L= " << as.getElementLength() << endl;
  RFCavity *elptr = dynamic_cast<RFCavity*>(as.clone());
  if (!elptr->hasAttribute("ELEMEDGE")){
    *gmsg << "RFCavity: no position of the element or no length of the field given!" << endl;
    return;
  }
  double startField = elptr->getAttribute("ELEMEDGE");
  double endField;

  elptr->initialise(itsBunch, startField, endField, 1.0);
  myElements.push_back(FieldListType1Entry(elptr,startField, endField));
}

void ParallelTTracker::visitTravelingWave(const TravelingWave &as)
{
  //   *gmsg << "In RFCavity; L= " << as.getElementLength() << endl;
  TravelingWave *elptr = dynamic_cast<TravelingWave*>(as.clone());
  if (!elptr->hasAttribute("ELEMEDGE")){
    *gmsg << "RFCavity: no position of the element or no length of the field given!" << endl;
    return;
  }
  double startField = elptr->getAttribute("ELEMEDGE");
  double endField;

  elptr->initialise(itsBunch, startField, endField, 1.0);
  myElements.push_back(FieldListType1Entry(elptr,startField, endField));
}


void ParallelTTracker::visitRFQuadrupole(const RFQuadrupole &rfq)
{
  //   *gmsg << "In RFQuadrupole; L= " << rfq.getElementLength() << endl;
  myElements.push_back(FieldListType1Entry(dynamic_cast<RFQuadrupole*>(rfq.clone()),0.,0.));
}

void ParallelTTracker::visitSBend(const SBend &bend)
{
  //   *gmsg << "In SBend; L= " << bend.getElementLength() << endl;
  myElements.push_back(FieldListType1Entry(dynamic_cast<SBend*>(bend.clone()),0.,0.));
}


void ParallelTTracker::visitSeparator(const Separator &sep)
{
  myElements.push_back(FieldListType1Entry(dynamic_cast<Separator*>(sep.clone()),0.,0.));

}


void ParallelTTracker::visitSeptum(const Septum &sept)
{
  myElements.push_back(FieldListType1Entry(dynamic_cast<Septum*>(sept.clone()),0.,0.));
}


void ParallelTTracker::visitSolenoid(const Solenoid &solenoid)
{
  Solenoid *elptr = dynamic_cast<Solenoid*>(solenoid.clone());
  if (!elptr->hasAttribute("ELEMEDGE"))
    {
      *gmsg << "Solenoid: no position of the element given!" << endl;    
      return;
    }
  
  double startField = elptr->getAttribute("ELEMEDGE");
  double endField;

  elptr->initialise(itsBunch, startField, endField, 1.0);

  myElements.push_back(FieldListType1Entry(elptr, startField, endField));
}

void ParallelTTracker::applyEntranceFringe(double angle, double curve,
                                           const BMultipoleField &field, double scale)
{
}


void ParallelTTracker::applyExitFringe(double angle, double curve,
                                       const BMultipoleField &field, double scale)
{
}

void ParallelTTracker::buildupFieldList()
{
  /** Build up a list of sections of fields: 
   *  for each section we get a list of pointers to the elements which contribute to
   *  the electromagnetic field in this section. Sections start and end where the fields
   *  of the elements start and end respectively. The field of an element can be conained 
   *  in one (at least) or more sections.
   */

  list<double> StartEnd;
  list<Component*> tmp;
  FieldListType1Iterator fld_it;
  list<double>::iterator pos_it, next_it, last_it;
  double tolerance = 1.e-8;
  
  myElements.sort(FieldListType1Entry::SortAscByStart);

  for (fld_it = myElements.begin(); fld_it != myElements.end(); fld_it++)
    {
      StartEnd.push_back((*fld_it).Start);
      StartEnd.push_back((*fld_it).End);
    }
  StartEnd.sort();
  next_it = StartEnd.begin(); next_it++;
  for (pos_it = StartEnd.begin(); next_it != StartEnd.end(); pos_it++, next_it++)
    if (*next_it - *pos_it < tolerance) *next_it = *pos_it;
  StartEnd.unique();
  
  next_it = StartEnd.begin(); next_it++;
  
  for (pos_it = StartEnd.begin(); next_it != StartEnd.end(); pos_it++, next_it++)
    {
      tmp.clear();
      for (fld_it = myElements.begin(); fld_it != myElements.end(); fld_it++)
        {
          if ((*fld_it).Start <= *pos_it + tolerance && (*fld_it).End >= *next_it - tolerance)
            tmp.push_back((*fld_it).Element);
          else if ((*fld_it).Start >= *next_it)
            break;
        }
      if (tmp.size() > 0)
        myFieldList.push_back(FieldListType2Entry(tmp,*pos_it,*next_it));
    }
  *gmsg << "--- BEGIN FIELD LIST ----------------------------------------------------------------------\n" << endl;
  for (FieldListType2Iterator fld2_it = myFieldList.begin(); fld2_it != myFieldList.end(); fld2_it++)
    {
      *gmsg << "--- " << (*fld2_it).Start << " m -- " << (*fld2_it).End << " m ---------------------------\n";
      for (list<Component*>::iterator el_it = (*fld2_it).Elements.begin(); el_it != (*fld2_it).Elements.end(); el_it++)
        *gmsg << (*el_it)->getName() << '\n';
    }
  *gmsg << "--- END   FIELD LIST ----------------------------------------------------------------------\n" << endl;

  /* there might be elements with length zero or extremely short ones. 
     we set them 'goneLive' and 'goneOff' such that they don't appear in 
     the simulation.
  */
  for (fld_it = myElements.begin(); fld_it != myElements.end(); fld_it++)
    if (fabs((*fld_it).End - (*fld_it).Start) < tolerance)
      {
        (*fld_it).goneLive = true;
        (*fld_it).goneOff = true;
      }

}
// 2007/04/19 CKR
void ParallelTTracker::execute(){
  // prepare the tracker for doing its job
  // prepare the elements of the line to be tracked, i.e. load fieldmaps and so on
  //  itsBeamline->arrange();

  double t = itsBunch->getT(); 
  double dt = itsBunch->getdT();
  const Vector_t vscaleFactor = Vector_t(scaleFactor_m);

  //  Fieldmap::getFieldmap(string("Bend.txt"),false);
  //  Fieldmap::readMap(string("Bend.txt"));

  double tEmission = itsBunch->getTEmission();
  int gunSubTimeSteps = 10;

  Vector_t um, a, s, externalE, externalB;
  bool EndOfLineReached;

  double tmp;
  double recpgamma, gamma;
  int emissionSteps = 0;
  long long step =  OPAL.inRestartRun() ? OPAL.getRestartStep() + 1 :lround(t/dt);
  
  bool partOutOfBounds;
  bool bends;               // flag which indicates wheter any particle is within the influence of bending element.
                            // if this is the case we track the reference particle as if it were a real particle, 
                            // otherwise the reference particle is defined as the centroid particle of the bunch

  Vector_t extE, extB;

  size_t totalParticles_i = itsBunch->getTotalNum();
  Inform msg ("bin- ",INFORM_ALL_NODES);

  ElementIterator element_it;
  FieldListType2Iterator fieldlist_it;
  Vector_t CentroidLocation;
  Vector_t CentroidMomentum;
  Vector_t rmin, rmax, rtmp;
  Vector_t EulerAngles;
  double AbsMomentum;
  double AbsMomentumProj; // projection of the momentum onto the z axis
  double RefPartT00, RefPartT01, RefPartT02, 
    RefPartT10, RefPartT11, RefPartT12, 
    RefPartT20, RefPartT21, RefPartT22;                 // defines a rotation matrix

  bool doWakeField = (Options::wakeCalcStep != -1);

  if (itsBunch->doEmission()) {
    emissionSteps = (int)itsBunch->pbin_m->getNBins()*gunSubTimeSteps;
    *gmsg << "Do emission for " << tEmission << " [s] using " << itsBunch->pbin_m->getNBins() << " energy bins " << endl;
    *gmsg << "Change dT from " <<  itsBunch->getdT() << " [s] to "; 
    
    //x    itsBunch->setdT(tEmission/emissionSteps);
    
    *gmsg <<  itsBunch->getdT() << " [s] during emission " << endl;; 
      
    if (Ippl::myNode() == 0) {
      //FIXME: the next step should be done properly when intialising the distribution; 
      //       right now we push the bunch back such that opal behaves like impactt
      //       which is a quite strange behaving...

      // if we dont want opal to behave like impactt move the bunch to the cathode; 
      // otherwise comment the next line out; 
      RealVariable *ar = dynamic_cast<RealVariable *>(OPAL.find("TSHIFT"));
      if (ar) {
        t = itsBunch->calcTimeDelay(ar->getReal());
        *gmsg << "The time delay is " << t << " [s]" << endl;
        itsBunch->setT(t);
      } 
    }
  }
  
  *gmsg << "executing ParallelTTracker, initial DT " << itsBunch->getdT() << " [s] max integration steps " << maxSteps_m << " step= " << step << endl;
  *gmsg << "the mass is: " << itsReference.getM() << ", its charge: " << itsReference.getQ() << endl;
  itsBeamline->accept(*this);

  // set dt for all particles allready in the simulation,
  // i.e. when doing a restarted simulation
  for(int i=0; i < itsBunch->getLocalNum(); i++)
    itsBunch->dt[i] = itsBunch->getdT();

  buildupFieldList();
  if (OPAL.inRestartRun())
    {
      for (int i = 0; i < itsBunch->getLocalNum(); i++) 
        {
          itsBunch->LastSection[i] = 0;
          for (int l = 0; l < myFieldList.size(); l++)
            {
              if (itsBunch->R[i](2) >= myFieldList[l].Start && itsBunch->R[i](2) <= myFieldList[l].End)
                {
                  itsBunch->LastSection[i] = l;
                  break;
                }
            }
        }      
    }
  
  // use scale factor to get dimensionless variables
  // OBSOLETE?: scaling is now down with particel dt
  // itsBunch->R /= vscaleFactor;
  // itsBunch->calcBeamParameters();
  
  // TODO: rescale P = beta_{x,y,z} gamma IF using a distribution
  // (NOT here)

  if (!(itsBunch->weHaveBins())) {
    IpplTimings::startTimer(BinRepartTimer_m);
    itsBunch->do_binaryRepart();
    IpplTimings::stopTimer(BinRepartTimer_m);
    Ippl::Comm->barrier();
  }

  itsBunch->calcBeamParameters();

  if (!OPAL.hasBunchAllocated()) {
    CentroidLocation = itsBunch->get_rmean();
    CentroidMomentum = itsBunch->get_pmean();

    /* build up rotation matrix such that after applying it to the positions and momenta of the particles
       the centroid momentum points in z direction. */
  
    AbsMomentum = sqrt(dot(CentroidMomentum,CentroidMomentum));
    AbsMomentumProj = sqrt(CentroidMomentum(0) * CentroidMomentum(0) + CentroidMomentum(2) * CentroidMomentum(2));
  
    RefPartT00 = CentroidMomentum(2) / AbsMomentumProj;
    RefPartT01 = -CentroidMomentum(0) * CentroidMomentum(1)/(AbsMomentum * AbsMomentumProj);
    RefPartT02 = CentroidMomentum(0) / AbsMomentum;
    RefPartT10 = 0.;
    RefPartT11 = AbsMomentumProj / AbsMomentum;
    RefPartT12 = CentroidMomentum(1) / AbsMomentum;
    RefPartT20 = -CentroidMomentum(0) / AbsMomentumProj;
    RefPartT21 = -CentroidMomentum(2) * CentroidMomentum(1) / (AbsMomentum * AbsMomentumProj);
    RefPartT22 = CentroidMomentum(2) / AbsMomentum;
    
    EulerAngles = Vector_t(-atan(CentroidMomentum(0)/CentroidMomentum(2)), \
                           -asin(CentroidMomentum(1)/AbsMomentum), \
                           0.0);
    // rotate the bunch
    itsBunch->rotateAbout(CentroidLocation,EulerAngles);   
    // move the bunch such that the new centroid location is at (0,0,z)
    itsBunch->moveBy(Vector_t(-CentroidLocation(0),-CentroidLocation(1),0.0));
    
    itsBunch->calcBeamParameters();
    itsBunch->RefPart_R = itsBunch->get_rmean();
    itsBunch->RefPart_P = itsBunch->get_pmean();
  }
  
  CentroidMomentum = itsBunch->get_pmean();
  CentroidLocation = itsBunch->get_rmean();


  /* Activate all elements which influence the particles when the simulation starts;
   *  mark all elements which are allready past.
   */
  itsBunch->get_bounds(rmin,rmax);
  double margin = 3. * itsBunch->RefPart_P(2) * scaleFactor_m / sqrt(1.0 + dot(CentroidMomentum, CentroidMomentum));
  margin = 0.01 > margin? 0.01: margin;

  for (FieldListType1Iterator fld_it = myElements.begin(); fld_it != myElements.end(); fld_it++)
    {
      if (!(*fld_it).goneLive && rmax(2) + margin > (*fld_it).Start  && (rmin(2) - margin) < (*fld_it).End)
        {
          (*fld_it).Element->goOnline();
          *gmsg << (*fld_it).Element->getName() << " of type " << (*fld_it).Element->getType() << " gone live at step #" << step << " rmax = " << rmax(2) << endl;
          (*fld_it).goneLive = true;
        }
      if ((rmin(2) - margin) > (*fld_it).End)
        {
          (*fld_it).goneLive = true;
          (*fld_it).goneOff = true;
        }
    }

  /** 
      Get values for some varibales
      for debug purposes
  */

  double minBinEmitted  = 10.0;
  RealVariable *ar = dynamic_cast<RealVariable *>(OPAL.find("MINBINEMITTED"));
  if (ar) {
    minBinEmitted = ar->getReal();
    *gmsg << "MINBINEMITTED " << minBinEmitted << endl;
  } 

  double minStepforReBin  = 200.0;
  RealVariable *br = dynamic_cast<RealVariable *>(OPAL.find("MINSTEPFORREBIN"));
  if (br) {
    minStepforReBin = br->getReal();
    *gmsg << "MINSTEPFORREBIN " << minStepforReBin << endl;
  } 

  for(step; step < maxSteps_m; step++){
    EndOfLineReached = true;  //check if any particle hasn't reached the end of the field from the last element
    bends = false;
    IpplTimings::startTimer(timeIntegrationTimer1_m);


    //reset E and B to Vector_t(0.0) for every step
    itsBunch->Ef = Vector_t(0.0);
    itsBunch->Bf = Vector_t(0.0);

    for (int i = 0; i < itsBunch->getLocalNum(); i++) {
      //scale each particle with own dt
      itsBunch->R[i] /= vscaleFactor;
      
      recpgamma = 1.0 / sqrt(1.0 + dot(itsBunch->P[i],itsBunch->P[i]));
      
      /** \f[ \vec{x}_{n+1/2} = \vec{x}_{n} + \frac{1}{2}\vec{v}_{n-1/2}\quad (= \vec{x}_{n} + \frac{\Delta t}{2} \frac{\vec{\beta}_{n-1/2}\gamma_{n-1/2}}{\gamma_{n-1/2}}) \f]
       *
       * \code
       * R[i] += 0.5 * P[i] * recpgamma; 
       * \endcode
       */
      itsBunch->R[i] += 0.5 * recpgamma * itsBunch->P[i];
    }
    
    if(totalParticles_i > minBinEmitted)
      itsBunch->boundp();
    
    IpplTimings::stopTimer(timeIntegrationTimer1_m);

    itsBunch->calcBeamParameters();
    
    /** \f[ Space Charge  \f]

    */

    if(itsBunch->hasFieldSolver() && totalParticles_i > minBinEmitted) {

      if (totalParticles_i > 1000) {
        IpplTimings::startTimer(BinRepartTimer_m);
        itsBunch->do_binaryRepart();
        IpplTimings::stopTimer(BinRepartTimer_m);
        Ippl::Comm->barrier();
      }
      if (itsBunch->weHaveBins()) {
        itsBunch->calcGammas();
        for (int b=0; b<=itsBunch->getLastemittedBin() && b<itsBunch->getNumBins()-1; b++) {
          itsBunch->setBinCharge(b,itsBunch->getChargePerParticle());
          itsBunch->computeSelfFields(b);
        } 
        itsBunch->Q = itsBunch->getChargePerParticle();
      }
      else {
        itsBunch->computeSelfFields();
      }
    }

    IpplTimings::startTimer(timeIntegrationTimer2_m);

    /* 
       transport and emit particles 
       that would pass the cathode in the
       first half-step or that already have
       passed it

       to make IPPL and the field solver happy
       make sure that at least 10 particles are emitted

       also remember that node 0 has 
       all the particles to be emitted

       this has to be done *after* the calculation of the
       space charges!
    */
    
    if ((itsBunch->weHaveBins())) { 
      int ne = 0;
      if (Ippl::myNode() == 0) {
        for (int b=0; b<itsBunch->getNumBins(); b++) 
          ne += itsBunch->emitParticles(b);

      }
      reduce(ne,ne,OpAddAssign());
      totalParticles_i += ne;

      itsBunch->updateBinStructure(); // we need this because only node 0 is emitting
      if (step > minStepforReBin) {
        itsBunch->calcGammas();
        const double maxdE = itsBunch->getMaxdEBins();
        if (maxdE < itsBunch->getRebinEnergy()) {
          *gmsg << "**********************************************************" << endl;
          *gmsg << "maxdE < " << maxdE << " [keV] we rebin" << endl;
          *gmsg << "**********************************************************" << endl;
          itsBunch->rebin();
        }
      }
    }

    //FIXME: rethink scaling!
    for (int i=0; i < itsBunch->getLocalNum(); i++)
      itsBunch->R[i] *= vscaleFactor;
    
    // push the reference particle by a half step
    recpgamma = 1.0 / sqrt(1.0 + dot(CentroidMomentum, CentroidMomentum));
    itsBunch->RefPart_R += itsBunch->RefPart_P * recpgamma / 2. * scaleFactor_m;

    /*
      This is a experimental implementation of Wake Fields
      Here we assume that we are in an RFCavity and many 
      parameters are hard wired.
    */

    if (doWakeField) {
      IpplTimings::startTimer(WakeFieldTimer_m);
      itsBunch->calcLineDensity();
      /*
        Convolution

      */
      IpplTimings::stopTimer(WakeFieldTimer_m);

    }
    for (int i = 0; i < itsBunch->getLocalNum(); i++) {
      
      partOutOfBounds = false;
      //reset helper vectors in each step since fields are added in getFieldstrenght
      externalB = Vector_t(0.0);
      externalE = Vector_t(0.0);

      //when calculating external Ef and Bf use dimensions
      //itsBunch->R[i] *= vscaleFactor;

      IpplTimings::startTimer(timeFieldEvaluation_m);
      if ( itsBunch->R[i](2) < myFieldList.back().End )
        {
          EndOfLineReached = false;
          if (itsBunch->R[i](2) >= myFieldList[itsBunch->LastSection[i]].Start)
            {
              while (itsBunch->R[i](2) >= myFieldList[itsBunch->LastSection[i]].End)
                itsBunch->LastSection[i] += 1;

              if (itsBunch->R[i](2) >= myFieldList[itsBunch->LastSection[i]].Start)
                {
                  for (ElementIterator el_it = myFieldList[itsBunch->LastSection[i]].Elements.begin(); el_it != myFieldList[itsBunch->LastSection[i]].Elements.end(); el_it++)
                    {

                      // if we take into account wake fields calculate their influences here
                      if (doWakeField && (i == 0) && ((*el_it)->getName() == string("LRF0")))
                        { 
                          IpplTimings::startTimer(WakeFieldTimer_m);		    
                          *gmsg << "applt wake field tp particle " << (*el_it)->getName() << endl;
                          IpplTimings::stopTimer(WakeFieldTimer_m);
                        }
                      
                      // the lines below should not be necessary if the margin which decides when an element goes online is big enough
                      //                       if (!(*el_it)->Online())
                      //                         {
                      //                           *gmsg << "* ************** W A R N I N G *****************************************************" << endl;
                      //                           *gmsg << "* trying to get field from element " << (*el_it)->getName() << " which is not online yet" << endl;
                      //                           *gmsg << "* z = " << itsBunch->R[i](2) << " m; Last section: " << itsBunch->LastSection[i] << endl;
                      //                           *gmsg << "* **********************************************************************************" << endl;
                      //                           (*el_it)->goOnline();
                      //                           *gmsg << (*el_it)->getName() << " gone live at step #" << step << " rmax = " << rmax(2) << endl;
                      //                         }
                      bends = bends || (*el_it)->bends(); // check if any element bends the beam; make sure that this is communicated to all processors!
                      partOutOfBounds = partOutOfBounds || (*el_it)->apply(i, t+itsBunch->dt[i]/2.0, externalE, externalB);
                    }
                }
            }
        }
      IpplTimings::stopTimer(timeFieldEvaluation_m);
      
      // skip rest of the particle push if the
      // particle is out of bounds i.e. does not see
      // a E or B field
      if (partOutOfBounds) {

        externalE = Vector_t(0.0);
        externalB = Vector_t(0.0);
        itsBunch->Ef[i] = Vector_t(0.0);
        itsBunch->Bf[i] = Vector_t(0.0);
        itsBunch->Bin[i] = -1;
      	continue;
      }
      
      itsBunch->Ef[i] += externalE;
      itsBunch->Bf[i] += externalB;

      itsBunch->R[i] /= vscaleFactor;
      

      /** Update the momenta using the \f$D_1\f$ algorithm as described in Birdsall and Langdon's book:
       *
       *  \f[ \vec{v}_{n+1/2} = \frac{1}{2}\overline{\vec{a}_{n}} \; \Delta t + \mathbb{R} \cdot (\vec{v}_{n-1/2} + \frac{1}{2} \overline{\vec{a}_{n}} \; \Delta t) \f]
       *
       *  where the operator \f$\mathbb{R}\f$ effects a rotation through angle \f$-2\tan^{-1}(\vec{\Omega} \; \Delta t/2)\f$ where \f$\vec{\Omega} = Q \vec{B}_{n}/(m_{e} c)\f$.
       *  \f$\mathbb{R}\f$ can be written as 
       *
       *  \f[\mathbb{R} = \frac{(1 - \Theta^2) \; \mathbb{I} + 2 \; \Theta \Theta^{T} - 2 \; \Theta \times \mathbb{I}}{1+\Theta}\f]
       *
       *  where \f$\vec{\Theta} = \vec{\Omega} \; \Delta t/2\f$ and \f$\mathbb{I}\f$ is the unit tensor.
       */

      /** \f[ \vec{u}_{m} = \vec{\beta}_{n} \gamma_{n} = \vec{\beta}_{n-1/2}\; \gamma_{n-1/2} + \frac{q}{m_{e} c} \vec{E} \frac{\Delta t}{2} \f]
       * \code
       * um = P[i] + 0.5 * Q * dt/M * c * E[i]; 
       * \endcode
       */
      um = itsBunch->P[i] + 0.5 * itsReference.getQ() * itsBunch->dt[i] / itsReference.getM() * c * itsBunch->Ef[i];

      /** \f[ \vec{a}_{n} = \frac{Q \; \Delta t}{2\gamma M c^2}\vec{B} \f]
       * \code
       * recpgamma = 1.0 / sqrt(1.0 + um[0]*um[0] + um[1]*um[1] + um[2]*um[2]); 
       * 
       * tmp = 0.5 * Q * recpgamma * dt/M  * c * c;
       * a = tmp * B[i];
       * \endcode
       */
      recpgamma = 1.0 / sqrt(1.0 + dot(um,um));
      
      tmp = 0.5 * itsReference.getQ() * recpgamma * itsBunch->dt[i] / itsReference.getM() * c * c;
      a = tmp * itsBunch->Bf[i];

      /** \f[ \vec{s} = \vec{\beta}_{n} \gamma_{n} + \frac{Q \; \Delta t}{2\gamma_{n} M c^2} \cdot \gamma_{n} \vec{\beta}_{n} \wedge \vec{B}_{n} \f]
       * \code
       * s = um + tmp * cross(um, B[i]);
       * \endcode
       */
      s = um + tmp * cross(um,itsBunch->Bf[i]);

      /** \f[\vec{u}_{m} = \frac{\mathbb{I} \vec{s} + \vec{a}_{n}\vec{a}_{n}^{T} \vec{s} - \vec{a}_{n} \wedge \vec{s}}{1 + a_{n}^{2}} \f]
       *
       * \code
       * tmp = 1.0 + dot(a,a);
       * um(0) = ((1.0 + a(0)*a(0))    * s(0) + (a(0) * a(1) + a(2)) * s(1) + (a(0) * a(2) - a(1)) * s(2)) / tmp;
       * um(1) = ((a(0) * a(1) - a(2)) * s(0) +    (1.0 + a(1)*a(1)) * s(1) + (a(1) * a(2) + a(0)) * s(2)) / tmp;
       * um(2) = ((a(0) * a(2) + a(1)) * s(0) + (a(1) * a(2) - a(0)) * s(1) +    (1.0 + a(2)*a(2)) * s(2)) / tmp;
       * \endcode
       */
      tmp = 1.0 + dot(a,a);

      // since there are no matrices and therefore the multiplication of a vector with the transpose of some other
      // vector does not exist we have to do the next step component wise...
      um(0) = ((1.0 + a(0)*a(0))    * s(0) + (a(0) * a(1) + a(2)) * s(1) + (a(0) * a(2) - a(1)) * s(2)) / tmp;
      um(1) = ((a(0) * a(1) - a(2)) * s(0) +    (1.0 + a(1)*a(1)) * s(1) + (a(1) * a(2) + a(0)) * s(2)) / tmp;
      um(2) = ((a(0) * a(2) + a(1)) * s(0) + (a(1) * a(2) - a(0)) * s(1) +    (1.0 + a(2)*a(2)) * s(2)) / tmp;

      /** \f[ \vec{\beta}_{n+1/2}\; \gamma_{n+1/2}  = \frac{\mathbb{I} \vec{s} + \vec{a}_{n}\vec{a}_{n}^{T} \vec{s} - \vec{a}_{n} \wedge \vec{s}}{1 + a_{n}^{2}} \\
       *                                            = \frac{(1 - (\frac{q}{m_{e} \gamma_{n}} \frac{\Delta t}{2} \vec{B}_{n})^{2}) \; \mathbb{I} + 2 \; \frac{q^{2}}{m_{e}^{2} \; \gamma_{n}^2}\frac{(\Delta t)^{2}}{4} \vec{B} \vec{B}^{T} - 2 \; \frac{q}{m_{e} \gamma_{n}} \frac{\Delta t}{2}\; \vec{B} \wedge \mathbb{I}}{1 + (\frac{q}{m_{e} \gamma_{n}} \frac{\Delta t}{2} \vec{B})^{2}}\; \vec{\beta}_{n} \gamma_{n} \\
       *                                            = \mathbb{R} \vec{\beta}_{n}\; \gamma_{n},\f]  
       * \code
       * P[i] = um + 0.5 * Q * dt / M * c * E[i]; 
       * \endcode
       */
      itsBunch->P[i] = um + 0.5 * itsReference.getQ()  * itsBunch->dt[i] / itsReference.getM() * c * itsBunch->Ef[i]; 
    } // end particle loop part 1

    reduce(bends,bends,OpOr());

    if (itsBunch->getLocalNum() == 0)
      EndOfLineReached = false;

    /**
       Delete marked particles
    */
    size_t totalParticles_f = totalParticles_i;
    int ne = 0;
    for (int i = 0; i < itsBunch->getLocalNum(); i++) 
      if(itsBunch->Bin[i] < 0)
        {
          itsBunch->destroy(1,i);
          ++ne;
        }
    reduce(ne,ne,OpAddAssign());
    totalParticles_f -= ne;

     
    if (totalParticles_f > 0)
      {
        if (!bends) // none of the particles is in a bending element
          {
            itsBunch->calcBeamParameters();
            CentroidLocation = itsBunch->get_rmean() * scaleFactor_m;
            CentroidMomentum = itsBunch->get_pmean();
             
            /* Update the position of the reference particle in ZXY-coordinates. The angle between the ZXY- and the SUV-coordinate
             *  system is determined by the momentum of the reference particle. We calculate the momentum of the reference
             *  particle by rotating the centroid momentum (= momentum of the reference particle in the SUV-coordinate system).
             *  Then we push the reference particle with this momentum for half a time step.
             */

            recpgamma = 1./sqrt(1.0 + dot(CentroidMomentum, CentroidMomentum));
            
            /* First update the momentum of the reference particle in zxy coordinate system, then update its position     */
            itsBunch->RefPart_P(0) = RefPartT00 * CentroidMomentum(0) + RefPartT01 * CentroidMomentum(1) + RefPartT02 * CentroidMomentum(2);
            itsBunch->RefPart_P(1) = RefPartT10 * CentroidMomentum(0) + RefPartT11 * CentroidMomentum(1) + RefPartT12 * CentroidMomentum(2);
            itsBunch->RefPart_P(2) = RefPartT20 * CentroidMomentum(0) + RefPartT21 * CentroidMomentum(1) + RefPartT22 * CentroidMomentum(2);

            itsBunch->RefPart_R += itsBunch->RefPart_P * recpgamma *  scaleFactor_m / 2.;
            CentroidLocation += CentroidMomentum * recpgamma * scaleFactor_m / 2.;
        
          }
        else /* at least one of the elements bends the beam; until all particles have left the bending elements we track the reference particle
              * as if it were a regular particle; from the moment when the reference particle has reached the bending field until it leaves
              * it again we rotate the bunch about the position of the reference particle such that the momentum of the reference particle points
              * in z direction
              */
          { 
            recpgamma = 1./sqrt(1.0 + dot(CentroidMomentum, CentroidMomentum));
            CentroidLocation += CentroidMomentum * recpgamma * scaleFactor_m / 2.;
            extE = Vector_t(0,0,0);
            extB = Vector_t(0,0,0);
        
            bool Cent_in_bend = false;
            // get field of reference particle
            for (fieldlist_it = myFieldList.begin(); fieldlist_it != myFieldList.end(); fieldlist_it++)
              {
                if (CentroidLocation(2) >= (*fieldlist_it).Start && CentroidLocation(2) < (*fieldlist_it).End)
                  {
                    for (element_it = (*fieldlist_it).Elements.begin(); element_it != (*fieldlist_it).Elements.end(); element_it++)
                      {
                        if (!(*element_it)->Online())
                          {
                            *gmsg << "* ************** W A R N I N G *****************************************************" << endl;
                            *gmsg << "* trying to get field from element " << (*element_it)->getName() << " which is not online yet" << endl;
                            *gmsg << "* **********************************************************************************" << endl;
                          }
                        else
                          {
                            Cent_in_bend = Cent_in_bend || (*element_it)->bends(); // is the reference particle in a bending field?
                            (*element_it)->apply(CentroidLocation, t, extE, extB);
                          }
                      }
                    break;
                  }
              }

            // FIXME: add self field 

            CentroidLocation /= vscaleFactor;
        
            // track reference particle
            um = CentroidMomentum + 0.5 * itsReference.getQ() * itsBunch->getdT() / itsReference.getM() * c * extE;

            recpgamma = 1.0 / sqrt(1.0 + dot(um,um));
        
            tmp = 0.5 * itsReference.getQ() * recpgamma * itsBunch->getdT() / itsReference.getM() * c * c;
            a = tmp * extB;
        
            s = um + tmp * cross(um,extB);
        
            tmp = 1.0 + dot(a,a);
        
            um(0) = ((1.0 + a(0)*a(0))    * s(0) + (a(0) * a(1) + a(2)) * s(1) + (a(0) * a(2) - a(1)) * s(2)) / tmp;
            um(1) = ((a(0) * a(1) - a(2)) * s(0) +    (1.0 + a(1)*a(1)) * s(1) + (a(1) * a(2) + a(0)) * s(2)) / tmp;
            um(2) = ((a(0) * a(2) + a(1)) * s(0) + (a(1) * a(2) - a(0)) * s(1) +    (1.0 + a(2)*a(2)) * s(2)) / tmp;
        
            CentroidMomentum = um + 0.5 * itsReference.getQ()  * itsBunch->getdT() / itsReference.getM() * c * extE;

            if (Cent_in_bend)
              {
                /* Update the position of the reference particle in ZXY-coordinates. The angle between the ZXY- and the SUV-coordinate
                 *  system is determined by the momentum of the reference particle. We calculate the momentum of the reference
                 *  particle by rotating the centroid momentum (= momentum of the reference particle in the SUV-coordinate system).
                 *  Then we push the reference particle with this momentum for half a time step.
                 */
                AbsMomentum = sqrt(itsBunch->RefPart_P(0) * itsBunch->RefPart_P(0) + itsBunch->RefPart_P(1) * itsBunch->RefPart_P(1) + itsBunch->RefPart_P(2) * itsBunch->RefPart_P(2));
                AbsMomentumProj = sqrt(itsBunch->RefPart_P(0) * itsBunch->RefPart_P(0) + itsBunch->RefPart_P(2) * itsBunch->RefPart_P(2));
            
                RefPartT00 = itsBunch->RefPart_P(2) / AbsMomentumProj;
                RefPartT01 = -itsBunch->RefPart_P(0) * itsBunch->RefPart_P(1)/(AbsMomentum * AbsMomentumProj);
                RefPartT02 = itsBunch->RefPart_P(0) / AbsMomentum;
                RefPartT10 = 0.;
                RefPartT11 = AbsMomentumProj / AbsMomentum;
                RefPartT12 = itsBunch->RefPart_P(1) / AbsMomentum;
                RefPartT20 = -itsBunch->RefPart_P(0) / AbsMomentumProj;
                RefPartT21 = -itsBunch->RefPart_P(2) * itsBunch->RefPart_P(1) / (AbsMomentum * AbsMomentumProj);
                RefPartT22 = itsBunch->RefPart_P(2) / AbsMomentum;
            
            
                /* Rotate the bunch about its centroid location by angles determined by its centroid momentum. The angles are such 
                 *  that \f$(p_s,p_u,p_v) = (p,0,0)\f$ after rotation.
                 */
                EulerAngles = Vector_t(-atan(CentroidMomentum(0)/CentroidMomentum(2)), \
                                       -asin(CentroidMomentum(1)/sqrt(dot(CentroidMomentum,CentroidMomentum))), \
                                       0.0);
                itsBunch->rotateAbout(CentroidLocation,EulerAngles);
                itsBunch->calcBeamParameters();

              }

            recpgamma = 1./sqrt(1.0 + dot(CentroidMomentum, CentroidMomentum));
            
            // First update the momentum of the reference particle in zxy coordinate system, then update its position    
            itsBunch->RefPart_P(0) = RefPartT00 * CentroidMomentum(0) + RefPartT01 * CentroidMomentum(1) + RefPartT02 * CentroidMomentum(2);
            itsBunch->RefPart_P(1) = RefPartT10 * CentroidMomentum(0) + RefPartT11 * CentroidMomentum(1) + RefPartT12 * CentroidMomentum(2);
            itsBunch->RefPart_P(2) = RefPartT20 * CentroidMomentum(0) + RefPartT21 * CentroidMomentum(1) + RefPartT22 * CentroidMomentum(2);
            itsBunch->RefPart_R += itsBunch->RefPart_P * recpgamma * scaleFactor_m / 2.;
        
            CentroidMomentum = Vector_t(0.0, 0.0, sqrt(dot(CentroidMomentum,CentroidMomentum)));
            CentroidLocation += CentroidMomentum * recpgamma / 2.;
            CentroidLocation *= vscaleFactor;
          }
      }

    // calculate the dimensions of the bunch and add a small margin to them; then decide which elements have to be triggered
    // when an element is triggered memory is allocated and the field map is read in
    itsBunch->get_bounds(rmin,rmax);

    margin = 3. * itsBunch->RefPart_P(2) * recpgamma;
    // trigger the elements
    margin = 0.01 > margin? 0.01: margin;
     
    for (FieldListType1Iterator fld_it = myElements.begin(); fld_it != myElements.end(); fld_it++)
      {
        //          if (rmax(2) > (*fld_it).Start && rmax(2) < (*fld_it).End)
        //            bends = bends || (*fld_it).Element->bends();
        if ((rmax(2) + margin) * scaleFactor_m > (*fld_it).Start  && !(*fld_it).goneLive)
          {
            (*fld_it).Element->goOnline();  // allocate memory and read in field map
            cout << (*fld_it).Element->getName() << " gone live at step #" << step << " rmax = " << rmax(2) * scaleFactor_m << endl;
            (*fld_it).goneLive = true;
          }
        if ((rmin(2) - margin) * scaleFactor_m > (*fld_it).End && !(*fld_it).goneOff)
          {
            (*fld_it).Element->goOffline();  // free memory if no other element uses the same field map
            cout << (*fld_it).Element->getName() << " gone off at step #" << step  << " rmin = " << rmin(2) * scaleFactor_m << endl;
            (*fld_it).goneOff = true;
          }
      }

    // start particle loop part 2
    for (int i = 0; i < itsBunch->getLocalNum(); i++) 
      {
        if (itsBunch->Bin[i] < 0)
          continue; // skip all particles which were out of bounds
        recpgamma = 1.0 / sqrt(1.0 + dot(itsBunch->P[i],itsBunch->P[i])); 

        /** \f[ \vec{x}_{n+1} = \vec{x}_{n+1/2} + \frac{1}{2}\vec{v}_{n+1/2}\quad (= \vec{x}_{n+1/2} + \frac{\Delta t}{2} \frac{\vec{\beta}_{n+1/2}\gamma_{n+1/2}}{\gamma_{n+1/2}}) \f] 
         * \code
         * R[i] += 0.5 * P[i] * recpgamma; 
         * \endcode
         */
        itsBunch->R[i] += 0.5 * itsBunch->P[i] * recpgamma;
        //and scale back to dimensions
        itsBunch->R[i] *= Vector_t(Physics::c*itsBunch->getdT(),Physics::c*itsBunch->getdT(),Physics::c*itsBunch->dt[i]);
      

      } // end particle loop part 2
    IpplTimings::stopTimer(timeIntegrationTimer2_m);

    if (itsBunch->weHaveBins())
      {
        for (int i = 0; i < itsBunch->getLocalNum(); ++i)
          {
            //reset time step if particle was emitted in the first half-step
            //the particle is now in sync with the simulation timestep
            itsBunch->dt[i] = itsBunch->getdT();
          }
      }

    if(totalParticles_f > minBinEmitted)
      itsBunch->boundp();

    if(totalParticles_i != totalParticles_f) {
      *gmsg << "======================================================================" << endl;
      *gmsg << "Particle loss at t= " << t << " n= " << totalParticles_i - totalParticles_f << " Ntot= " << totalParticles_f << endl;  
      *gmsg << "======================================================================" << endl;
      totalParticles_i = totalParticles_f;
    } 
   
    t += /*0.5*/ itsBunch->getdT(); //t after a full global timestep with dT "synchronization point" for simulation time
    
    itsBunch->setT(t);

    //IFF: cheap step dump regulation
    if ( totalParticles_f > 0)
      { 
        double sposRef = itsBunch->get_sPos();
        if (totalParticles_f <= minBinEmitted)
          {
            *gmsg << "Step " << step << "; only " << totalParticles_f << " particles emitted; t= " << t << " [s]" << endl;
          }
        else if (isnan(sposRef) || isinf(sposRef))
          {
            *gmsg << "Step " << step << "; there seems to be something wrong with the position of the bunch!" << endl;
          }
        else
          {
            *gmsg << "Step " << step << " at " << sposRef << " [m] t= " << t << " [s]" << endl;

            if (step % Options::psDumpFreq == 0 ) 
              {

                //FDext = {BHead, EHead, BRef, ERef, BTail, ETail}
                Vector_t FDext[6];
                double sposHead, sposTail; //, meanEnergy, energy[itsBunch->getTotalNum()];

                for(int k=0; k<6; k++) FDext[k] = Vector_t(0.0,0.0,0.0);
                //sample fields at (0,0,rmin), (0,0,rmax) and the centroid location
                itsBunch->get_bounds(rmin,rmax);

                sposHead = rmax(2);
                sposTail = rmin(2);

                Vector_t pos[3];
                pos[0] = Vector_t(rmax(0),rmax(1),sposHead);
                pos[1] = Vector_t(0.00,0.00,sposRef);
                pos[2] = Vector_t(rmin(0),rmin(1),sposTail);

                for(int k=0; k < 3; k++) 
                  {
                    extE = Vector_t(0,0,0);
                    extB = Vector_t(0,0,0);

                    if ( pos[k](2) <= myFieldList.back().End && pos[k](2) >= myFieldList.front().Start )
                      for (FieldListType2Iterator fl_it = myFieldList.begin(); fl_it != myFieldList.end(); fl_it++)
                        {
                          if (pos[k](2) >= (*fl_it).Start && pos[k](2) < (*fl_it).End)
                            {
                              for (ElementIterator el_it = (*fl_it).Elements.begin(); el_it != (*fl_it).Elements.end(); el_it++)
                                (*el_it)->apply(pos[k], t, extE, extB);
                              break;
                            }
                        }
	
                    FDext[2*k]   = extB;
                    FDext[2*k+1] = extE;
                  }
           
                itsDataSink->writePhaseSpace(*itsBunch, FDext, sposHead, sposRef, sposTail);
                *gmsg << *itsBunch << endl;
              }
            
          }
      }
    else 
      {
        *gmsg << "Step " << step << " no emission yet "  << " t= " << t << " [s]" << endl;
      }

    if(step > emissionSteps) 
      {
        //TEST the non-len reduce: reduce(&EndOfLineReached, &EndOfLineReached, OpBitwiseOrAssign());
        reduce(&EndOfLineReached, &EndOfLineReached + 1, &EndOfLineReached, OpBitwiseOrAssign());
        if (EndOfLineReached) break;
      }
  }
  for (FieldListType1Iterator fld_it = myElements.begin(); fld_it != myElements.end(); fld_it++)
    {
      if ((*fld_it).goneLive && !(*fld_it).goneOff)
        {
          (*fld_it).Element->goOffline();  // free memory if no other element uses the same field map
          cout << (*fld_it).Element->getName() << " gone off at step #" << step  << " rmin = " << rmin(2) * scaleFactor_m << endl;
          (*fld_it).goneOff = true;
        }
    }
  *gmsg << "done executing ParallelTTracker" << endl;
}

// 2007/04/19 CKR
void ParallelTTracker::visitBeamline(const Beamline & bl){
  // or maybe here?
  // for (int step = 0; step < maxstep; ++step){
  itsBeamline->iterate(*dynamic_cast<BeamlineVisitor*>(this),false);//, from, to);
  // }
}


static Vector implicitIntStep(const Vector &zin, const VSeries &f, const MxSeries gradf, double ds, int nx)
{
  //   //std::cerr << "==> In implicitIntStep(zin,f,gradf,ds,nx) ..." << std::endl;
  //   //std::cerr << "  ds = " << ds << std::endl;
  //   //std::cerr << " zin =\n" << zin << std::endl;
  //   // This routine integrates the N-dimensional autonomous differential equation
  //   // z' = f(z) for a single step of size ds, using Newton's method to solve the
  //   // implicit equation zf = zin + ds*f((zin+zf)/2).  For reasons of efficiency,
  //   // its arguments include the matrix gradf = grad(f).  The (optional) argument
  //   // nx limits the number of Newton iterations.  This routine returns a result
  //   // zf accurate through second-order in the step-size ds.  When f derives from
  //   // a Hamiltonian---i.e., f=J.grad(H)---then this routine performs symplectic
  //   // integration.

  //   // Set up flags, etc., for convergence (bounce) test.
  //   FVector<bool,PSdim> bounce(false);
  //   Vector dz,dz_old;
  //   int bcount=0;
  //   static const double thresh=1.e-8;

  //   // Use second-order Runge-Kutta integration to determine a good initial guess.
  //   double ds2=0.5*ds;
  //   Vector z=f.constantTerm(zin);
  //   z=zin+ds2*(z+f.constantTerm(zin+ds*z));

  //   // Newton iterations:
  //   //   z :-> [I-ds/2.grad(f)]^{-1}.[zin+ds.f((zin+z)/2)-ds/2.grad(f).z]
  //   // (A possible method for speeding up this computation would
  //   //  be to recompute grad(f) every n-th step, where n > 1!)
  Vector zf;
  //   int ni=0;
  //   while(bcount<PSdim){
  //     if(ni==nx){
  //       string msg = "Convergence not achieved within " + NumToStr(nx) + " iterations!";
  //       throw ConvergenceError("ParallelTTracker::implicitIntStep()",msg);
  //     }

  //     // Build gf = -ds/2.grad(f)[(zin+z)/2] and idgf_inv = [I-ds/2.grad(f)]^{-1}[(zin+z)/2].
  //     Vector zt=0.5*(zin+z);
  //     Matrix gf,idgf,idgf_inv;
  //     for(int i=0;i<PSdim;++i)
  //       for(int j=0;j<PSdim;++j)
  //         gf[i][j]=-ds2*gradf[i][j].evaluate(zt);
  //     idgf=gf;
  //     for(int i=0;i<PSdim;++i) idgf[i][i]+=1.;
  //     FLUMatrix<double,PSdim> lu(idgf);
  //     idgf_inv = lu.inverse();

  //     // Execute Newton step.
  //     zf = idgf_inv*(zin+ds*f.constantTerm(zt)+gf*z);

  //     //std::cerr << " -(ds/2)grad(f) =\n" << gf << std::endl;
  //     //std::cerr << " f =\n" << f.constantTerm(zt) << std::endl;
  //     //std::cerr << "zk =\n" << zf << std::endl;

  //     // Test for convergence ("bounce" test).
  //     dz_old = dz;
  //     dz = zf - z;
  //     if(ni){// (we need at least two iterations before testing makes sense)
  //       for(int i=0;i<PSdim;++i){
  //         if(!bounce[i] && (dz[i]==0. || (abs(dz[i])<thresh && abs(dz[i])>=abs(dz_old[i]))))
  //           {bounce[i]=true; ++bcount;}
  //       }
  //     }
  //     z=zf;
  //     ++ni;
  //   }

  //   //std::cerr << "  zf =\n" << zf << std::endl;
  //   //std::cerr << "==> Leaving implicitIntStep(zin,f,gradf,ds,nx)" << std::endl;
  return zf;
}


static Vector implicitInt2(const Vector &zin, const VSeries &f, double s, double ds, int nx, int cx)
{
  //   //std::cerr << "==> In implicitInt2(zin,f,s,ds,nx,cx) ..." << std::endl;
  //   //std::cerr << " zin =\n" << zin << std::endl;
  //   // Default: nx = 20, cx = 4

  //   // This routine integrates the N-dimensional autonomous differential equation
  //   // z' = f(z) for a distance s, in steps of size ds, using Newton's method to
  //   // solve the implicit equation zt = zin + ds*f((zin+zt)/2).  The optional
  //   // argument nx limits the number of Newton iterations.  The other optional
  //   // argument limits this routine to cutting the step size in half at most cx
  //   // times, so that the step size may shrink to a size no smaller than ds/(2^cx).
  //   // This routine returns a result zf accurate through second-order in the
  //   // step-size ds.  When f derives from a Hamiltonian---f = J.grad(H)---then this
  //   // routine performs symplectic integration.

  //   // Convergence warning flag.
  //   static bool cnvWarn = false;

  //   // Build matrix grad(f).
  //   MxSeries gradf;
  //   for(int i=0;i<PSdim;++i)
  //     for(int j=0;j<PSdim;++j)
  //       gradf[i][j]=f[i].derivative(j);
  //   //std::cerr << " Series f =\n" << f << std::endl;
  //   //std::cerr << " MxSeries gf =\n" << gradf << std::endl;

  //   // Initialize accumulated length, current step-size, and number of cuts.
  //   double as=abs(s), st=0., dsc=abs(ds);
  //   if(s<0.) dsc = -dsc;
  //   int ci=0;

  //   // Integrate each step.
  Vector zf=zin;
  //   while(abs(st)<as){
  //     Vector zt;
  //     bool ok=true;
  //     try{
  //       if(abs(st+dsc)>as) dsc=s-st;
  //       zt=implicitIntStep(zf,f,gradf,dsc,nx);
  //     }
  //     catch (ConvergenceError &cnverr) {
  //       if(++ci>cx) {
  //         string msg = "Convergence not achieved within " + NumToStr(cx) + " cuts of step-size!";
  //         throw ConvergenceError("ParallelTTracker::implicitInt2()",msg);
  //       }
  //       if (!cnvWarn) {
  //         std::cerr << " <***WARNING***> [ParallelTTracker::implicitInt2()]:\n"
  //                   << "   Cutting step size, a probable violation of the symplectic condition."
  //                   << std::endl;
  //         cnvWarn=true;
  //       }
  //       dsc*=0.5; ok=false;
  //     }
  //     if(ok){zf=zt; st+=dsc;}
  //   }

  //   //std::cerr << "  zf =\n" << zf  << std::endl;
  //   //std::cerr << "==> Leaving implicitInt2(zin,f,s,ds,nx,cx)" << std::endl;
  return zf;
}

static Vector implicitInt4(const Vector &zin, const VSeries &f, double s, double ds, int nx, int cx)
{
  //   //std::cerr << "==> In implicitInt4(zin,f,s,ds,nx,cx) ..." << std::endl;
  //   // Default: nx = 20, cx = 4

  //   // This routine integrates the N-dimensional autonomous differential equation
  //   // z' = f(z) for a distance s, in steps of size ds.  It uses a "Yoshida-fied"
  //   // version of implicitInt2 to obtain zf accurate through fourth-order in the
  //   // step-size ds.  When f derives from a Hamiltonian---i.e., f = J.grad(H)---
  //   // then this routine performs symplectic integration.  The optional arguments
  //   // nx and cx have the same meaning as in implicitInt2().

  //   // Convergence warning flag.
  //   static bool cnvWarn = false;

  //   // The Yoshida constants: 2ya+yb=1; 2ya^3+yb^3=0.
  //   static const double yt=pow(2.,1/3.);
  //   static const double ya=1/(2.-yt);
  //   static const double yb=-yt*ya;

  //   // Build matrix grad(f).
  //   MxSeries gradf;
  //   for(int i=0;i<PSdim;++i)
  //     for(int j=0;j<PSdim;++j)
  //       gradf[i][j]=f[i].derivative(j);

  //   // Initialize accumulated length, current step-size, and number of cuts.
  //   double as=abs(s),st=0., dsc=abs(ds);
  //   if(s<0.) dsc = -dsc;
  //   int ci=0;

  //   // Integrate each step.
  Vector zf=zin;
  //   while(abs(st)<as){
  //     Vector zt;
  //     bool ok=true;
  //     try{
  //       if(abs(st+dsc)>as) dsc=s-st;
  //       zt=implicitIntStep(zf,f,gradf,ya*dsc,nx);
  //       zt=implicitIntStep(zt,f,gradf,yb*dsc,nx);
  //       zt=implicitIntStep(zt,f,gradf,ya*dsc,nx);
  //     }
  //     catch (ConvergenceError &cnverr) {
  //       if(++ci>cx) {
  //         string msg = "Convergence not achieved within " + NumToStr(cx) + " cuts of step-size!";
  //         throw ConvergenceError("ParallelTTracker::implicitInt4()",msg);
  //       }
  //       if (!cnvWarn) {
  //         std::cerr << " <***WARNING***> [ParallelTTracker::implicitInt4()]:\n"
  //                   << "   Cutting step size, a probable violation of the symplectic condition."
  //                   << std::endl;
  //         cnvWarn=true;
  //       }
  //       dsc*=0.5; ok=false;
  //     }
  //     if(ok){zf=zt; st+=dsc;}
  //   }

  //   //std::cerr << "==> Leaving implicitInt4(...)" << std::endl;
  return zf;
}


static Vector fixedPointInt2(const Vector &zin, const VSeries &f, double ds, int nx)
{
  //   //std::cerr << "==> In fixedPointInt2(zin,f,ds,nx) ..." << std::endl;
  //   // Default: nx = 50

  //   //std::cerr << "  ds = " << ds << std::endl;
  //   //std::cerr << " zin =\n" << zin << std::endl;
  //   // This routine integrates the N-dimensional autonomous differential equation
  //   // z' = f(z) for a single step of size ds by iterating the equation
  //   //         z = zin + ds * f((zin+z)/2)
  //   // to find a fixed-point zf for z.  It is accurate through second order in the
  //   // step size ds.

  //   // Set up flags, etc., for convergence (bounce) test.
  //   FVector<bool,PSdim> bounce(false);
  //   Vector dz,dz_old;
  //   int bcount=0;
  //   static const double thresh=1.e-8;

  //   // Iterate z :-> zin + ds * f( (zin + z)/2 ).
  Vector zf;
  //   Vector z=zin;
  //   int ni=0;
  //   while(bcount<PSdim){
  //     if(ni==nx){
  //       string msg = "Convergence not achieved within " + NumToStr(nx) + " iterations!";
  //       throw ConvergenceError("ParallelTTracker::fixedPointInt2()",msg);
  //     }

  //     // Do iteration.
  //     zf = zin + ds * f.constantTerm((zin+z)/2.0);

  //     // Test for convergence.
  //     dz_old = dz;
  //     dz = zf - z;
  //     if(ni){// (we need at least two iterations before testing makes sense)
  //       for(int i=0;i<PSdim;++i){
  //         if(!bounce[i] && (dz[i]==0. || (abs(dz[i])<thresh && abs(dz[i])>=abs(dz_old[i]))))
  //           {bounce[i]=true; ++bcount;}
  //       }
  //     }
  //     z=zf;
  //     ++ni;
  //   }
  //   //std::cerr << "  zf =\n" << zf << std::endl;
  //   //std::cerr << "==> Leaving fixedPointInt2(...)" << std::endl;
  return zf;
}

static Vector fixedPointInt4(const Vector &zin, const VSeries &f, double s, double ds, int nx, int cx)
{
  //   //std::cerr << "==> In fixedPointInt4(zin,f,s,ds,nx,cx) ..." << std::endl;
  //   // Default: nx = 50, cx = 4

  //   // This routine integrates the N-dimensional autonomous differential equation
  //   // z' = f(z) for a distance s, in steps of size ds.  It uses a "Yoshida-fied"
  //   // version of fixedPointInt2 to obtain zf accurate through fourth-order in the
  //   // step-size ds.  The optional arguments nx and cx have the same meaning as in
  //   // implicitInt2().

  //   // Convergence warning flag.
  //   static bool cnvWarn = false;

  //   // The Yoshida constants: 2ya+yb=1; 2ya^3+yb^3=0.
  //   static const double yt=pow(2.,1/3.);
  //   static const double ya=1/(2.-yt);
  //   static const double yb=-yt*ya;

  //   // Initialize accumulated length, current step-size, and number of cuts.
  //   double as=abs(s),st=0., dsc=abs(ds);
  //   if(s<0.) dsc = -dsc;
  //   int ci=0;

  //   // Integrate each step.
  Vector zf=zin;
  //   while(abs(st)<as){
  //     Vector zt;
  //     bool ok=true;
  //     try{
  //       if(abs(st+dsc)>as) dsc=s-st;
  //       zt=fixedPointInt2(zf,f,ya*dsc,nx);
  //       zt=fixedPointInt2(zt,f,yb*dsc,nx);
  //       zt=fixedPointInt2(zt,f,ya*dsc,nx);
  //     }
  //     catch (ConvergenceError &cnverr) {
  //       if(++ci>cx) {
  //         string msg = "Convergence not achieved within " + NumToStr(cx) + " cuts of step-size!";
  //         throw ConvergenceError("ParallelTTracker::fixedPointInt4()",msg);
  //       }
  //       if (!cnvWarn) {
  //         std::cerr << " <***WARNING***> [ParallelTTracker::fixedPointInt4()]:\n"
  //                   << "   Cutting step size, a probable violation of the symplectic condition."
  //                   << std::endl;
  //         cnvWarn=true;
  //       }
  //       dsc*=0.5; ok=false;
  //     }
  //     if(ok){zf=zt; st+=dsc;}
  //   }

  //   //std::cerr << "==> Leaving fixedPointInt4(...)" << std::endl;
  return zf;
}
