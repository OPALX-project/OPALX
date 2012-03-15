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


void ParallelTTracker::visitMonitor(const Monitor &corr)
{
  //   *gmsg << "In Monitor; L= " << corr.getElementLength() << endl;
  myElements.push_back(FieldListType1Entry(dynamic_cast<Monitor*>(corr.clone()),0.,0.));
}


void ParallelTTracker::visitMultipole(const Multipole &mult)
{
  //   *gmsg << "In Multipole; L= " << mult.getElementLength() << endl;
  myElements.push_back(FieldListType1Entry(dynamic_cast<Multipole*>(mult.clone()),0.,0.));
}


void ParallelTTracker::visitRBend(const RBend &bend)
{
  //   *gmsg << "In RBend; L= " << bend.getElementLength() << endl;
  myElements.push_back(FieldListType1Entry(dynamic_cast<RBend*>(bend.clone()),0.,0.));

}


void ParallelTTracker::visitRFCavity(const RFCavity &as)
{
  //   *gmsg << "In RFCavity; L= " << as.getElementLength() << endl;
  Component *elptr = dynamic_cast<RFCavity*>(as.clone());

  if (!elptr->hasAttribute("ELEMEDGE")){
    *gmsg << "RFCavity: no position of the element or no length of the field given!" << endl;
    return;
  }
  double startField = elptr->getAttribute("ELEMEDGE");
  double endField;
  double scaleField = elptr->getAttribute("VOLT");

  elptr->readFieldMap(startField, endField, 1.0);
  myElements.push_back(FieldListType1Entry(elptr,startField, endField));
}

void ParallelTTracker::visitTravelingWave(const TravelingWave &as)
{
  //   *gmsg << "In RFCavity; L= " << as.getElementLength() << endl;
  Component *elptr = dynamic_cast<TravelingWave*>(as.clone());
  if (!elptr->hasAttribute("ELEMEDGE")){
    *gmsg << "RFCavity: no position of the element or no length of the field given!" << endl;
    return;
  }
  double startField = elptr->getAttribute("ELEMEDGE");
  double endField;
  double scaleField = elptr->getAttribute("VOLT");

  elptr->readFieldMap(startField, endField, 1.0);
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
  Component *elptr = dynamic_cast<Solenoid*>(solenoid.clone());

  if (!elptr->hasAttribute("ELEMEDGE"))
    {
      *gmsg << "Solenoid: no position of the element given!" << endl;    
      return;
    }
  
  double startField = elptr->getAttribute("ELEMEDGE");
  double endField;

  elptr->readFieldMap(startField, endField, 1.0);

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

}
// 2007/04/19 CKR
void ParallelTTracker::execute(){
  // prepare the tracker for doing its job
  // prepare the elements of the line to be tracked, i.e. load fieldmaps and so on
  //  itsBeamline->arrange();

  double t = itsBunch->getT(), dt = itsBunch->getdT();

  /// in case of a gun save initial timestep
  double dtSave = itsBunch->getdT();

  double tEmission = itsBunch->getTEmission();
  int gunSubTimeSteps = 10;

  Vector_t um, a, s, externalE, externalB;
  bool EndOfLineReached;

  double tmp;
  double recpgamma;
  int emissionSteps = 0;
  int step = 0;

  bool partOutOfBounds;

  bool doWakeField = (Options::wakeCalcStep != -1);

  if(OPAL.inRestartRun())
    {
      step = OPAL.getRestartStep()+1;
    }
  else if ( itsBunch->doEmission()) 
    {
      emissionSteps = (int)itsBunch->pbin_m->getNBins()*gunSubTimeSteps;
      
      *gmsg << "Do emission for " << tEmission << " [s] using " << itsBunch->pbin_m->getNBins() << " energy bins " << endl;
      *gmsg << "Change dT from " <<  itsBunch->getdT() << " [s] to "; 
      
      //x    itsBunch->setdT(tEmission/emissionSteps);

      *gmsg <<  itsBunch->getdT() << " [s] during emission " << endl;; 
    }

  *gmsg << "executing ParallelTTracker, initial DT " << itsBunch->getdT() << " [s] max integration steps " << maxSteps_m << endl;
  //  *gmsg << "the mass is: " << itsReference.getM() << ", its charge: " << itsReference.getQ() << endl;
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
  // `use scale factor to get dimensionless variables
  // OBSOLETE?: scaling is now down with particel dt
  // itsBunch->R /= Vector_t(scaleFactor_m);
  // itsBunch->calcBeamParameters();
  
  // TODO: rescale P = beta_{x,y,z} gamma IF using a distribution
  // (NOT here)

  if (!(itsBunch->weHaveBins())) {
    IpplTimings::startTimer(BinRepartTimer_m);
    itsBunch->do_binaryRepart();
    IpplTimings::stopTimer(BinRepartTimer_m);
    Ippl::Comm->barrier();
  }

  size_t nEmit = 0;
  Inform msg ("bin- ",INFORM_ALL_NODES);


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
    *gmsg << "MINSTEPFORREBIN " << minBinEmitted << endl;
  } 

  for(step; step < maxSteps_m; step++){
    EndOfLineReached = true;  //check if any particle hasn't reached the end of the field from the last element

    IpplTimings::startTimer(timeIntegrationTimer1_m);
    /* 
       transport and emit particles 
       that would pass the cathode in the
       first half-step

       to make IPPL and the field solver happy
       make sure that at least 10 particles are emitted

       also remember that node 0 has 
       all the particles to be emitted
    */
    if ((itsBunch->weHaveBins())) { 
      if ((itsBunch->getTotalNum() == 0)) {
        while (nEmit < minBinEmitted) {
          size_t ne = 0;
          if (Ippl::myNode() == 0) {
            for (int b=0; b<itsBunch->getNumBins(); b++) 
              ne += itsBunch->emitFirstHalfStep(b);
          }
          step++;
          reduce(ne,ne,OpAddAssign());
          nEmit += ne;
        }
      }
      else {
        if (Ippl::myNode() == 0) {
          size_t ne = 0;
          for (int b=0; b<itsBunch->getNumBins(); b++) 
            ne += itsBunch->emitFirstHalfStep(b);
        }
      }
      itsBunch->updateBinStructure(); // we need this because only node 0 is emitting

      if (step>minStepforReBin) { 
	// check if we should rebin
	double maxdE = itsBunch->getMaxdEBins();
	if (maxdE < itsBunch->getRebinEnergy()) {
	  *gmsg << "**********************************************************" << endl;
	  *gmsg << "maxdE < " << maxdE << " [keV] we rebin" << endl;
	  *gmsg << "**********************************************************" << endl;
	  itsBunch->rebin();
	}
      }
    }
    itsBunch->boundp();

        
    //reset E and B to Vector_t(0.0) for every step
    itsBunch->Ef = Vector_t(0.0);
    itsBunch->Bf = Vector_t(0.0);

    //IFF: inc global simulation time at the end of the loop
    //simulation time uses the global big timestep
    //t += 0.5 * itsBunch->getdT(); 

    for (int i = 0; i < itsBunch->getLocalNum(); i++) {
      //scale each particle with own dt
      itsBunch->R[i] /= Vector_t(Physics::c * itsBunch->dt[i]);
      
      recpgamma = 1.0 / sqrt(1.0 + dot(itsBunch->P[i],itsBunch->P[i]));
      
      /** \f[ \vec{x}_{n+1/2} = \vec{x}_{n} + \frac{1}{2}\vec{v}_{n-1/2}\quad (= \vec{x}_{n} + \frac{\Delta t}{2} \frac{\vec{\beta}_{n-1/2}\gamma_{n-1/2}}{\gamma_{n-1/2}}) \f]
       *
       * \code
       * R[i] += 0.5 * P[i] * recpgamma; 
       * \endcode
       */
      itsBunch->R[i] += 0.5 * recpgamma * itsBunch->P[i];
    }

    itsBunch->boundp();
    
    IpplTimings::stopTimer(timeIntegrationTimer1_m);

    itsBunch->calcBeamParameters();

    if(itsBunch->hasFieldSolver() && itsBunch->getTotalNum() > 100) {

      if (itsBunch->getTotalNum() > 1000) {
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
      //itsBunch->boundp();
    }

    IpplTimings::startTimer(timeIntegrationTimer2_m);

    //FIXME: rethink scaling!
    for (int i=0; i < itsBunch->getLocalNum(); i++)
      itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->dt[i]);
    
    /* 
       transport and emit particles 
       that would pass the cathode in the
       second half-step

       also remember that node 0 has 
       all the particles to be emitted
    */
    if (itsBunch->weHaveBins()) {
      if (Ippl::myNode() == 0) {
        size_t ne = 0;
        for (int b=0; b<itsBunch->getNumBins(); b++) 
          ne +=  itsBunch->emitSecondHalfStep(b);
      }
      itsBunch->boundp();
      itsBunch->updateBinStructure(); // we need this because only node 0 is emitting
    }

    //for (int i=0; i < itsBunch->getLocalNum(); i++)
    //itsBunch->R[i] /= Vector_t(Physics::c * itsBunch->dt[i]);

    size_t totalParticles_i = itsBunch->getTotalNum();

    /*
      This is a experimental implementation of Wake Fields
      Here we assume that we are in an RFCavity and many 
      parameters are hard wired.
    */

    if (doWakeField) {
      IpplTimings::startTimer(WakeFieldTimer_m);
      itsBunch->calcLineDensity();
      /*
	Concolution

      */
      IpplTimings::stopTimer(WakeFieldTimer_m);

    }

    for (int i = 0; i < itsBunch->getLocalNum(); i++) {
      
      //dont push particles that have just been emitted
      if(itsBunch->R[i](2) == 0.0)
        continue;
    
      partOutOfBounds = false;
      //reset helper vectors in each step since fields are added in getFieldstrenght
      externalB = Vector_t(0.0);
      externalE = Vector_t(0.0);

      //when calculating external Ef and Bf use dimensions
      //itsBunch->R[i] *= Vector_t(scaleFactor_m);
      //itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->dt[i]);

      IpplTimings::startTimer(timeFieldEvaluation_m);

      if ( itsBunch->R[i](2) < myFieldList.back().End )
        {
          EndOfLineReached = false;
          if (itsBunch->R[i](2) >= myFieldList[itsBunch->LastSection[i]].Start)
            {
              //               if (itsBunch->R[i](2) >= myFieldList[itsBunch->LastSection[i]].End)
              while (itsBunch->R[i](2) >= myFieldList[itsBunch->LastSection[i]].End)
                itsBunch->LastSection[i] += 1;
              //               if (i == 0)
              //                 cout << itsBunch->R[i] << "\t";
              for (ElementIterator el_it = myFieldList[itsBunch->LastSection[i]].Elements.begin(); el_it != myFieldList[itsBunch->LastSection[i]].Elements.end(); el_it++)
                {
		  if ((i == 0) && ((*el_it)->getName() == string("LRF0"))){
		    IpplTimings::startTimer(WakeFieldTimer_m);		    
		    *gmsg << "applt wake field tp particle " << (*el_it)->getName() << endl;
		    IpplTimings::stopTimer(WakeFieldTimer_m);
		  }
                  partOutOfBounds = partOutOfBounds || (*el_it)->getFieldstrength(itsBunch->R[i], t+itsBunch->dt[i]/2.0, externalE, externalB);
                }
              //               if (i == 0)
              //                 cout << externalE << endl;
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
      
      //itsBunch->R[i] /= Vector_t(scaleFactor_m);
      itsBunch->R[i] /= Vector_t(Physics::c * itsBunch->dt[i]);
      
      itsBunch->Ef[i] += externalE;
      itsBunch->Bf[i] += externalB;

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
      um = itsBunch->P[i] + 0.5 * itsReference.getQ() * itsBunch->dt[i] /*itsBunch->getdT()*/ / itsReference.getM() * c * itsBunch->Ef[i];

      /** \f[ \vec{a}_{n} = \frac{Q \; \Delta t}{2\gamma M c^2}\vec{B} \f]
       * \code
       * recpgamma = 1.0 / sqrt(1.0 + um[0]*um[0] + um[1]*um[1] + um[2]*um[2]); 
       * 
       * tmp = 0.5 * Q * recpgamma * dt/M  * c * c;
       * a = tmp * B[i];
       * \endcode
       */
      recpgamma = 1.0 / sqrt(1.0 + dot(um,um));
      
      tmp = 0.5 * itsReference.getQ() * recpgamma * itsBunch->dt[i] /*itsBunch->getdT()*/ / itsReference.getM() * c * c;
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

      itsBunch->P[i] = um + 0.5 * itsReference.getQ()  * itsBunch->dt[i] /*itsBunch->getdT()*/ / itsReference.getM() * c * itsBunch->Ef[i]; 

      recpgamma = 1.0 / sqrt(1.0 + dot(itsBunch->P[i],itsBunch->P[i])); 

      /** \f[ \vec{x}_{n+1} = \vec{x}_{n+1/2} + \frac{1}{2}\vec{v}_{n+1/2}\quad (= \vec{x}_{n+1/2} + \frac{\Delta t}{2} \frac{\vec{\beta}_{n+1/2}\gamma_{n+1/2}}{\gamma_{n+1/2}}) \f] 
       * \code
       * R[i] += 0.5 * P[i] * recpgamma; 
       * \endcode
       */
      itsBunch->R[i] += 0.5 * itsBunch->P[i] * recpgamma;

      //and scale back to dimensions
      itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->dt[i]); 
      
      //reset time step if particle was emitted in the first half-step
      //the particle is now in sync with the simulation timestep
      if(itsBunch->R[i](2) > 0.0 && itsBunch->dt[i] != itsBunch->getdT()) 
        itsBunch->dt[i] = itsBunch->getdT();

    } // end particle loop

    IpplTimings::stopTimer(timeIntegrationTimer2_m);

    itsBunch->boundp();

    /**
       Delete marked particles
    */
    for (int i = 0; i < itsBunch->getLocalNum(); i++) 
      if(itsBunch->Bin[i] < 0)
      	itsBunch->destroy(1,i);
    itsBunch->update();
    
    size_t totalParticles_f = itsBunch->getTotalNum();

    if(totalParticles_i != totalParticles_f) {
      *gmsg << "======================================================================" << endl;
      *gmsg << "Particle loss at t= " << t << " n= " << totalParticles_i - totalParticles_f << " Ntot= " << itsBunch->getTotalNum() << endl;  
      *gmsg << "======================================================================" << endl;
    } 
   
    t += /*0.5*/ itsBunch->getdT(); //t after a full global timestep with dT "synchronization point" for simulation time
    
    itsBunch->setT(t);

    //IFF: cheap step dump regulation
    if ( totalParticles_f > 0 && step % Options::psDumpFreq == 0 ) {
      //itsBunch->R *= Vector_t(scaleFactor_m);
      //itsBunch->boundp();
      double sposRef = itsBunch->get_sPos();
      if (isnan(sposRef))
        *gmsg << "Step " << step << " no emission yet "  << " t= " << itsBunch->getdT()*step << " [s]" << endl;
      else
        *gmsg << "Step " << step << " at " << sposRef << " [m] t= " << itsBunch->getdT()*step << " [s]" << endl;
      *gmsg << *itsBunch << endl;

      //FDext = {BHead, EHead, BRef, ERef, BTail, ETail}
      Vector_t rmin, rmax, FDext[6], extE, extB;
      double sposHead, sposTail; //, meanEnergy, energy[itsBunch->getTotalNum()];
      //float zCenter;

      for(int k=0; k<6; k++) FDext[k] = Vector_t(0.0,0.0,0.0);
      //sample fields at (0,0,rmin), (0,0,rmax) and the centroid location
      itsBunch->get_bounds(rmax,rmin);

      sposHead = rmax(2);
      sposTail = rmin(2);

      Vector_t pos[3];
      pos[0] = Vector_t(0.002,0.002,sposHead);
      pos[1] = Vector_t(0.002,0.002,sposRef);
      pos[2] = Vector_t(0.002,0.002,sposTail);

      for(int k=0; k < 3; k++) {
        extE = Vector_t(0,0,0);
        extB = Vector_t(0,0,0);

        if ( pos[k](2) <= myFieldList.back().End && pos[k](2) >= myFieldList.front().Start )
          for (FieldListType2Iterator fl_it = myFieldList.begin(); fl_it != myFieldList.end(); fl_it++)
            {
              if (pos[k](2) >= (*fl_it).Start && pos[k](2) < (*fl_it).End)
                {
                  for (ElementIterator el_it = (*fl_it).Elements.begin(); el_it != (*fl_it).Elements.end(); el_it++)
                    (*el_it)->getFieldstrength(pos[k], t, extE, extB);
                  break;
                }
            }
	
        FDext[2*k]   = extB;
        FDext[2*k+1] = extE;
      }
      if (!(isnan(sposRef)))      
        itsDataSink->writePhaseSpace(*itsBunch, FDext, sposHead, sposRef, sposTail);

      //itsBunch->R /= Vector_t(scaleFactor_m);
    } else {
      //itsBunch->R *= Vector_t(scaleFactor_m);
      
      //FDext = {BHead, EHead, BRef, ERef, BTail, ETail}
      if (isnan(itsBunch->get_sPos()))
        *gmsg << "Step " << step << " no emission yet "  << " t= " << itsBunch->getdT()*step << " [s]" << endl;
      else
        *gmsg << "Step " << step << " at " << itsBunch->get_sPos() << " [m] t= " << itsBunch->getdT()*step << " [s]" << endl;
      //itsBunch->R /= Vector_t(scaleFactor_m);
    }
    
    //TEST the non-len reduce: reduce(&EndOfLineReached, &EndOfLineReached, OpBitwiseOrAssign());
    reduce(&EndOfLineReached, &EndOfLineReached + 1, &EndOfLineReached, OpBitwiseOrAssign());
    
    if(step > emissionSteps) 
      if (EndOfLineReached) break;
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
