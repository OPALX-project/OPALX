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
  allElements(),
  allSections()
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
  allElements(),
  allSections(),
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
  for (FieldListIterator compindex = allElements.begin(); compindex != allElements.end(); ++compindex)
    {
      delete (*compindex).Element;
    }
  allElements.clear();
  allSections.clear();
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
  *gmsg << "BeamBeam not implemented yet!" << endl;
  //  *gmsg << "In BeamBeam; "<< endl;
}


void ParallelTTracker::visitCollimator(const Collimator &coll)
{
  *gmsg << "Collimator not implemented yet!" << endl;
  //  *gmsg << "In Collimator; L= " << coll.getElementLength() << endl;
  //   allElements.push_back(FieldListEntry(dynamic_cast<Collimator*>(coll.clone()),0.,0.));
}


void ParallelTTracker::visitCorrector(const Corrector &corr)
{
  *gmsg << "Corrector not implemented yet!" << endl;
  //    *gmsg << "In Corrector; L= " << corr.getElementLength() << endl;
  //   allElements.push_back(FieldListEntry(dynamic_cast<Corrector*>(corr.clone()),0.,0.));
}


void ParallelTTracker::visitDiagnostic(const Diagnostic &diag)
{
  //     *gmsg << "In Diagnostic; L= " << diag.getElementLength() << endl;
  //   allElements.push_back(FieldListEntry(dynamic_cast<Diagnostic*>(diag.clone()),0.,0.));
}


void ParallelTTracker::visitDrift(const Drift &drift)
{
  //   *gmsg << "In drift L= " << drift.getElementLength() << endl;
  Drift *elptr = dynamic_cast<Drift*>(drift.clone()->removeWrappers());
  if (!elptr->hasAttribute("ELEMEDGE"))
    {
      *gmsg << "Drift: no position of the element given!" << endl;
      return;
    }
  double startField = elptr->getAttribute("ELEMEDGE");
  double endField;
  elptr->initialise(itsBunch, startField, endField, 1.0);
  allElements.push_back(FieldListEntry(elptr, startField, endField));

}


void ParallelTTracker::visitLambertson(const Lambertson &lamb)
{
  *gmsg << "Lambertson not implemented yet!" << endl;
  //     *gmsg << "In Lambertson; L= " << lamb.getElementLength() << endl;
  //   allElements.push_back(FieldListEntry(dynamic_cast<Lambertson*>(lamb.clone()),0.,0.));
}


void ParallelTTracker::visitMarker(const Marker &marker)
{
  //     *gmsg << "In Marker; L= " << marker.getElementLength() << endl;
  //   allElements.push_back(FieldListEntry(dynamic_cast<Marker*>(marker.clone()),0.,0.));
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
  allElements.push_back(FieldListEntry(elptr,startField, endField));
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
  allElements.push_back(FieldListEntry(elptr,startField, endField));
}


void ParallelTTracker::visitRBend(const RBend &bend)
{
  RBend *elptr = dynamic_cast<RBend*>(bend.clone()->removeWrappers());
  if (!elptr->hasAttribute("ELEMEDGE"))
    {
      *gmsg << "RBend: no position of the element given!" << endl;
      return;
    }

  if (bend.getElementLength() < 1e-3)
    {
      *gmsg << "RBend: no length of the element given!" << endl;
      return;
    }

  double startField = elptr->getAttribute("ELEMEDGE");
  double endField;
  elptr->setLength(bend.getElementLength());
  elptr->initialise(itsBunch, startField, endField, 1.0);
  allElements.push_back(FieldListEntry(elptr,startField, endField));

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
  allElements.push_back(FieldListEntry(elptr,startField, endField));
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
  allElements.push_back(FieldListEntry(elptr,startField, endField));
}


void ParallelTTracker::visitRFQuadrupole(const RFQuadrupole &rfq)
{
  *gmsg << "RFQuadrupole not implemented yet!" << endl;
  //     *gmsg << "In RFQuadrupole; L= " << rfq.getElementLength() << endl;
  //   allElements.push_back(FieldListEntry(dynamic_cast<RFQuadrupole*>(rfq.clone()),0.,0.));
}

void ParallelTTracker::visitSBend(const SBend &bend)
{
  *gmsg << "SBend not implemented yet!" << endl;
  //     *gmsg << "In SBend; L= " << bend.getElementLength() << endl;
  //   allElements.push_back(FieldListEntry(dynamic_cast<SBend*>(bend.clone()),0.,0.));
}


void ParallelTTracker::visitSeparator(const Separator &sep)
{
  *gmsg << "Separator not implemented yet!" << endl;
  //  *gmsg << "In Separator; L= " << sep.getElementLength() << endl;
  //   allElements.push_back(FieldListEntry(dynamic_cast<Separator*>(sep.clone()),0.,0.));

}


void ParallelTTracker::visitSeptum(const Septum &sept)
{
  *gmsg << "Septum not implemented yet!" << endl;
  //  *gmsg << "In Septum; L= " << sept.getElementLength() << endl;
  //   allElements.push_back(FieldListEntry(dynamic_cast<Septum*>(sept.clone()),0.,0.));
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

  allElements.push_back(FieldListEntry(elptr, startField, endField));
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
   *  of the elements start and end respectively. The field of an element can be contained
   *  in one (at least) or more sections.
   */

  list<double> StartEnd;
  list<Component*> tmp;
  FieldListIterator fld_it;
  const FieldListIterator fldend_it = allElements.end();
  SectionListIterator sec_it, secend_it;
  list<double>::iterator pos_it, next_it, last_it;
  double tolerance = 1.e-8;

  /* there might be elements with length zero or extremely short ones.
     we delete them such that they don't appear in the simulation
  */
  allElements.sort(FieldListEntry::SortAscByStart);

  for (fld_it = allElements.begin();
       fld_it != fldend_it;
       ++fld_it)
    {
      while (fabs((*fld_it).End - (*fld_it).Start) < tolerance && fld_it != fldend_it)
        {
          FieldListIterator temp = fld_it;
          ++temp;
          allElements.erase(fld_it);
          fld_it = temp;
        }
      if (fld_it != fldend_it)
        {
          StartEnd.push_back((*fld_it).Start);
          StartEnd.push_back((*fld_it).End);
        }
    }
  StartEnd.sort();

  next_it = StartEnd.begin(); ++next_it;
  last_it = StartEnd.end();
  for (pos_it = StartEnd.begin();
       next_it != last_it;
       ++pos_it, ++next_it)
    if (*next_it - *pos_it < tolerance)
      *next_it = *pos_it;

  StartEnd.unique();  // remove duplicate entries

  next_it = StartEnd.begin(); ++next_it;
  last_it = StartEnd.end();
  for (pos_it = StartEnd.begin();
       next_it != last_it;
       ++pos_it, ++next_it)
    {
      tmp.clear();
      for (fld_it = allElements.begin();
           fld_it != fldend_it;
           ++fld_it)
        {
          if ((*fld_it).Start <= *pos_it  + tolerance
              &&
              (*fld_it).End   >= *next_it - tolerance)
            tmp.push_back((*fld_it).Element);
          else if ((*fld_it).Start >= *next_it)
            break;
        }
      if (tmp.size() > 0)
        allSections.push_back(SectionListEntry(tmp,*pos_it,*next_it));
    }

  secend_it = allSections.end();
  *gmsg << "--- BEGIN FIELD LIST ----------------------------------------------------------------------\n" << endl;
  for (sec_it = allSections.begin();
       sec_it != secend_it;
       ++sec_it)
    {
      *gmsg << "--- " << (*sec_it).Start << " m -- " << (*sec_it).End << " m -- alpha = " << (*sec_it).Orientation(0) << " -- beta = " << (*sec_it).Orientation(1) << "---------------------------\n";

      for (ElementIterator el_it = (*sec_it).Elements.begin();
           el_it != (*sec_it).Elements.end();
           ++el_it)
        *gmsg << (*el_it)->getName() << '\n';

    }
  *gmsg << "--- END   FIELD LIST ----------------------------------------------------------------------\n" << endl;


}

// 2007/04/19 CKR
void ParallelTTracker::execute()
{
  double tmp;
  double recpgamma, gamma;
  double t = itsBunch->getT();
  double dt = itsBunch->getdT();
  double tEmission = itsBunch->getTEmission();
  BorisPusher pusher(itsReference);

  int gunSubTimeSteps = 10;
  int emissionSteps = 0;


  Vector_t um, a, s;
  Vector_t externalE, externalB;
  Vector_t rmin, rmax, rtmp;
  const Vector_t vscaleFactor = Vector_t(scaleFactor_m);

  bool EndOfLineReached;
  bool partOutOfBounds;
  bool bends;               // flag which indicates whether any particle is within the influence of bending element.
                            // if this is the case we track the reference particle as if it were a real particle,
                            // otherwise the reference particle is defined as the centroid particle of the bunch

  bool hasWake = false;     // flag which indicates whether any particle is within the influence of a wake field


  size_t totalParticles_i = itsBunch->getTotalNum();

  WakeFunction *wf;

  FieldListIterator fld_it;
  FieldListIterator fldend_it;

  long long step;
  if (OPAL.inRestartRun())
    step = OPAL.getRestartStep() + 1;
  else if (OPAL.hasBunchAllocated())
    step = 1;
  else
    step = lround(t/dt);

  if (itsBunch->doEmission())
    {
      emissionSteps = (int)itsBunch->pbin_m->getNBins()*gunSubTimeSteps;
      *gmsg << "Do emission for " << tEmission << " [s] using " << itsBunch->pbin_m->getNBins() << " energy bins " << endl
            << "Change dT from " <<  itsBunch->getdT() << " [s] to "
            <<  itsBunch->getdT() << " [s] during emission " << endl;;

      if (Ippl::myNode() == 0)
        {
          // if we don't want opal to behave like impactt move the bunch to the cathode;
          // otherwise comment the next line out;
          RealVariable *ar = dynamic_cast<RealVariable *>(OPAL.find("TSHIFT"));
          if (ar)
            {
              t = itsBunch->calcTimeDelay(ar->getReal());
              *gmsg << "The time delay is " << t << " [s]" << endl;
              itsBunch->setT(t);
            }
        }
    }

  // set dt for all particles already in the simulation,
  // i.e. when doing a restarted simulation
  for(int i=0; i < itsBunch->getLocalNum(); ++i)
    {
      itsBunch->dt[i] = itsBunch->getdT();
    }

  *gmsg << "executing ParallelTTracker, initial DT " << itsBunch->getdT()
        << " [s]; max integration steps " << maxSteps_m << " step= " << step << endl
        << "the mass is: " << itsReference.getM()*1e-6 << " MeV, its charge: " << itsReference.getQ() << endl;

  itsBeamline->accept(*this);  // fill list allElements
  buildupFieldList();          // fill list allSections

  fldend_it = allElements.end();

  if (OPAL.inRestartRun())   // the bunch is in the middle of the lattice.
                             // therefore assign to each particle the number of the section in which it is at the moment.
                             // THIS IS PROBABELY NOT NEEDED ANYMORE SINCE THIS IS ALSO DONE BEFORE WE EVALUATE THE FIELDS
    {
      for (int i = 0; i < itsBunch->getLocalNum(); ++i)
        {
          itsBunch->LastSection[i] = 0;
          for (int l = 0; l < allSections.size(); ++l)
            {
              if (itsBunch->R[i](2) >= allSections[l].Start
                  &&
                  itsBunch->R[i](2) <= allSections[l].End)
                {
                  itsBunch->LastSection[i] = l;
                  break;
                }
            }
        }
    }


  if (!(itsBunch->weHaveBins()))
    {
      IpplTimings::startTimer(BinRepartTimer_m);
      itsBunch->do_binaryRepart();
      IpplTimings::stopTimer(BinRepartTimer_m);
      Ippl::Comm->barrier();
    }

  // Check if there are any particles in simulation. If there are,
  // as in a restart, use the usual function to calculate beam
  // parameters. If not, calculate beam parameters of the initial
  // beam distribution.
  if (totalParticles_i == 0)
    {
      itsBunch->calcBeamParametersInitial();
    }
  else
    {
      itsBunch->calcBeamParameters();
    }

 // if (!OPAL.hasBunchAllocated()) // update the orientation of the bunch such that
                                 // vec{p} = (0,0,p_z) and vec{r} = (0,0,r_z)
    {
      RefPartR_suv_m = RefPartR_zxy_m = itsBunch->get_rmean();
      RefPartP_suv_m = RefPartP_zxy_m = itsBunch->get_pmean();

      updateSpaceOrientation(true);
    }

  RefPartR_suv_m = RefPartR_zxy_m = itsBunch->get_rmean();
  RefPartP_suv_m = RefPartP_zxy_m = itsBunch->get_pmean();

  // set local coordinate system
  for (int i = 0; i < itsBunch->getLocalNum(); ++i)
    itsBunch->ResetLocalCoordinateSystem(i, Vector_t(0.0), 0.0);


  /* Activate all elements which influence the particles when the simulation starts;
   *  mark all elements which are already past.
   */
  itsBunch->get_bounds(rmin,rmax);
  double margin = 3. * RefPartP_suv_m(2) * scaleFactor_m / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
  margin = 0.01 > margin? 0.01: margin;
  for (fld_it = allElements.begin(); fld_it != fldend_it; ++fld_it)
    {
      if (!(*fld_it).goneLive
          &&
          rmax(2) + margin > (*fld_it).Start
          &&
          (rmin(2) - margin) < (*fld_it).End)
        {
          (*fld_it).Element->goOnline();
          *gmsg << (*fld_it).Element->getName() << " gone live at step #" << step << " rmax = " << rmax(2) << endl;
          (*fld_it).goneLive = true;
          ///// TODO:
          ///// if field bends inform the element about the current orientation of the reference particle
          ///// update the local coordinates of all particles
          /////
        }
      if ((rmin(2) - margin) > (*fld_it).End)
        {
          (*fld_it).goneLive = true;
          (*fld_it).goneOff = true;
        }
    }

  double minBinEmitted  = 10.0;
  RealVariable *ar = dynamic_cast<RealVariable *>(OPAL.find("MINBINEMITTED"));
  if (ar)
    {
      minBinEmitted = ar->getReal();  // the space charge solver crashes if we use less than ~10 particles.
                                      // This variable controls the number of particles to be emitted before we use
                                      // the space charge solver.
      *gmsg << "MINBINEMITTED " << minBinEmitted << endl;
    }

  double minStepforReBin  = 200.0;
  RealVariable *br = dynamic_cast<RealVariable *>(OPAL.find("MINSTEPFORREBIN"));
  if (br)
    {
      minStepforReBin = br->getReal();  // this variable controls the minimal number of steps of emission (using bins)
                                        // before we can merge the bins
      *gmsg << "MINSTEPFORREBIN " << minStepforReBin << endl;
    }

  
  int repartFreq = 1000;

  for(step; step < maxSteps_m; ++step)
    {
      EndOfLineReached = true;  //check if any particle hasn't reached the end of the field from the last element
      bends = false;
      hasWake = false;
      for (SectionListIterator sec_it = allSections.begin(); sec_it != allSections.end(); ++sec_it)
        (*sec_it).live = false;

          IpplTimings::startTimer(timeIntegrationTimer1_m);

      //reset E and B to Vector_t(0.0) for every step
      itsBunch->Ef = Vector_t(0.0);
      itsBunch->Bf = Vector_t(0.0);

      for (int i = 0; i < itsBunch->getLocalNum(); ++i)
        {
          //scale each particle with c*dt
          itsBunch->R[i] /= vscaleFactor;

          pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);

          // update local coordinate system of particle
          itsBunch->X[i] /= vscaleFactor;
          pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i], allSections[itsBunch->LastSection[i]].Orientation), itsBunch->getdT());
          itsBunch->X[i] *= vscaleFactor;
        }

      if(totalParticles_i > minBinEmitted)
        itsBunch->boundp();

      IpplTimings::stopTimer(timeIntegrationTimer1_m);

      itsBunch->calcBeamParameters();

      /** \f[ Space Charge  \f]
      */
      if(itsBunch->hasFieldSolver() && totalParticles_i > minBinEmitted && fabs(itsBunch->getChargePerParticle()) > 0.0)
        {
    	  // Do repartition if we have enough particles.
          if (totalParticles_i > 1000 && ((step%repartFreq) == 0) )
            {
	      *gmsg << "**********************************************************" << endl;
	      *gmsg << "do repartition because of repartFreq" << endl;
	      *gmsg << "**********************************************************" << endl;
              IpplTimings::startTimer(BinRepartTimer_m);
              itsBunch->do_binaryRepart();
              IpplTimings::stopTimer(BinRepartTimer_m);
              Ippl::Comm->barrier();
	      *gmsg << "**********************************************************" << endl;
	      *gmsg << "do repartition done" << endl;
	      *gmsg << "**********************************************************" << endl;
            }

          // Calculate space charge.
          if (itsBunch->weHaveBins())
        	// When we have energy bins.
            {
              itsBunch->calcGammas();
              for (int binNumber = 0; binNumber <= itsBunch->getLastemittedBin() && binNumber < itsBunch->getNumBins(); ++binNumber)
                {
                  itsBunch->setBinCharge(binNumber, itsBunch->getChargePerParticle());
                  itsBunch->computeSelfFields(binNumber);
                }
              // the next line is not needed anymore
              itsBunch->Q = itsBunch->getChargePerParticle();
            }
          else
        	// When we don't.
            {
              itsBunch->computeSelfFields();
            }
        }

      IpplTimings::startTimer(timeIntegrationTimer2_m);

      /*
         transport and emit particles
         that passed the cathode in the first
         half-step or that would pass it in the
         second half-step.

         to make IPPL and the field solver happy
         make sure that at least 10 particles are emitted

         also remember that node 0 has
         all the particles to be emitted

         this has to be done *after* the calculation of the
         space charges! thereby we neglect space charge effects
         in the very first step of a new-born particle.
      */

      if ((itsBunch->weHaveBins()))
        {
          int ne = 0;
          if (Ippl::myNode() == 0)
            {
              for (int b=0; b<itsBunch->getNumBins(); ++b)
                {
                  ne += itsBunch->emitParticles(b);
                }

            }
          reduce(ne,ne,OpAddAssign());
          totalParticles_i += ne;

          itsBunch->updateBinStructure(); // we need this because only node 0 is emitting
          if (step > minStepforReBin)
            {
              itsBunch->calcGammas();
              const double maxdE = itsBunch->getMaxdEBins();

              if (maxdE < itsBunch->getRebinEnergy())
                {
                  *gmsg << "**********************************************************" << endl;
                  *gmsg << "maxdE < " << maxdE << " [keV] we rebin" << endl;
                  *gmsg << "**********************************************************" << endl;
                  itsBunch->rebin();
                  *gmsg << "**********************************************************" << endl;
                  *gmsg << "do repartition " << endl;
                  *gmsg << "**********************************************************" << endl;
		  IpplTimings::startTimer(BinRepartTimer_m);
		  itsBunch->do_binaryRepart();
		  IpplTimings::stopTimer(BinRepartTimer_m);
		  Ippl::Comm->barrier();
		  *gmsg << "**********************************************************" << endl;
                  *gmsg << "do repartition done " << endl;
                  *gmsg << "**********************************************************" << endl;
		}
            }
        }

      // push the reference particle by a half step
      recpgamma = 1.0 / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
      RefPartR_zxy_m += RefPartP_zxy_m * recpgamma / 2. * scaleFactor_m;

      //
      // get external fields for all particles
      //
      IpplTimings::startTimer(timeFieldEvaluation_m);
      for (int i = 0; i < itsBunch->getLocalNum(); ++i)
        {
          //FIXME: rethink scaling!
          itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i]);
	  //          itsBunch->R[i] *= vscaleFactor;
          partOutOfBounds = false;

          //reset helper vectors in each step since fields are added in getFieldstrength
          externalB = Vector_t(0.0);
          externalE = Vector_t(0.0);

          if ( itsBunch->R[i](2) < allSections.back().getEnd(itsBunch->R[i](0), itsBunch->R[i](1)))
            {
              EndOfLineReached = false;
              int lastSection = itsBunch->LastSection[i];
              while (itsBunch->R[i](2) >= allSections[lastSection].getEnd(itsBunch->R[i](0), itsBunch->R[i](1)))
                ++lastSection;
              if (lastSection > itsBunch->LastSection[i])
                {
                  itsBunch->LastSection[i] = lastSection;
                  itsBunch->ResetLocalCoordinateSystem(i,allSections[lastSection].Orientation, allSections[lastSection].Start);
                }

              allSections[lastSection].live = true;

              if (itsBunch->R[i](2) >= allSections[lastSection].getStart(itsBunch->R[i](0), itsBunch->R[i](1)))
                {
                  ElementIterator last_el_it = allSections[lastSection].Elements.end();
                  bends = bends     || allSections[lastSection].bends; // check if any element bends the beam; make sure that this is communicated to all processors!

                  if (!hasWake && allSections[lastSection].hasWake)
                    {
                      hasWake = true;
                      wf = allSections[lastSection].Wake; // DOES NOT WORK IN PARALLEL!
                    }
                  for (ElementIterator el_it = allSections[lastSection].Elements.begin(); el_it != last_el_it; ++el_it)
                   {
                      partOutOfBounds = partOutOfBounds || (*el_it)->apply(i, t+itsBunch->dt[i]/2.0, externalE, externalB);
                    }
                }
            }

          // skip rest of the particle push if the
          // particle is out of bounds i.e. does not see
          // a E or B field
          if (partOutOfBounds)
            {
              itsBunch->Bin[i] = -1;
              continue;
            }

          if (fabs(allSections[itsBunch->LastSection[i]].Orientation(0)) > 1e-6
              ||
              fabs(allSections[itsBunch->LastSection[i]].Orientation(1)) > 1e-6)
            {
//               externalE = TransformTo(externalE, allSections[itsBunch->LastSection[i]].Orientation);
//               externalB = TransformTo(externalB, allSections[itsBunch->LastSection[i]].Orientation);
              externalE = TransformBack(externalE, allSections[itsBunch->LastSection[i]].Orientation);
              externalB = TransformBack(externalB, allSections[itsBunch->LastSection[i]].Orientation);
            }

          itsBunch->Ef[i] += externalE;
          itsBunch->Bf[i] += externalB;

          itsBunch->R[i] /= Vector_t(Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i]);
	  //          itsBunch->R[i] /= vscaleFactor;
        }

      reduce(bends,bends,OpOr());
      reduce(hasWake,hasWake,OpOr());

      IpplTimings::stopTimer(timeFieldEvaluation_m);

      if (itsBunch->getLocalNum() == 0)
        EndOfLineReached = false;

      /**
         Delete marked particles
      */
      int ne = 0;
      for (int i = 0; i < itsBunch->getLocalNum(); ++i)
        {
          if(itsBunch->Bin[i] < 0)
            {
              itsBunch->destroy(1,i);
              ++ne;
            }
        }
      reduce(ne,ne,OpAddAssign());

      size_t totalParticles_f = totalParticles_i - ne;

      if (hasWake)
        {
          IpplTimings::startTimer(WakeFieldTimer_m);
          if (!wf)
            INFOMSG("no wakefunction attached" << endl);
          wf->apply(*itsBunch);
          IpplTimings::stopTimer(WakeFieldTimer_m);

        }

      kickParticles(pusher);

      if (totalParticles_f > 0)
        {
          if (!bends) // none of the particles is in a bending element
            updateReferenceParticle();

          else /* at least one of the elements bends the beam; until all particles have left the bending elements we track the reference particle
                * as if it were a regular particle; from the moment when the reference particle has reached the bending field until it leaves
                * it again we rotate the bunch about the position of the reference particle such that the momentum of the reference particle points
                * in z direction
                */
            {
              recpgamma = 1./sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
              RefPartR_suv_m += RefPartP_suv_m * recpgamma * scaleFactor_m / 2.;

              // get field of reference particle
              bool RefPartinBend = getExternalField(RefPartR_suv_m, t, externalE, externalB);
              // FIXME: add self field

              RefPartR_suv_m /= vscaleFactor;

              kickReferenceParticle(externalE, externalB);

              if (RefPartinBend)
                {
                  updateSpaceOrientation();
                }

              // First update the momentum of the reference particle in zxy coordinate system, then update its position
              RefPartP_zxy_m(0) = SpaceOrientation_m[0] * RefPartP_suv_m(0) + SpaceOrientation_m[1] * RefPartP_suv_m(1) + SpaceOrientation_m[2] * RefPartP_suv_m(2);
              RefPartP_zxy_m(1) = SpaceOrientation_m[3] * RefPartP_suv_m(0) + SpaceOrientation_m[4] * RefPartP_suv_m(1) + SpaceOrientation_m[5] * RefPartP_suv_m(2);
              RefPartP_zxy_m(2) = SpaceOrientation_m[6] * RefPartP_suv_m(0) + SpaceOrientation_m[7] * RefPartP_suv_m(1) + SpaceOrientation_m[8] * RefPartP_suv_m(2);
              RefPartR_zxy_m += RefPartP_zxy_m * recpgamma * scaleFactor_m / 2.;

              RefPartP_suv_m = Vector_t(0.0, 0.0, sqrt(dot(RefPartP_suv_m,RefPartP_suv_m)));
              RefPartR_suv_m += RefPartP_suv_m * recpgamma / 2.;
              RefPartR_suv_m *= vscaleFactor;


            }
        }

      itsBunch->RefPart_R = RefPartR_zxy_m;
      itsBunch->RefPart_P = RefPartP_zxy_m;

      // calculate the dimensions of the bunch and add a small margin to them; then decide which elements have to be triggered
      // when an element is triggered memory is allocated and the field map is read in
      itsBunch->get_bounds(rmin,rmax);

      // trigger the elements
      margin = 3. * RefPartP_zxy_m(2) * recpgamma;
      margin = 0.01 > margin? 0.01: margin;
      for (fld_it = allElements.begin(); fld_it != fldend_it; ++fld_it)
        {
          if ((rmax(2) + margin) * scaleFactor_m > (*fld_it).Start  && !(*fld_it).goneLive)
            {
              (*fld_it).Element->goOnline();  // allocate memory and read in field map
              *gmsg << (*fld_it).Element->getName() << " gone live at step #" << step << " rmax = " << rmax(2) * scaleFactor_m << endl;
              (*fld_it).goneLive = true;
            }
          if ((rmin(2) - margin) * scaleFactor_m > (*fld_it).End && !(*fld_it).goneOff)
            {
              (*fld_it).Element->goOffline();  // free memory if no other element uses the same field map
              *gmsg << (*fld_it).Element->getName() << " gone off at step #" << step  << " rmin = " << rmin(2) * scaleFactor_m << endl;
              (*fld_it).goneOff = true;
            }
        }


      // start particle loop part 2
      for (int i = 0; i < itsBunch->getLocalNum(); ++i)
        {
          /** \f[ \vec{x}_{n+1} = \vec{x}_{n+1/2} + \frac{1}{2}\vec{v}_{n+1/2}\quad (= \vec{x}_{n+1/2} + \frac{\Delta t}{2} \frac{\vec{\beta}_{n+1/2}\gamma_{n+1/2}}{\gamma_{n+1/2}}) \f]
           * \code
           * R[i] += 0.5 * P[i] * recpgamma;
           * \endcode
           */
          pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);
          //and scale back to dimensions
          itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i]);
          // update local coordinate system
          itsBunch->X[i] /= vscaleFactor;
          pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i], allSections[itsBunch->LastSection[i]].Orientation), itsBunch->getdT());
          itsBunch->X[i] *= vscaleFactor;
          //reset time step if particle was emitted in the first half-step
          //the particle is now in sync with the simulation timestep
          itsBunch->dt[i] = itsBunch->getdT();
        }
      IpplTimings::stopTimer(timeIntegrationTimer2_m);

      if(totalParticles_f > minBinEmitted)
        itsBunch->boundp();

      if(totalParticles_i != totalParticles_f) {
        *gmsg << "======================================================================" << endl;
        *gmsg << "Particle loss at t= " << t << " n= " << totalParticles_i - totalParticles_f << " Ntot= " << totalParticles_f << endl;
        *gmsg << "======================================================================" << endl;
        totalParticles_i = totalParticles_f;
      }

      t += itsBunch->getdT(); //t after a full global timestep with dT "synchronization point" for simulation time

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
                writePhaseSpace(sposRef);
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
          if (EndOfLineReached)
            break;
        }
    }

  /**
     If use if ((*fld_it).goneLive) && !(*fld_it).goneOff)
     we prevents all cores to participate in putting the Screen
     offline. However all cores must participate in offline call i.e. 
     the write, if not the code hangs!

     OLD: code
     if ((*fld_it).goneLive) && !(*fld_it).goneOff)
      {
        (*fld_it).Element->goOffline();  // free memory if no other element uses the same field map
        *gmsg << (*fld_it).Element->getName() << " gone off at step #" << step  << " rmin = " << rmin(2) * scaleFactor_m << endl;
        (*fld_it).goneOff = true;
      }

  */


  for (fld_it = allElements.begin(); fld_it != allElements.end(); ++fld_it)
    {
      (*fld_it).Element->goOffline();  // free memory if no other element uses the same field map
      *gmsg << (*fld_it).Element->getName() << " gone off at step #" << step  << " rmin = " << rmin(2) * scaleFactor_m << endl;
      (*fld_it).goneOff = true;
    }
  
  *gmsg << "done executing ParallelTTracker" << endl;
}

// 2007/04/19 CKR
void ParallelTTracker::visitBeamline(const Beamline & bl){
  itsBeamline->iterate(*dynamic_cast<BeamlineVisitor*>(this),false);
}


static Vector implicitIntStep(const Vector &zin, const VSeries &f, const MxSeries gradf, double ds, int nx)
{
  Vector zf;
  return zf;
}


static Vector implicitInt2(const Vector &zin, const VSeries &f, double s, double ds, int nx, int cx)
{
  Vector zf=zin;
  return zf;
}

static Vector implicitInt4(const Vector &zin, const VSeries &f, double s, double ds, int nx, int cx)
{
  Vector zf=zin;
  return zf;
}


static Vector fixedPointInt2(const Vector &zin, const VSeries &f, double ds, int nx)
{
  Vector zf;
  return zf;
}

static Vector fixedPointInt4(const Vector &zin, const VSeries &f, double s, double ds, int nx, int cx)
{
  Vector zf=zin;
  return zf;
}
