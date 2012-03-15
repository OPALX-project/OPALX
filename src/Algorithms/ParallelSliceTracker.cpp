#ifdef HAVE_ENVELOPE_SOLVER
// ------------------------------------------------------------------------
// $RCSfile: ParallelSliceTracker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelSliceTracker
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
#include "Algorithms/ParallelSliceTracker.h"

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

// Class ParallelSliceTracker
// ------------------------------------------------------------------------

ParallelSliceTracker::ParallelSliceTracker(const Beamline &beamline,
                                   const PartData &reference,
		                           bool revBeam, 
                                   bool revTrack):
  Tracker(beamline, reference, revBeam, revTrack),
  allElements(),
  allSections()
{
  itsBeamline = dynamic_cast<Beamline*>(beamline.clone());
}


ParallelSliceTracker::ParallelSliceTracker(const Beamline &beamline,
                                   SLPartBunch &bunch,
                                   SLDataSink &ds,
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


ParallelSliceTracker::~ParallelSliceTracker()
{
  for (FieldListIterator compindex = allElements.begin(); compindex != allElements.end(); ++compindex)
    {
      delete (*compindex).Element;
    }
  allElements.clear();
  allSections.clear();
  delete itsBeamline;
}

void ParallelSliceTracker::visitAlignWrapper(const AlignWrapper &wrap)
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

void ParallelSliceTracker::visitBeamBeam(const BeamBeam &)
{
  *gmsg << "BeamBeam not implemented yet!" << endl;
  //  *gmsg << "In BeamBeam; "<< endl;
}


void ParallelSliceTracker::visitCollimator(const Collimator &coll)
{
  *gmsg << "Collimator not implemented yet!" << endl;
  //  *gmsg << "In Collimator; L= " << coll.getElementLength() << endl;
//   allElements.push_back(FieldListEntry(dynamic_cast<Collimator*>(coll.clone()),0.,0.));
}


void ParallelSliceTracker::visitCorrector(const Corrector &corr)
{
  *gmsg << "Corrector not implemented yet!" << endl;
  //    *gmsg << "In Corrector; L= " << corr.getElementLength() << endl;
//   allElements.push_back(FieldListEntry(dynamic_cast<Corrector*>(corr.clone()),0.,0.));
}


void ParallelSliceTracker::visitDiagnostic(const Diagnostic &diag)
{
  //     *gmsg << "In Diagnostic; L= " << diag.getElementLength() << endl;
//   allElements.push_back(FieldListEntry(dynamic_cast<Diagnostic*>(diag.clone()),0.,0.));
}


void ParallelSliceTracker::visitDrift(const Drift &drift)
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


void ParallelSliceTracker::visitLambertson(const Lambertson &lamb)
{
  *gmsg << "Lambertson not implemented yet!" << endl;
  //     *gmsg << "In Lambertson; L= " << lamb.getElementLength() << endl;
//   allElements.push_back(FieldListEntry(dynamic_cast<Lambertson*>(lamb.clone()),0.,0.));
}


void ParallelSliceTracker::visitMarker(const Marker &marker)
{
  //     *gmsg << "In Marker; L= " << marker.getElementLength() << endl;
//   allElements.push_back(FieldListEntry(dynamic_cast<Marker*>(marker.clone()),0.,0.));
}


void ParallelSliceTracker::visitMonitor(const Monitor &mon)
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


void ParallelSliceTracker::visitMultipole(const Multipole &mult)
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


void ParallelSliceTracker::visitRBend(const RBend &bend)
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
  allElements.push_back(FieldListEntry(elptr,startField, endField));

}


void ParallelSliceTracker::visitRFCavity(const RFCavity &as)
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

void ParallelSliceTracker::visitTravelingWave(const TravelingWave &as)
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


void ParallelSliceTracker::visitRFQuadrupole(const RFQuadrupole &rfq)
{
  *gmsg << "RFQuadrupole not implemented yet!" << endl;
  //     *gmsg << "In RFQuadrupole; L= " << rfq.getElementLength() << endl;
//   allElements.push_back(FieldListEntry(dynamic_cast<RFQuadrupole*>(rfq.clone()),0.,0.));
}

void ParallelSliceTracker::visitSBend(const SBend &bend)
{
  *gmsg << "SBend not implemented yet!" << endl;
  //     *gmsg << "In SBend; L= " << bend.getElementLength() << endl;
//   allElements.push_back(FieldListEntry(dynamic_cast<SBend*>(bend.clone()),0.,0.));
}


void ParallelSliceTracker::visitSeparator(const Separator &sep)
{
  *gmsg << "Separator not implemented yet!" << endl;
  //  *gmsg << "In Separator; L= " << sep.getElementLength() << endl;
//   allElements.push_back(FieldListEntry(dynamic_cast<Separator*>(sep.clone()),0.,0.));

}


void ParallelSliceTracker::visitSeptum(const Septum &sept)
{
  *gmsg << "Septum not implemented yet!" << endl;
  //  *gmsg << "In Septum; L= " << sept.getElementLength() << endl;
//   allElements.push_back(FieldListEntry(dynamic_cast<Septum*>(sept.clone()),0.,0.));
}


void ParallelSliceTracker::visitSolenoid(const Solenoid &solenoid)
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

void ParallelSliceTracker::applyEntranceFringe(double angle, double curve,
                                           const BMultipoleField &field, double scale)
{
}


void ParallelSliceTracker::applyExitFringe(double angle, double curve,
                                       const BMultipoleField &field, double scale)
{
}

void ParallelSliceTracker::buildupFieldList()
{
  /** Build up a list of sections of fields: 
   *  for each section we get a list of pointers to the elements which contribute to
   *  the electromagnetic field in this section. Sections start and end where the fields
   *  of the elements start and end respectively. The field of an element can be conained 
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
      if (fabs((*fld_it).End - (*fld_it).Start) < tolerance)
        {
          allElements.erase(fld_it);
        }
      else
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
      *gmsg << "--- " << (*sec_it).Start << " m -- " << (*sec_it).End << " m ---------------------------\n";

      for (ElementIterator el_it = (*sec_it).Elements.begin(); 
           el_it != (*sec_it).Elements.end(); 
           ++el_it)
        *gmsg << (*el_it)->getName() << '\n';

    }
  *gmsg << "--- END   FIELD LIST ----------------------------------------------------------------------\n" << endl;


}

// 2007/04/19 CKR
void ParallelSliceTracker::execute()
{

  double recpgamma, gamma;
  double t = itsBunch->getT(); 
  double dt = itsBunch->getdT();
  double tEmission = itsBunch->getTEmission();

  long long step =  OPAL.inRestartRun() ? OPAL.getRestartStep() + 1 :lround(t/dt);

  Vector_t um, a, s;
  Vector_t externalE, externalB;
  Vector_t rmin, rmax;
  const Vector_t vscaleFactor = Vector_t(scaleFactor_m);

  bool EndOfLineReached;
  bool partOutOfBounds;
  bool bends;               // flag which indicates wheter any particle is within the influence of bending element.
                            // if this is the case we track the reference particle as if it were a real particle, 
                            // otherwise the reference particle is defined as the centroid particle of the bunch

  bool hasWake = false;     // flag which indicates wheter any particle is within the influence of a wake field

  FieldListIterator fld_it;
  FieldListIterator fldend_it;

  *gmsg << "executing ParallelSliceTracker, initial DT " << itsBunch->getdT() 
        << " [s]; max integration steps " << maxSteps_m << " step= " << step << endl
        << "the mass is: " << itsReference.getM()*1e-6 << " MeV, its charge: " << itsReference.getQ() << endl;
  
  itsBeamline->accept(*this);  // fill list allElements
  buildupFieldList();          // fill list allSections

  fldend_it = allElements.end();

  if (OPAL.inRestartRun()) {
    for (int i = 0; i < itsBunch->getLocalNumSlices(); ++i) {
      itsBunch->LastSection[i] = 0;
      for (int l = 0; l < allSections.size(); ++l) {
	if (itsBunch->Z[i](2) >= allSections[l].Start 
	    && 
	    itsBunch->Z[i](2) <= allSections[l].End) {
	  itsBunch->LastSection[i] = l;
	  break;
	}
      }
    }      
  }
  
  itsBunch->calcBeamParameters();

  for(step; step < maxSteps_m; ++step) {
    EndOfLineReached = true;  //check if any particle hasn't reached the end of the field from the last element
    bends = false;
    hasWake = false;
    
    double margin = 1e-6;

    itsBunch->calcBeamParameters();
    
    for (int i = 0; i < itsBunch->getLocalNumSlices(); ++i) {
      //reset helper vectors in each step since fields are added in getFieldstrenght
      externalB = Vector_t(0.0);
      externalE = Vector_t(0.0);

      if ( itsBunch->Z[i](2) < allSections.back().End ) {
	EndOfLineReached = false;
	int lastSection = itsBunch->LastSection[i];
	
	/// go ev. to the next section
	while (itsBunch->Z[i](2) >= allSections[lastSection].End)
	  lastSection = ++itsBunch->LastSection[i];
	    
	if (itsBunch->Z[i](2) >= allSections[lastSection].Start) {
	  bends = bends     || allSections[lastSection].bends; // check if any element bends the beam; make sure that this is communicated to all processors!
	  hasWake = hasWake || allSections[lastSection].hasWake;
	  
	  ElementIterator last_el_it = allSections[lastSection].Elements.end();
	  for (ElementIterator el_it = allSections[lastSection].Elements.begin(); el_it != last_el_it; ++el_it) {
	    partOutOfBounds = partOutOfBounds || (*el_it)->apply(i, itsBunch->get_rmean(), t+itsBunch->dt[i]/2.0, externalE, externalB);
	  }
	}
      }
	  
      /**
	 Field evaluation on element done
      */
	  
      for (fld_it = allElements.begin(); fld_it != fldend_it; ++fld_it)
        {
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
    }
    
    t += itsBunch->getdT(); //t after a full global timestep with dT "synchronization point" for simulation time
    
    itsBunch->setT(t);

    double sposRef = 1.;
    if (step % Options::psDumpFreq == 0 ) 
      writePhaseSpace(sposRef);           
  }    
  for (fld_it = allElements.begin(); fld_it != fldend_it; ++fld_it) {
    if ((*fld_it).goneLive && !(*fld_it).goneOff) {
      (*fld_it).Element->goOffline();  // free memory if no other element uses the same field map
      cout << (*fld_it).Element->getName() << " gone off at step #" << step  << " rmin = " << rmin(2) * scaleFactor_m << endl;
      (*fld_it).goneOff = true;
    }
  }
  *gmsg << "done executing ParallelSliceTracker" << endl;
}

// 2007/04/19 CKR
void ParallelSliceTracker::visitBeamline(const Beamline & bl){
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
#endif
