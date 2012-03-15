#ifdef HAVE_ENVELOPE_SOLVER
// ------------------------------------------------------------------------
// $RCSfile: SLDataSink.hh,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SLDataSink
//
// ------------------------------------------------------------------------
//
// Revision History:
// $Date: 2003/01/23 13:29:44 $
// $Author: adelmann $
// $Log: SLDataSink.hh,v $

/////////////////////////////////////////////////////////////////////////////
// Class: SLDataSink
// Observer generates diagnostic output of the accelerator beam.
// The class computes statistical descriptors of particles' positions
// and momemta and writes those in files.
// It also writes particles' positions and momenta to files.
// prints initial and final particle positions into files
// accumulates beam statistics in Inform objects (files) during run
/////////////////////////////////////////////////////////////////////////////

#ifndef SLDataSink_H_
#define SLDataSink_H_

#include <fstream>
#include <vector>
#include <iostream>
#include "AbstractObjects/OpalData.h"
#include "Algorithms/bet/SLPartBunch.h"
#include "Ippl.h"

using namespace std;
/////////////////////////////////////////////////////////////////////////////

class SLDataSink
{
public:
  SLDataSink();    
  SLDataSink(int restartStep);    

  ~SLDataSink();

private:
 
  SLDataSink( const SLDataSink & ) { }
  SLDataSink & operator = ( const SLDataSink & ) { return *this; }

public:
 

  /** Write Stat Data 
  */
  void writeStatData(SLPartBunch &beam,
		     const string fname, 
		     const string el);

  void writePhaseSpace(SLPartBunch &beam);

private:
      
  bool firstWriteToStat_m;

  IpplTimings::TimerRef StatMarkerTimer_m;
  
  // to compensate IMPACT-T wrong spos
  double sshift_m;

};

#endif // SLDataSink_H_
#endif
/***************************************************************************
 * $RCSfile: SLDataSink.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 13:29:44 $
 ***************************************************************************/

