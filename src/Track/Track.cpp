// ------------------------------------------------------------------------
// $RCSfile: Track.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Struct: Track
//   This structure holds all data for tracking.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:46 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Track/Track.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"
// Class Track
// ------------------------------------------------------------------------

Track *Track::block = 0;


/**

 Track is asking the dictionary if already a 
 particle bunch was allocated. If that is the
 case Track is using the allready allocted bunch,
 otherwise a new bunch ich allocated in the dictionary.
*/


Track::Track(BeamSequence *u, const PartData &ref, double dt, int maxtsteps, int stepsperturn):
  reference(ref), 
  use(u), 
  dT(dt),
  maxTSteps(maxtsteps),
  stepsPerTurn(stepsperturn),
  parser()
{
  if (Options::bet) {
#ifdef HAVE_ENVELOPE_SOLVER
    if (!OPAL.hasSLBunchAllocated()) 
      OPAL.setSLPartBunch(new SLPartBunch(&ref));
    
    slbunch = OPAL.getSLPartBunch();
    INFOMSG("* ********************************************************************************** "<< endl);
    INFOMSG("  Selected Tracking Method == PARALLEL-SLICE, slbunch initialized                   " << endl);
    INFOMSG("* ********************************************************************************** " << endl);
#else
    ERRORMSG("* ********************************************************************************** "<< endl);
    ERRORMSG("  Selected Tracking Method == PARALLEL-SLICE, but the code is not compiled for the " << endl);
    ERRORMSG("  Use --envelope-solver in the configure make process   " << endl);
    ERRORMSG("* ********************************************************************************** " << endl);
    exit(1);
#endif
  }
  else {
    if (!OPAL.hasBunchAllocated()) 
      OPAL.setPartBunch(new PartBunch(&ref));
    
    bunch = OPAL.getPartBunch();
  }
}


Track::~Track()
{}
