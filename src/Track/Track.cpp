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

// Class Track
// ------------------------------------------------------------------------

Track *Track::block = 0;


/**

 Track is asking the dictionary if already a 
 particle bunch was allocated. If that is the
 case Track is using the allready allocted bunch,
 otherwise a new bunch ich allocated in the dictionary.
*/


Track::Track(BeamSequence *u, const PartData &ref, double dt, int maxtsteps):
  reference(ref), 
  use(u), 
  dT(dt),
  maxTSteps(maxtsteps),
  parser()
{
  if (!OPAL.hasBunchAllocated()) 
    OPAL.setPartBunch(new PartBunch(&ref));
  
  bunch = OPAL.getPartBunch();
}


Track::~Track()
{}
