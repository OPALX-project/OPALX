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


// Class Track
// ------------------------------------------------------------------------

Track *Track::block = 0;


Track::Track(BeamSequence *u, const PartData &ref, double dt, int maxtsteps):
  bunch(), 
  reference(ref), 
  use(u), 
  dT(dt),
  maxTSteps(maxtsteps),
  parser()
{}


Track::~Track()
{}
