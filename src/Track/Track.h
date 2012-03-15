#ifndef OPAL_Track_HH
#define OPAL_Track_HH

// ------------------------------------------------------------------------
// $RCSfile: Track.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Struct: Track
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:46 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/PartBunch.h"
#include "Algorithms/PartData.h"
#include "Algorithms/Particle.h"
#include "Track/TrackParser.h"

class BeamSequence;
class TrackParser;


// Class Track
// ------------------------------------------------------------------------
/// Hold data for tracking.
//  Acts as a communication area between the various tracking commands.

class Track {

public:

  Track(BeamSequence *, const PartData &, double dt, int maxtsteps);
  ~Track();

  /// The particle bunch to be tracked.
  PartBunch *bunch;

  /// The reference data.
  PartData reference;

  /// The lattice to be tracked through.
  BeamSequence *use;

  /// The parser used during tracking.
  TrackParser parser;

  /// The block of track data.
  static Track *block;

  /// The intial timestep
  double dT;
  
  /// Maximal number of timesteps 
  int maxTSteps;
  
private:

  // Not implemented.
  Track();
  Track(const Track &);
  void operator=(const Track &);
};

#endif // OPAL_Track_HH
