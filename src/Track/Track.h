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

#include "Algorithms/bet/EnvelopeBunch.h"
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

    Track(BeamSequence *, const PartData &, double dt, int maxtsteps, int stepsperturn, double zStop, int timeintegrator, int nslices);
    ~Track();

    /// The particle bunch to be tracked.
    PartBunch *bunch;

    EnvelopeBunch *slbunch;

    /// The reference data.
    PartData reference;

    /// The lattice to be tracked through.
    BeamSequence *use;

    /// The parser used during tracking.
    TrackParser parser;

    /// The block of track data.
    static Track *block;

    /// The initial timestep
    double dT;

    /// Maximal number of timesteps
    int maxTSteps;

    /// The timsteps per revolution period. ONLY available for OPAL-cycl.
    int stepsPerTurn;

    /// The location at which the simulation stops
    double zstop;

    /// The ID of time integrator
    // 0 --- RK-4(default)
    // 1 --- LF-2
    int timeIntegrator;

private:

    // Not implemented.
    Track();
    Track(const Track &);
    void operator=(const Track &);
};

#endif // OPAL_Track_HH
