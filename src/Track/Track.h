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

#include "Algorithms/PartData.h"
#include "Track/TrackParser.h"

class BeamSequence;
class TrackParser;
class PartBunch;
class EnvelopeBunch;

// Class Track
// ------------------------------------------------------------------------
/// Hold data for tracking.
//  Acts as a communication area between the various tracking commands.

class Track {

public:

    Track(BeamSequence *, const PartData &, const std::vector<double> & dt,
          const std::vector<unsigned long long> & maxtsteps, int stepsperturn,
          const std::vector<double> & zStop, int timeintegrator, int nslices,
          double t0, double dtScInit, double deltaTau);
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
    std::vector<double> dT;

    // For AMTS integrator in OPAL-T
    double dtScInit, deltaTau;

    /// The ellapsed time of the beam can be used to propper
    /// start the beam when created in a cavity i.e. without emission
    double t0_m;

    /// Maximal number of timesteps
    std::vector<unsigned long long> localTimeSteps;

    /// The timsteps per revolution period. ONLY available for OPAL-cycl.
    int stepsPerTurn;

    /// The location at which the simulation stops
    std::vector<double> zstop;

    /// The ID of time integrator
    // 0 --- RK-4(default)
    // 1 --- LF-2
    // 2 --- MTS
    // 3 --- AMTS
    int timeIntegrator;

private:

    // Not implemented.
    Track();
    Track(const Track &);
    void operator=(const Track &);
};

#endif // OPAL_Track_HH