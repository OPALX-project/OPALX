#ifndef OPAL_OPTIONS_HH
#define OPAL_OPTIONS_HH

#include "Utilities/Random.h"

namespace Options {

    // CKR: nowhere used
    // // true if in bet mode
    // extern bool bet;

    /// Trace flag.
    //  If true, print CPU time before and after each command.
    extern bool mtrace;

    /// Verify flag.
    //  If true, print warning about undefined variables.
    extern bool verify;

    /// Warn flag.
    //  If true, print warning messages.
    extern bool warn;

    /// Random generator.
    //  The global random generator.
    extern Random rangen;

    /// The current random seed.
    extern int seed;

    /// The frequency to dump the phase space, i.e.dump data when step%psDumpFreq==0
    extern int psDumpFreq;

    // CKR: nowhere used
    // /// Dump centroid when R >rDump
    // extern double rDump;

    /// The frequency to dump statistical values, e.e. dump data when step%statDumpFreq==0
    extern int statDumpFreq;

    /// The frequency to dump single particle trajectory of particles with ID = 0 & 1
    extern int sptDumpFreq;

    /// The frequency to do particles repartition for better load balance between nodes
    extern int repartFreq;

    /// The frequency to reset energy bin ID for all particles
    extern int rebinFreq;

    /// phase space dump flag for OPAL-cycl
    //  if true, dump phase space after each turn
    extern bool psDumpEachTurn;

    // Governs how often boundp_destroy is called to destroy lost particles
    // Mainly used in the CyclotronTracker as of now -DW
    extern int boundpDestroyFreq;

    /// flag to decide in which coordinate frame the phase space will be dumped for OPAL-cycl
    // if true, in local Cartesian frame, otherwise in global Cartesian frame
    extern bool psDumpLocalFrame;

    /// The frequency to solve space charge fields.
    extern int scSolveFreq;

    // How many small timesteps are inside the large timestep used in multiple time stepping (MTS) integrator
    extern int mtsSubsteps;

    /// this allows to repeat tracks starting always at the begining of the lattice and
    /// generates a new distribution

    extern bool scan;

    extern bool rhoDump;

    extern bool ebDump;

    // CKR: nowhere used
    // extern bool efDump;

    // if true opal find the phases in the cavities, such that the energy gain is at maximum
    extern int autoPhase;

    /// The frequency to dump the particle-geometry surface interation data.
    extern int surfDumpFreq;

    /// RCG: cycle length
    extern int numBlocks;

    /// RCG: number of recycle blocks
    extern int recycleBlocks;

    /// number of old left hand sides used to extrapolate a new start vector
    extern int nLHS;

    extern std::string rngtype;

    /// if true
    extern bool schottkyCorrection;

    ///
    extern double schottkyRennormalization;

    // CKR: nowhere used
    // /// If true change the time step during emision from the cathode so that that one bin
    // /// of the time histogram that describes the longitudinal beam distribution is emitted
    // /// during each time step. If false then the time step during emission is set so that one
    // /// energy bin of the beam is emitted during each time step.
    // extern bool fineEmission;
}

#endif
