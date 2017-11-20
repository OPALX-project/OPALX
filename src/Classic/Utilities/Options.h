#ifndef OPTIONS_HH
#define OPTIONS_HH

// ------------------------------------------------------------------------
// $RCSfile: Options.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Struct: Options
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:48 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------


// Namespace Options.
// ------------------------------------------------------------------------
/// The global OPAL option flags.
//  This namespace contains the global option flags.

#include "OptionTypes.h"
#include "Utilities/ClassicRandom.h"

namespace Options {
    /// Echo flag.
    //  If true, print an input echo.
    extern bool echo;

    /// Info flag.
    //  If true, print informative messages.
    extern bool info;

    extern bool csrDump;

    /// ppdebug flag.
    //  If true, use special initial velocity distribution for parallel plate and print special debug output .
    extern bool ppdebug;

    /// if true create symmetric distribution
    extern bool cZero;

    /// If true HDF5 files are written
    extern bool enableHDF5;

    extern bool asciidump;

    // If the distance of a particle to bunch mass larger than remotePartDel times of the rms size of the bunch in any dimension,
    // the particle will be deleted artifically to hold the accuracy of space charge calculation. The default setting of -1 stands for no deletion.
    extern double remotePartDel;

    extern double beamHaloBoundary;

    extern bool writeBendTrajectories;

    extern OPENMODE openMode;

    extern bool idealized;

    // CKR: nowhere used
    // // true if in bet mode
    // extern bool bet;

    /// Trace flag.
    //  If true, print CPU time before and after each command.
    extern bool mtrace;

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
    //  - GLOBAL, in Cartesian frame of the global particle
    //  - BUNCH_MEAN, in Cartesian frame of the bunch mean
    //  - REFERENCE, in Cartesian frame of the reference (0) particle
    extern DumpFrame psDumpFrame;

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

    /// Do closed orbit and tune calculation only.
    extern bool cloTuneOnly;

    extern int version;
    
#ifdef ENABLE_AMR
    extern bool amr;
#endif
    
    extern bool memoryDump;
}

#endif // OPAL_Options_HH