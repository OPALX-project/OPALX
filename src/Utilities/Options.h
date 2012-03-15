#ifndef OPAL_Options_HH
#define OPAL_Options_HH

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

#include "Utilities/Random.h"


// Namespace Options.
// ------------------------------------------------------------------------
/// The global OPAL option flags.
//  This namespace contains the global option flags.

namespace Options {

  /// Echo flag.
  //  If true, print an input echo.
  extern bool echo;

  /// Info flag.
  //  If true, print informative messages.
  extern bool info;

  /// MAD-8 flag.
  //  If true, give output in MAD-8 format.
  extern bool mad8;
  extern bool opal8;

  // true if in bet mode
  extern bool bet;


  /// Trace flag.
  //  If true, print CPU time before and after each command.
  extern bool mtrace;

  /// Verify flag.
  //  If true, print warning about undefined variables.
  extern bool verify;

  /// Warn flag.
  //  If true, print warning messages.
  extern bool warn;

  /// The table format flag.
  //  If true, print tables in TFS format.
  extern bool tfsFormat;

  /// Random generator.
  //  The global random generator.
  extern Random rangen;

  /// The current random seed.
  extern int seed;

  /// If we do a phase scan
  extern bool doPhaseScan;

  /// The frequency to dump the phase space, i.e.dump data when step%psDumpFreq==0
  extern int psDumpFreq;

  /// Dump centroid when R >rDump
  extern double rDump;

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

  /// flag to decide in which coordinate frame the phase space will be dumped for OPAL-cycl
  // if true, in local Cartesian frame, otherwise in global Cartesian frame
  extern bool psDumpLocalFrame;


  /// The frequency to solve space charge fields.
  extern int scSolveFreq;
  
  /// this allows to repeat tracks starting always at the begining of the lattice and 
  /// generates a new distribution
 
  extern bool scan;

  extern bool rhoDump;

  extern bool efDunp;

  // if true opal find the phases in the cavities, such that the energy gain is at maximum 
  extern bool autoPhase;	
}

#endif // OPAL_Options_HH
