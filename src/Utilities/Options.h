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
  extern bool opal8;
  extern bool mad8;
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

  /// The frequency to dump the phase space i.e. timestep%psDumpFreq
  extern int psDumpFreq;

  /// The frequency to dump single particle trajectory of particles with ID = 0 & 1 
  extern int sptDumpFreq;

  /// The frequency to do particles repartition for better load balance between nodes 
  extern int repartFreq;

  /// Defines if an when we do the eake field calculation
  extern int wakeCalcStep;

  /// Defines the number of bins for the line density calculation
  extern int wakeBins;
}

#endif // OPAL_Options_HH
