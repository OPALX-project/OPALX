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
    extern int remotePartDel;
}

#endif // OPAL_Options_HH
