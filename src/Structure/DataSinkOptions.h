//
// Struct DataSinkOptions
//   The OPTION values DataSink and its writers read.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPAL_DATA_SINK_OPTIONS_H
#define OPAL_DATA_SINK_OPTIONS_H

/**
 * @brief The OPTION values of one RUN that DataSink and its writers read, fixed when the sink
 *        and the phase-space files are made.
 */
struct DataSinkOptions {
    bool enableHDF5 = false;  ///< Write the phase space to H5 files.
    int psDumpFreq  = 0;      ///< Stored as the "dump frequency" attribute of the H5 files.
    bool rankDump   = false;  ///< H5 dataset "rank": the MPI rank of each particle.
    bool ebDump     = false;  ///< H5 datasets of the E and B fields at the particles.
    bool rhoDump    = false;  ///< H5 unit attribute of the charge density.
    /// The stat column partsOutside counts particles beyond this many rms beam sizes; 0 skips it.
    double beamHaloBoundary = 0;
};

#endif  // OPAL_DATA_SINK_OPTIONS_H
