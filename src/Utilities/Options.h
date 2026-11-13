//
// Namespace Options
//   The global OPAL option flags.
//   This namespace contains the global option flags.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPTIONS_HH
#define OPTIONS_HH

#include "Utilities/ExpressionRandom.h"

#include <array>
#include <string>

namespace Options {
    /// Echo flag.
    //  If true, print an input echo.
    extern bool echo;

    /// Info flag.
    //  If true, print informative messages.
    extern bool info;

    /// Requested FFT2D5 batch size; zero uses serial slices on one rank and batch one with MPI.
    extern int fft2d5BatchSize;

    /// Trace flag.
    //  If true, print CPU time before and after each command.
    extern bool mtrace;

    /// Warn flag.
    //  If true, print warning messages.
    extern bool warn;

    //  The global random generator.
    extern Random rangen;

    /// The current random seed.
    extern int seed;

    /// The frequency to dump the phase space, i.e. dump data when
    /// step%psDumpFreq==0; 0 disables phase space dumps.
    extern int psDumpFreq;

    /// The frequency to dump statistical values, i.e. dump data when
    /// step%statDumpFreq==0; 0 disables statistical dumps.
    extern int statDumpFreq;

    /// The frequency to write a restart checkpoint after completed integration steps;
    /// 0 disables checkpointing.
    extern int checkpointFreq;

    /// The frequency to print per-step tracking status lines; 0 disables them.
    extern int stepInfoFreq;

    /// The frequency to print per-rank particle distribution tables; 0 disables them.
    extern int printRankDistrFreq;

    /// If true, store `Q`/`M` as per-particle attributes (accessible via
    /// `getQView()`/`getMView()`). If false (default), use a single shared value per container to
    /// save memory.
    extern bool useQMAttributes;

    /// The frequency to do particles repartition for better load balance between nodes
    extern int repartFreq;

    /// The number of bins that have to be emitted before the bin are squashed into a single bin
    extern int minBinEmitted;

    extern bool rhoDump;

    extern bool ebDump;

    extern bool rankDump;

    // the number of refinements of the search range for the phase with maximum energy
    // if eq 0 then no autophase
    // if true opal find the phases in the cavities, such that the energy gain is at maximum
    extern int autoPhase;

    /// If true HDF5 files are written
    extern bool enableHDF5;

    /// If true VTK files are written
    extern bool enableVTK;

    /// Cells of the lab-frame element grid along x, y and z (default 8 each); see
    /// GeneralElementMapping.
    extern std::array<unsigned, 3> elementMappingGrid;
    /// Default minimum transverse half-width [m] of an element box (default 0.5); an element's
    /// REACH attribute overrides it.
    extern double elementMappingReach;

    extern bool asciidump;

    /// Governs how many sigmas away particles are deleted
    extern double boundpDestroy;

    extern double beamHaloBoundary;

    /// opal version of input file
    extern int version;

    /// The constant parameter C to shift halo, by < w^4 > / < w^2 > ^2 - C (w=x,y,z)
    extern double haloShift;

    extern bool computePercentiles;

    /// If true, every `BunchStateHandler` setter performs an `ippl::Comm->allreduce`
    /// so that the flag is guaranteed to be consistent across MPI ranks, even if the
    /// caller set it on a subset of ranks. Off by default because the extra
    /// collective on every state change has a noticeable cost in tight loops; enable
    /// it only for debugging or in contexts where rank-local divergence is possible.
    extern bool aggressiveStateSync;

    /// The threshold for triggering load balancing. If the ratio difference of particles in a rank
    /// exceeds this threshold, load balancing will be triggered. Default is 0.05 (5%). This
    /// threshold is only tested every `repartFreq` steps.
    extern double loadBalancingThreshold;

    /// Share of the ORB weight that follows the P3M pair count instead of the particle count, in
    /// [0, 1]. Default is 0 (particle counts only).
    extern double loadBalancingPairWeight;
}  // namespace Options

#endif  // OPAL_Options_HH
