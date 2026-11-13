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
#include "Utilities/Options.h"
#include "Utilities/ExpressionRandom.h"

#include <string>

namespace Options {
    // The global program options.
    bool echo = false;

    bool info           = true;
    int fft2d5BatchSize = 0;

    bool mtrace = false;

    bool warn = true;

    Random rangen;
    int seed = 123456789;

    int psDumpFreq = 10;

    int statDumpFreq = 10;

    int checkpointFreq = 0;

    int stepInfoFreq = 1;

    int printRankDistrFreq = 0;

    int repartFreq = 10;

    int minBinEmitted = 10;

    bool rhoDump = false;

    bool ebDump = false;

    bool rankDump = false;

    int autoPhase = 6;

    bool enableHDF5 = true;

    bool enableVTK = true;

    std::array<unsigned, 3> elementMappingGrid{8, 8, 8};
    double elementMappingReach = 0.5;

    bool asciidump = false;

    double boundpDestroy = 10;

    double beamHaloBoundary = 0;

    int version = 10000;

    double haloShift = 0.0;

    bool computePercentiles = false;

    bool useQMAttributes = false;  // corresponds to QM_MODE=SINGLE

    bool aggressiveStateSync = false;

    double loadBalancingThreshold = 0.05;  // Default threshold for triggering load balancing

    double loadBalancingPairWeight = 0.0;

}  // namespace Options
