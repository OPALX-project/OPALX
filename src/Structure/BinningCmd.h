//
// Class BinningCmd
//   The class for the OPAL BINNING command.
//   A BINNING definition is used to configure adaptive
//   binning parameters that can later be attached to a
//   field solver.
//
// Copyright (c) 200x - 2026, Paul Scherrer Institut,
// Villigen PSI, Switzerland
//
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

#ifndef OPALX_BinningCmd_HH
#define OPALX_BinningCmd_HH

#include <string>

#include "AbstractObjects/Definition.h"
#include "Attributes/Attributes.h"
#include "Track/ObjectValues.h"

class Inform;

/// The parameter that is used for binning.
enum class BinningParameter : short { VELOCITYZ = 0, POSITIONZ = 1, PZ = 2, GAMMAZ = 3 };

// The attributes of class BinningCmd.
namespace BINNING {
    enum {
        MAXBINS,          // Maximum number of bins for adaptive binning
        DESIREDWIDTH,     // Target / bias for bin width
        BINNINGALPHA,     // Aggressiveness of bin-number reduction
        BINNINGBETA,      // Aggressiveness of using wider bins
        PARAMETER,        // Which bunch attribute is used for binning
        ADAPTIVEBINNING,  // Enable adaptive bin merging after uniform rebinning
        DUMPBINSFILE,     // File name for dumping bins
        DUMPBINSFREQ,     // Frequency of dumping bins to a file (only used if DUMPBINSFILE is set)
        TABLEPRINTFREQ,   // Frequency of printing bin stats table to console (binned mode only)
        SIZE
    };
}

class BinningCmd : public Definition {
public:
    /// Exemplar constructor.
    BinningCmd();

    virtual ~BinningCmd();

    /// Find named BINNING command.
    static BinningCmd* find(const std::string& name);

    static int getMaxBins(const ObjectValues& values);
    static double getDesiredWidth(const ObjectValues& values);
    static double getBinningAlpha(const ObjectValues& values);
    static double getBinningBeta(const ObjectValues& values);
    static bool getAdaptiveBinning(const ObjectValues& values);

    static std::string getParameter(const ObjectValues& values);

    /// The PARAMETER as an enum. @throws OpalException for an unknown PARAMETER.
    static BinningParameter getParameterType(const ObjectValues& values);

    /// Get the file name for dumping bins to a file; ".json" is added if it has no extension.
    static std::string getDumpBinsFileName(const ObjectValues& values);

    /// Get the frequency of dumping bins to a file.
    static int getDumpBinsFrequency(const ObjectValues& values);

    /// Get frequency (in global timesteps) of printing bin stats to console.
    /// Returns an integer >= 0. If 0, printing is disabled.
    static int getTablePrintFrequency(const ObjectValues& values);

    /// Check if dumping bins to a file is enabled.
    static bool dumpBinsToFile(const ObjectValues& values);

    /// Check PARAMETER, DUMPBINSFREQ and TABLEPRINTFREQ when the statement is read.
    virtual void execute();

    static Inform& printInfo(Inform& os, const ObjectValues& values);

private:
    // Not implemented.
    BinningCmd(const BinningCmd&);
    void operator=(const BinningCmd&);
};

#endif  // OPALX_BinningCmd_HH
