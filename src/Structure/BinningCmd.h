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

class Inform;

/// The parameter that is used for binning.
enum class BinningParameter : short { VELOCITYZ = 0, POSITIONZ = 1, PZ = 2, GAMMAZ = 3 };

// The attributes of class BinningCmd.
namespace BINNING {
    enum {
        MAXBINS,         // Maximum number of bins for adaptive binning
        DESIREDWIDTH,    // Target / bias for bin width
        BINNINGALPHA,    // Aggressiveness of bin-number reduction
        BINNINGBETA,     // Aggressiveness of using wider bins
        PARAMETER,       // Which bunch attribute is used for binning
        DUMPBINSFILE,    // File name for dumping bins
        DUMPBINSFREQ,    // Frequency of dumping bins to a file (only used if DUMPBINSFILE is set)
        TABLEPRINTFREQ,  // Frequency of printing bin stats table to console (binned mode only)
        SIZE
    };
}

class BinningCmd : public Definition {
public:
    /// Exemplar constructor.
    BinningCmd();

    virtual ~BinningCmd();

    /// Make clone.
    virtual BinningCmd* clone(const std::string& name);

    /// Find named BINNING command.
    static BinningCmd* find(const std::string& name);

    int getMaxBins() const;
    double getDesiredWidth() const;
    double getBinningAlpha() const;
    double getBinningBeta() const;

    std::string getParameter();
    BinningParameter getParameterType() const;

    /// Get the file name for dumping bins to a file.
    std::string getDumpBinsFileName() const;

    /// Get the frequency of dumping bins to a file.
    int getDumpBinsFrequency() const;

    /// Get frequency (in global timesteps) of printing bin stats to console.
    /// Returns an integer >= 0. If 0, printing is disabled.
    int getTablePrintFrequency() const;

    /// Check if dumping bins to a file is enabled.
    bool dumpBinsToFile() const;

    /// Update the binning data (internal cache of attributes).
    virtual void update();

    /// Execute the BINNING command (currently a thin wrapper around update()).
    virtual void execute();

    Inform& printInfo(Inform& os) const;

private:
    // Not implemented.
    BinningCmd(const BinningCmd&);
    void operator=(const BinningCmd&);

    // Clone constructor.
    BinningCmd(const std::string& name, BinningCmd* parent);

    /// Recompute the enum-valued parameter from the string attribute.
    void setParameterType();

    std::string parameterName_m;
    BinningParameter parameterType_m;
};

inline Inform& operator<<(Inform& os, const BinningCmd& bc) { return bc.printInfo(os); }

#endif  // OPALX_BinningCmd_HH
