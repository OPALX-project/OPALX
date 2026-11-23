//
// Class StatWriter
//   This class writes bunch statistics (*.stat).
//
// Copyright (c) 2019, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
//                     Christof Metzger-Kraus, Open Sourcerer
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
#ifndef OPAL_STAT_WRITER_H
#define OPAL_STAT_WRITER_H

#include "StatBaseWriter.h"

#include "Moments/BeamMoments.h"

class StatWriter : public StatBaseWriter {
public:
    /// @param inputFileName Input file named in the header.
    StatWriter(const std::string& fname, const std::string& inputFileName, bool restart);

    /** \brief Write statistical data.
     *
     * Writes statistical beam data to proper output file. This is information such as RMS beam
     * parameters etc.
     *
     * \param beam The beam.
     * \param FDext The external fields at the reference particle: FDext[0] = B (T),
     * FDext[1] = E (V/m), written in MV/m.
     */
    void write(
            PartBunch_t& beam, const opalx::moments::BeamMoments& moments,
            Vector_t<double, 3> FDext[], int reportedBinCount, const size_t npOutside = 0,
            size_t particleContainerIndex = 0);

    /// Moment groups write() reads for a non-empty container; Polarization only with spin.
    static constexpr opalx::moments::Request requiredMoments(bool hasSpin) {
        return hasSpin ? opalx::moments::allMoments
                       : opalx::moments::allMoments.without(opalx::moments::Moment::Polarization);
    }

private:
    void fillHeader(const std::string& species = "", bool hasSpin = false);

    std::string inputFileName_m;
};

#endif
