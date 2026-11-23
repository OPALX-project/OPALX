//
// Class DataSink
//   This class acts as an observer during the calculation. It generates diagnostic
//   output of the accelerated beam such as statistical beam descriptors of particle
//   positions, momenta, beam phase space (emittance) etc. These are written to file
//   at periodic time steps during the calculation.
//
//   This class also writes the full beam phase space to an H5 file at periodic time
//   steps in the calculation (this period is different from that of the statistical
//   numbers).
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Structure/DataSink.h"

#include "Structure/H5PartWrapper.h"
#include "Structure/H5PartWrapperForPT.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <algorithm>

std::string DataSink::diagnosticStemForContainer(
        const std::string& inputBasename, size_t numContainers, size_t index) {
    if (numContainers <= 1) {
        return inputBasename;
    }
    return inputBasename + "_c" + std::to_string(index);
}

DataSink::DataSink(
        const DataSinkOptions& options, const std::vector<H5PartWrapper*>& h5wrappers,
        size_t numParticleContainers, const std::string& outputBasename,
        const std::string& inputFileName, bool appendStatistics)
    : options_m(options), StatMarkerTimer_m(IpplTimings::getTimer("Write Stat")) {
    for (size_t i = 0; i < numParticleContainers; ++i) {
        const std::string stem =
                diagnosticStemForContainer(outputBasename, numParticleContainers, i);
        statWriters_m.push_back(
                std::make_unique<StatWriter>(stem + ".stat", inputFileName, appendStatistics));
    }

    if (options_m.enableHDF5) {
        if (h5wrappers.size() != numParticleContainers) {
            throw OpalException(
                    "DataSink::DataSink",
                    "HDF5 enabled: expected " + std::to_string(numParticleContainers)
                            + " H5PartWrapper(s), got " + std::to_string(h5wrappers.size()) + ".");
        }
        for (H5PartWrapper* w : h5wrappers) {
            h5Writers_m.push_back(std::make_unique<H5Writer>(w));
        }
    } else if (!h5wrappers.empty()) {
        throw OpalException("DataSink::DataSink", "H5 wrappers passed while HDF5 is disabled.");
    }
}

opalx::moments::Request DataSink::requiredMoments(bool psDump, bool statDump, bool hasSpin) {
    opalx::moments::Request request;
    if (statDump) {
        request |= StatWriter::requiredMoments(hasSpin);
    }
    if (psDump && Options::enableHDF5) {
        request |= H5PartWrapperForPT::requiredMoments;
    }
    return request;
}

opalx::moments::Request DataSink::momentRequest(bool psDump, bool statDump, bool hasSpin) const {
    opalx::moments::Request request;
    if (statDump) {
        request |= StatWriter::requiredMoments(hasSpin);
    }
    if (psDump && options_m.enableHDF5) {
        request |= H5PartWrapperForPT::requiredMoments;
    }
    return request;
}

void DataSink::dumpH5(
        PartBunch_t& beam, const std::vector<std::array<Vector_t<double, 3>, 2>>& fdextPerContainer,
        std::span<const opalx::moments::BeamMoments> moments) const {
    if (!options_m.enableHDF5) {
        return;
    }

    const size_t n = beam.getNumParticleContainers();
    if (fdextPerContainer.size() < n || moments.size() < n) {
        throw OpalException(
                "DataSink::dumpH5", "fdextPerContainer size ("
                                            + std::to_string(fdextPerContainer.size())
                                            + ") < num containers (" + std::to_string(n) + ").");
    }

    for (size_t i = 0; i < n; ++i) {
        if (i >= h5Writers_m.size()) {
            break;
        }
        auto pc = beam.getParticleContainer(i);
        if (!pc || pc->getTotalNum() == 0) {
            continue;
        }
        Vector_t<double, 3> fd[2] = {fdextPerContainer[i][0], fdextPerContainer[i][1]};
        h5Writers_m[i]->writePhaseSpace(beam, moments[i], fd, i);
    }
}

void DataSink::dumpSDDS(
        PartBunch_t& beam, const std::vector<std::array<Vector_t<double, 3>, 2>>& fdextPerContainer,
        std::span<const opalx::moments::BeamMoments> moments, int reportedBinCount) const {
    const size_t n = beam.getNumParticleContainers();
    if (fdextPerContainer.size() < n || moments.size() < n) {
        throw OpalException(
                "DataSink::dumpSDDS", "fdextPerContainer size ("
                                              + std::to_string(fdextPerContainer.size())
                                              + ") < num containers (" + std::to_string(n) + ").");
    }

    IpplTimings::startTimer(StatMarkerTimer_m);

    for (size_t i = 0; i < n; ++i) {
        if (i >= statWriters_m.size()) {
            break;
        }
        auto pc = beam.getParticleContainer(i);
        if (!pc) {
            continue;
        }

        size_t npOutside = 0;
        if (pc->getTotalNum() > 0 && options_m.beamHaloBoundary > 0) {
            npOutside = beam.calcNumPartsOutside(options_m.beamHaloBoundary * moments[i].rmsR());
        }

        Vector_t<double, 3> fd[2] = {fdextPerContainer[i][0], fdextPerContainer[i][1]};
        statWriters_m[i]->write(beam, moments[i], fd, reportedBinCount, npOutside, i);
    }

    IpplTimings::stopTimer(StatMarkerTimer_m);
}

void DataSink::dumpBinConfig(
        long long step, double time, bool preMerge, const std::vector<std::size_t>& binCounts,
        const std::vector<double>& binWidths, double xMin, const std::string& fileName) {
    if (ippl::Comm->rank() != 0) {
        return;
    }

    Inform m("DataSink::dumpBinConfig");

    if (binCounts.size() != binWidths.size()) {
        m << level4 << "Invalid bin configuration: binCounts.size() = " << binCounts.size()
          << ", binWidths.size() = " << binWidths.size() << endl;
        throw OpalException(
                "DataSink::dumpBinConfig", "binCounts and binWidths must have the same length.");
    }

    if (!binConfigWriter_m) {
        m << level4 << "Creating BinConfigWriter for JSON binning output file \"" << fileName
          << "\"." << endl;
        binConfigWriter_m = std::make_unique<BinConfigWriter>(fileName);
    }

    m << level5 << "Dumping bin configuration snapshot: step=" << step << ", time=" << time
      << ", preMerge=" << (preMerge ? 1 : 0) << ", nBins=" << binCounts.size() << ", xMin=" << xMin
      << endl;

    binConfigWriter_m->writeEntry(step, time, preMerge, binCounts, binWidths, xMin);
}

void DataSink::rewindToCheckpoint(PartBunch_t& beam) {
    for (const auto& writer : h5Writers_m) {
        writer->rewindToTime(beam.getT(), beam.getdT());
    }

    const std::size_t count = std::min(statWriters_m.size(), beam.getNumParticleContainers());
    for (std::size_t i = 0; i < count; ++i) {
        const auto pc = beam.getParticleContainer(i);
        if (!pc || !statWriters_m[i]->exists()) {
            continue;
        }
        statWriters_m[i]->rewindToSpos(pc->get_sPos());
    }
}
