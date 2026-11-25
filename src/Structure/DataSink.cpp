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
#include "Utilities/Timer.h"
#include "Utilities/Util.h"

#include <algorithm>
#include <cmath>
#include <iomanip>

extern Inform* gmsg;

namespace {
    /// Whether an output with frequency @p freq (0: never) is due at global step @p step: the
    /// last step of every cycle of @p freq steps.
    bool outputDue(int freq, long long step) { return freq > 0 && (step % freq) + 1 == freq; }
}  // namespace

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

DataSink::StepOutputs DataSink::outputsDue(long long step) const {
    StepOutputs due;
    due.phaseSpace = outputDue(options_m.psDumpFreq, step);
    due.statistics = outputDue(options_m.statDumpFreq, step);
    due.stepLine   = options_m.stepInfoFreq > 0 && step % options_m.stepInfoFreq == 0;
    due.rankLoad   = outputDue(options_m.printRankDistrFreq, step);
    return due;
}

opalx::moments::Request DataSink::requiredMoments(const StepOutputs& due, bool hasSpin) const {
    opalx::moments::Request request;
    if (due.stepLine) {
        request |= opalx::moments::Moment::MeanMomentum;
    }
    if (due.statistics) {
        request |= StatWriter::requiredMoments(hasSpin);
    }
    if (due.phaseSpace && options_m.enableHDF5) {
        request |= H5PartWrapperForPT::requiredMoments;
    }
    return request;
}

std::vector<opalx::moments::BeamMoments> DataSink::momentsOfContainers(
        const PartBunch_t& bunch, const StepOutputs& due,
        std::span<const opalx::moments::BeamMoments> computed) const {
    const auto& containers = bunch.getParticleContainers();
    std::vector<opalx::moments::BeamMoments> moments(containers.size());
    for (size_t i = 0; i < containers.size(); ++i) {
        if (!containers[i] || containers[i]->getTotalNum() == 0) {
            continue;
        }
        const opalx::moments::Request request = requiredMoments(due, containers[i]->hasSpin());
        if (request.empty()) {
            continue;
        }
        const bool reuse = i < computed.size() && computed[i].available().contains(request);
        moments[i]       = reuse ? computed[i] : opalx::moments::compute(*containers[i], request);
    }
    return moments;
}

void DataSink::writeStep(
        PartBunch_t& bunch, const StepOutputs& due,
        std::span<const opalx::moments::BeamMoments> computed,
        const std::function<ReferenceFields()>& referenceFields, int reportedBinCount) {
    OPALTimer::Timer timer;
    const size_t totalAll      = bunch.getTotalNumAllContainers();
    const long long globalStep = bunch.getGlobalTrackStep();

    if (totalAll == 0 && due.stepLine) {
        *gmsg << level1 << "* " << timer.time() << " "
              << "Step " << std::setw(6) << globalStep << "; "
              << "   -- no emission yet --     "
              << "t= " << Util::getTimeString(bunch.getT()) << endl;
    }

    // One request per container covers the step line and every dump due now.
    const auto moments = momentsOfContainers(bunch, due, computed);

    bool anyNonEmpty       = false;
    const auto& containers = bunch.getParticleContainers();
    for (size_t ci = 0; ci < containers.size(); ++ci) {
        const auto& pc = containers[ci];
        if (!pc || pc->getTotalNum() == 0) {
            continue;
        }
        const double sPos = pc->get_sPos();
        if (std::isnan(sPos) || std::isinf(sPos)) {
            throw OpalException(
                    "DataSink::writeStep",
                    "invalid path length s for particle container " + std::to_string(ci));
        }
        if (due.stepLine) {
            *gmsg << level1 << "* " << timer.time() << " "
                  << "Step " << std::setw(6) << globalStep << " "
                  << "container[" << ci << "] ";
            *gmsg << "at " << Util::getLengthString(sPos) << ", ";
            *gmsg << "t= " << Util::getTimeString(bunch.getT()) << ", "
                  << "E=" << Util::getEnergyString(moments[ci].meanKineticEnergy()) << endl;
        }
        anyNonEmpty = true;
    }

    if (anyNonEmpty || due.statistics) {
        writeDumps(bunch, due, moments, referenceFields, reportedBinCount);
    }

    if (due.rankLoad) {
        for (size_t i = 0; i < containers.size(); ++i) {
            const auto& pc = containers[i];
            if (!pc || !bunch.isPcActive(i)) {
                continue;
            }
            pc->printRankLoadInfo("container[" + std::to_string(i) + "]");
        }
    }
}

void DataSink::writeLastStep(
        PartBunch_t& bunch, const std::function<ReferenceFields()>& referenceFields,
        int reportedBinCount) {
    // The step loop wrote the last step already if a dump was due there.
    const long long lastStep = bunch.getGlobalTrackStep() - 1;
    StepOutputs due;
    due.phaseSpace = options_m.psDumpFreq > 0 && !outputDue(options_m.psDumpFreq, lastStep);
    due.statistics = options_m.statDumpFreq > 0 && !outputDue(options_m.statDumpFreq, lastStep);

    writeDumps(bunch, due, momentsOfContainers(bunch, due, {}), referenceFields, reportedBinCount);

    if (due.phaseSpace) {
        *gmsg << level2 << "* Dump phase space of last step" << endl;
    }
}

void DataSink::writeDumps(
        PartBunch_t& bunch, const StepOutputs& due,
        std::span<const opalx::moments::BeamMoments> moments,
        const std::function<ReferenceFields()>& referenceFields, int reportedBinCount) {
    if (!due.phaseSpace && !due.statistics) {
        return;
    }
    const ReferenceFields fields = referenceFields();

    if (due.statistics) {
        dumpSDDS(bunch, fields, moments, reportedBinCount);
        *gmsg << level3 << "* Wrote beam statistics." << endl;
    }

    if (due.phaseSpace && bunch.getTotalNumAllContainers() > 0) {
        dumpH5(bunch, fields, moments);
        *gmsg << level3 << "* Wrote beam phase space." << endl;
    }
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

        Vector_t<double, 3> fd[2] = {fdextPerContainer[i][0], fdextPerContainer[i][1]};
        // Particles outside BEAMHALOBOUNDARY are not counted; partsOutside keeps the file layout.
        statWriters_m[i]->write(beam, moments[i], fd, reportedBinCount, 0, i);
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
