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

//   Class also writes processor load balancing data to file to track parallel
//   calculation efficiency.
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

#include "AbstractObjects/OpalData.h"
#include "Fields/Fieldmap.h"
#include "Physics/Units.h"
#include "Structure/BoundaryGeometry.h"
#include "Structure/H5PartWrapper.h"
#include "Structure/LBalWriter.h"
#include "Utilities/Options.h"
#include "Utilities/Timer.h"
#include "Utilities/Util.h"

#include <algorithm>
#include <sstream>

std::string DataSink::diagnosticStemForContainer(
    const std::string& inputBasename, size_t numContainers, size_t index) {
    if (numContainers <= 1) {
        return inputBasename;
    }
    return inputBasename + "_c" + std::to_string(index);
}

DataSink::DataSink() {
    this->init(false, {}, 1);
}

DataSink::DataSink(const std::vector<H5PartWrapper*>& h5wrappers, bool restart, size_t numParticleContainers) {
    if (restart && !Options::enableHDF5) {
        throw OpalException("DataSink::DataSink()", "Can not restart when HDF5 is disabled");
    }

    this->init(restart, h5wrappers, numParticleContainers);

    if (restart) {
        rewindLines();
    }
}

DataSink::DataSink(H5PartWrapper* h5wrapper, bool restart)
    : DataSink(std::vector<H5PartWrapper*>{h5wrapper}, restart, 1) {
}

DataSink::DataSink(H5PartWrapper* h5wrapper) : DataSink(h5wrapper, false) {
}

void DataSink::dumpH5(const std::shared_ptr<PartBunch_t>& beam, Vector_t<double, 3> FDext[]) const {
    const size_t n = beam->getNumParticleContainers();
    std::vector<std::array<Vector_t<double, 3>, 2>> v(n);
    for (size_t i = 0; i < n; ++i) {
        v[i][0] = FDext[0];
        v[i][1] = FDext[1];
    }
    dumpH5(beam, v);
}

void DataSink::dumpH5(
    const std::shared_ptr<PartBunch_t>& beam,
    const std::vector<std::array<Vector_t<double, 3>, 2>>& fdextPerContainer) const {
    if (!Options::enableHDF5) {
        return;
    }

    const size_t n = beam->getNumParticleContainers();
    if (fdextPerContainer.size() < n) {
        throw OpalException(
            "DataSink::dumpH5",
            "fdextPerContainer size (" + std::to_string(fdextPerContainer.size()) + ") < num containers ("
                + std::to_string(n) + ").");
    }

    for (size_t i = 0; i < n; ++i) {
        if (i >= h5Writers_m.size()) {
            break;
        }
        auto pc = beam->getParticleContainer(i);
        if (!pc || pc->getTotalNum() == 0) {
            continue;
        }
        Vector_t<double, 3> fd[2] = {fdextPerContainer[i][0], fdextPerContainer[i][1]};
        h5Writers_m[i]->writePhaseSpace(beam.get(), fd, i);
    }
}

int DataSink::dumpH5(
    const std::shared_ptr<PartBunch_t>& beam, Vector_t<double, 3> FDext[], double meanEnergy,
    double refPr, double refPt,
    double refPz, double refR, double refTheta, double refZ, double azimuth, double elevation,
    bool local) const {
    if (!Options::enableHDF5 || h5Writers_m.empty()) {
        return -1;
    }

    return h5Writers_m[0]->writePhaseSpace(
        beam.get(), FDext, meanEnergy, refPr, refPt, refPz, refR, refTheta, refZ, azimuth, elevation,
        local, 0);
}

void DataSink::dumpSDDS(
    const std::shared_ptr<PartBunch_t>& beam,
    const std::vector<std::array<Vector_t<double, 3>, 2>>& fdextPerContainer,
    const double& azimuth) const {
    dumpSDDS(beam, fdextPerContainer, losses_t(), azimuth);
}

void DataSink::dumpSDDS(
    const std::shared_ptr<PartBunch_t>& beam,
    const std::vector<std::array<Vector_t<double, 3>, 2>>& fdextPerContainer,
    const losses_t& losses,
    const double& azimuth) const {

    const size_t n = beam->getNumParticleContainers();
    if (fdextPerContainer.size() < n) {
        throw OpalException(
            "DataSink::dumpSDDS",
            "fdextPerContainer size (" + std::to_string(fdextPerContainer.size()) + ") < num containers ("
                + std::to_string(n) + ").");
    }

    IpplTimings::startTimer(StatMarkerTimer_m);

    for (size_t i = 0; i < n; ++i) {
        if (i >= statWriters_m.size()) {
            break;
        }
        auto pc = beam->getParticleContainer(i);
        if (!pc || pc->getTotalNum() == 0) {
            continue;
        }

        pc->updateMoments();

        size_t npOutside = 0;
        if (Options::beamHaloBoundary > 0) {
            npOutside = beam->calcNumPartsOutside(
                Options::beamHaloBoundary * pc->getRmsR());
        }

        Vector_t<double, 3> fd[2] = {fdextPerContainer[i][0], fdextPerContainer[i][1]};
        statWriters_m[i]->write(beam.get(), fd, losses, azimuth, npOutside, i);
    }

    beam->gatherLoadBalanceStatistics();
    beam->calcBeamParameters();

    IpplTimings::stopTimer(StatMarkerTimer_m);
}

void DataSink::storeCavityInformation() {
    if (!Options::enableHDF5 || h5Writers_m.empty()) {
        return;
    }

    h5Writers_m[0]->storeCavityInformation();
}

void DataSink::changeH5Wrappers(const std::vector<H5PartWrapper*>& h5wrappers) {
    if (!Options::enableHDF5) {
        return;
    }

    if (h5wrappers.size() != h5Writers_m.size()) {
        throw OpalException(
            "DataSink::changeH5Wrappers",
            "Expected " + std::to_string(h5Writers_m.size()) + " wrappers, got "
                + std::to_string(h5wrappers.size()) + ".");
    }

    for (size_t i = 0; i < h5wrappers.size(); ++i) {
        h5Writers_m[i]->changeH5Wrapper(h5wrappers[i]);
    }
}

void DataSink::writeGeomToVtk(BoundaryGeometry& bg, std::string fn) {
    if (ippl::Comm->rank() == 0 && Options::enableVTK) {
        bg.writeGeomToVtk(fn);
    }
}

void DataSink::writeImpactStatistics(
    const PartBunch_t* beam, long long& step, size_t& impact, double& sey_num,
    size_t numberOfFieldEmittedParticles, bool nEmissionMode, std::string fn) {
    double charge  = 0.0;
    size_t Npart   = 0;
    double Npart_d = 0.0;
    if (!nEmissionMode) {
        const auto pc = beam->getParticleContainers()[0];
        charge = -1.0 * pc->getTotalCharge();
        // reduce(charge, charge, OpAddAssign());
        Npart_d = -1.0 * charge / pc->getChargePerParticle();
    } else {
        Npart = beam->getParticleContainers()[0]->getTotalNum();
    }
    if (ippl::Comm->rank() == 0) {
        std::string ffn = fn + std::string(".dat");

        std::unique_ptr<Inform> ofp(new Inform(nullptr, ffn.c_str(), Inform::APPEND, 0));
        Inform& fid = *ofp;
        fid.precision(6);
        fid << std::setiosflags(std::ios::scientific);
        double t = beam->getT() * Units::s2ns;
        if (!nEmissionMode) {
            if (step == 0) {
                fid << "#Time/ns" << std::setw(18) << "#Geometry impacts" << std::setw(18)
                    << "tot_sey" << std::setw(18) << "TotalCharge" << std::setw(18) << "PartNum"
                    << " numberOfFieldEmittedParticles " << endl;
            }
            fid << t << std::setw(18) << impact << std::setw(18) << sey_num << std::setw(18)
                << charge << std::setw(18) << Npart_d << std::setw(18)
                << numberOfFieldEmittedParticles << endl;
        } else {
            if (step == 0) {
                fid << "#Time/ns" << std::setw(18) << "#Geometry impacts" << std::setw(18)
                    << "tot_sey" << std::setw(18) << "ParticleNumber"
                    << " numberOfFieldEmittedParticles " << endl;
            }
            fid << t << std::setw(18) << impact << std::setw(18) << sey_num << std::setw(18)
                << double(Npart) << std::setw(18) << numberOfFieldEmittedParticles << endl;
        }
    }
}

void DataSink::dumpBinConfig(long long step,
                             double time,
                             bool preMerge,
                             const std::vector<std::size_t>& binCounts,
                             const std::vector<double>& binWidths,
                             double xMin,
                             const std::string& fileName) {
    if (ippl::Comm->rank() != 0) {
        return;
    }

    Inform m("DataSink::dumpBinConfig");

    if (binCounts.size() != binWidths.size()) {
        m << level4 << "Invalid bin configuration: binCounts.size() = " << binCounts.size()
          << ", binWidths.size() = " << binWidths.size() << endl;
        throw OpalException("DataSink::dumpBinConfig",
                            "binCounts and binWidths must have the same length.");
    }

    if (!binConfigWriter_m) {
        m << level4 << "Creating BinConfigWriter for JSON binning output file \""
          << fileName << "\"." << endl;
        binConfigWriter_m = std::make_unique<BinConfigWriter>(fileName);
    }

    m << level5 << "Dumping bin configuration snapshot: step=" << step
      << ", time=" << time
      << ", preMerge=" << (preMerge ? 1 : 0)
      << ", nBins=" << binCounts.size()
      << ", xMin=" << xMin << endl;

    binConfigWriter_m->writeEntry(step, time, preMerge, binCounts, binWidths, xMin);
}

void DataSink::rewindLines() {
    unsigned int linesToRewind = 0;

    for (size_t i = 0; i < h5Writers_m.size(); ++i) {
        double spos = h5Writers_m[i]->getLastPosition();
        if (i < statWriters_m.size() && statWriters_m[i]->exists()) {
            unsigned int li = statWriters_m[i]->rewindToSpos(spos);
            statWriters_m[i]->replaceVersionString();
            linesToRewind = std::max(linesToRewind, li);
        }
        h5Writers_m[i]->close();
    }

    if (linesToRewind > 0) {
        for (size_t i = 0; i < sddsWriter_m.size(); ++i) {
            sddsWriter_m[i]->rewindLines(linesToRewind);
            sddsWriter_m[i]->replaceVersionString();
        }
    }
}

void DataSink::init(bool restart, const std::vector<H5PartWrapper*>& h5wrappers, size_t numParticleContainers) {
    std::string fn = OpalData::getInstance()->getInputBasename();

    lossWrCounter_m   = 0;
    StatMarkerTimer_m = IpplTimings::getTimer("Write Stat");

    statWriters_m.clear();
    for (size_t i = 0; i < numParticleContainers; ++i) {
        std::string stem = diagnosticStemForContainer(fn, numParticleContainers, i);
        statWriters_m.push_back(statWriter_t(new StatWriter(stem + std::string(".stat"), restart)));
    }

    sddsWriter_m.clear();
    sddsWriter_m.push_back(sddsWriter_t(new LBalWriter(fn + std::string(".lbal"), restart)));

    h5Writers_m.clear();
    if (Options::enableHDF5) {
        if (h5wrappers.size() != numParticleContainers) {
            throw OpalException(
                "DataSink::init",
                "HDF5 enabled: expected " + std::to_string(numParticleContainers) + " H5PartWrapper(s), got "
                    + std::to_string(h5wrappers.size()) + ".");
        }
        for (H5PartWrapper* w : h5wrappers) {
            h5Writers_m.push_back(h5Writer_t(new H5Writer(w, restart)));
        }
    } else if (!h5wrappers.empty()) {
        throw OpalException("DataSink::init", "H5 wrappers passed while HDF5 is disabled.");
    }
}
