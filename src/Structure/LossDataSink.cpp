//
// Class LossDataSink
//   This class writes file attributes to describe phase space of loss files
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
#include "Structure/LossDataSink.h"

#include "AbstractObjects/OpalData.h"
#include "BuildInfo.h"
#include "Utilities/GeneralOpalException.h"
#include "Utilities/GSLHistogram.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"
#include "Utility/IpplInfo.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <functional>
#include <limits>
#include <numeric>
#include <sstream>
#include <tuple>


extern Inform* gmsg;

void LossDataSink::writeFileAttribString(const char* attribute, const char* value) {
    const h5_int64_t h5err = H5WriteFileAttribString(H5file_m, attribute, value);

    if (h5err <= H5_ERR) {
        std::stringstream ss;
        ss << "failed to write string attribute " << attribute << " to file " << fileName_m;
        throw GeneralOpalException("LossDataSink::writeFileAttribString", ss.str());
    }
}

void LossDataSink::writeStepAttribFloat64(
    const char* attribute,
    const h5_float64_t* value,
    h5_int64_t size) {
    const h5_int64_t h5err = H5WriteStepAttribFloat64(H5file_m, attribute, value, size);

    if (h5err <= H5_ERR) {
        std::stringstream ss;
        ss << "failed to write float64 attribute " << attribute << " to file " << fileName_m;
        throw GeneralOpalException("LossDataSink::writeStepAttribFloat64", ss.str());
    }
}

void LossDataSink::writeStepAttribInt64(
    const char* attribute,
    const h5_int64_t* value,
    h5_int64_t size) {
    const h5_int64_t h5err = H5WriteStepAttribInt64(H5file_m, attribute, value, size);

    if (h5err <= H5_ERR) {
        std::stringstream ss;
        ss << "failed to write int64 attribute " << attribute << " to file " << fileName_m;
        throw GeneralOpalException("LossDataSink::writeStepAttribInt64", ss.str());
    }
}

void LossDataSink::writeDataFloat64(const char* name, const h5_float64_t* value) {
    const h5_int64_t h5err = H5PartWriteDataFloat64(H5file_m, name, value);

    if (h5err <= H5_ERR) {
        std::stringstream ss;
        ss << "failed to write float64 data " << name << " to file " << fileName_m;
        throw GeneralOpalException("LossDataSink::writeDataFloat64", ss.str());
    }
}

void LossDataSink::writeDataInt64(const char* name, const h5_int64_t* value) {
    const h5_int64_t h5err = H5PartWriteDataInt64(H5file_m, name, value);

    if (h5err <= H5_ERR) {
        std::stringstream ss;
        ss << "failed to write int64 data " << name << " to file " << fileName_m;
        throw GeneralOpalException("LossDataSink::writeDataInt64", ss.str());
    }
}

void LossDataSink::setStep() {
    const h5_int64_t h5err = H5SetStep(H5file_m, H5call_m);

    if (h5err <= H5_ERR) {
        std::stringstream ss;
        ss << "failed to set step " << H5call_m << " in file " << fileName_m;
        throw GeneralOpalException("LossDataSink::setStep", ss.str());
    }
}

void LossDataSink::getNumSteps() {
    H5call_m = H5GetNumSteps(H5file_m);

    if (H5call_m <= H5_ERR) {
        std::stringstream ss;
        ss << "failed to get number of steps of file " << fileName_m;
        throw GeneralOpalException("LossDataSink::getNumSteps", ss.str());
    }
}

void LossDataSink::setNumParticles(h5_int64_t num) {
    const h5_int64_t h5err = H5PartSetNumParticles(H5file_m, num);

    if (h5err <= H5_ERR) {
        std::stringstream ss;
        ss << "failed to set number of particles to " << num << " in step " << H5call_m
           << " in file " << fileName_m;
        throw GeneralOpalException("LossDataSink::setNumParticles", ss.str());
    }
}

void LossDataSink::openFile(const char* fname, h5_int32_t mode, h5_prop_t props) {
    H5file_m = H5OpenFile(fname, mode, props);

    if (H5file_m == static_cast<h5_file_t>(H5_ERR)) {
        std::stringstream ss;
        ss << "failed to open file " << fileName_m;
        throw GeneralOpalException("LossDataSink::openFile", ss.str());
    }
}

void LossDataSink::closeFile() {
    const h5_int64_t h5err = H5CloseFile(H5file_m);

    if (h5err <= H5_ERR) {
        std::stringstream ss;
        ss << "failed to close file " << fileName_m;
        throw GeneralOpalException("LossDataSink::closeFile", ss.str());
    }

    H5file_m = 0;
}

namespace {
    constexpr double percentileOneSigmaNormalDist =
        0.6826894921370859;

    constexpr double percentileTwoSigmasNormalDist =
        0.9544997361036416;

    constexpr double percentileThreeSigmasNormalDist =
        0.9973002039367398;

    constexpr double percentileFourSigmasNormalDist =
        0.9999366575163338;

    int percentileParticleCount(unsigned long nTotal, double fraction) {
        const double value =
            std::floor(static_cast<double>(nTotal) * fraction + 0.5);

        if (value > static_cast<double>(std::numeric_limits<int>::max())) {
            std::stringstream ss;
            ss << "Percentile particle count " << value
               << " does not fit into int.";

            throw GeneralOpalException(
                "LossDataSink::computeSetPercentiles",
                ss.str());
        }

        return static_cast<int>(value);
    }

    using OneDPhaseSpace_t = std::array<double, 2>;
    using OneDIterator_t   = std::vector<OneDPhaseSpace_t>::iterator;

    void f64transform(
        const std::vector<OpalParticle>& particles,
        size_t startIdx,
        size_t numParticles,
        h5_float64_t* buffer,
        std::function<h5_float64_t(const OpalParticle&)> select) {
        std::transform(
            particles.begin() + startIdx, 
            particles.begin() + startIdx + numParticles, 
            buffer,
            select);
    }

    void i64transform(
        const std::vector<OpalParticle>& particles,
        size_t startIdx,
        size_t numParticles,
        h5_int64_t* buffer,
        std::function<h5_int64_t(const OpalParticle&)> select) {
        std::transform(
            particles.begin() + startIdx, 
            particles.begin() + startIdx + numParticles, 
            buffer,
            select);
    }

    std::pair<double, OneDIterator_t> determinePercentileRange(
        const OneDIterator_t& begin,
        const OneDIterator_t& end,
        const std::vector<int>& globalAccumulatedHistogram,
        const std::vector<int>& localAccumulatedHistogram,
        unsigned int dimension,
        int numRequiredParticles,
        double meanR) {
        const unsigned int numBins = globalAccumulatedHistogram.size() / 3;

        double percentile = 0.0;
        OneDIterator_t endPercentile = end;

        for (unsigned int i = 1; i < numBins; ++i) {
            const unsigned int idx = dimension * numBins + i;

            if (globalAccumulatedHistogram[idx] > numRequiredParticles) {
                OneDIterator_t beginBin = begin + localAccumulatedHistogram[idx - 1];
                OneDIterator_t endBin   = begin + localAccumulatedHistogram[idx];

                int numMissingParticles =
                    numRequiredParticles - globalAccumulatedHistogram[idx - 1];

                unsigned int shift = 2;
                while (numMissingParticles == 0 && idx >= shift) {
                    beginBin = begin + localAccumulatedHistogram[idx - shift];
                    numMissingParticles =
                        numRequiredParticles - globalAccumulatedHistogram[idx - shift];
                    ++shift;
                }

                std::vector<unsigned int> numParticlesInBin(ippl::Comm->size() + 1, 0);
                numParticlesInBin[ippl::Comm->rank() + 1] =
                    static_cast<unsigned int>(endBin - beginBin);

                ippl::Comm->allreduce(
                    &(numParticlesInBin[1]),
                    ippl::Comm->size(),
                    std::plus<unsigned int>());

                std::partial_sum(
                    numParticlesInBin.begin(),
                    numParticlesInBin.end(),
                    numParticlesInBin.begin());

                std::vector<double> positions(numParticlesInBin.back(), 0.0);

                std::transform(
                    beginBin,
                    endBin,
                    positions.begin() + numParticlesInBin[ippl::Comm->rank()],
                    [meanR](const OneDPhaseSpace_t& particle) {
                        return std::abs(particle[0] - meanR);
                    });

                if (!positions.empty()) {
                    ippl::Comm->allreduce(
                        positions.data(),
                        positions.size(),
                        std::plus<double>());

                    std::sort(positions.begin(), positions.end());

                    const auto rhs = static_cast<std::size_t>(
                        std::min<int>(
                            numMissingParticles,
                            static_cast<int>(positions.size()) - 1));

                    const auto lhs = rhs == 0 ? 0 : rhs - 1;

                    percentile = 0.5 * (positions[lhs] + positions[rhs]);
                }

                for (OneDIterator_t it = beginBin; it != endBin; ++it) {
                    if (std::abs((*it)[0] - meanR) > percentile) {
                        return std::make_pair(percentile, it);
                    }
                }

                endPercentile = endBin;
                break;
            }
        }

        return std::make_pair(percentile, endPercentile);
    }

    double computeNormalizedEmittance(
        const OneDIterator_t& begin,
        const OneDIterator_t& end) {
        double localStatistics[] = {0.0, 0.0, 0.0, 0.0};

        localStatistics[0] = static_cast<double>(end - begin);

        for (OneDIterator_t it = begin; it < end; ++it) {
            const OneDPhaseSpace_t& rp = *it;

            localStatistics[1] += rp[0];
            localStatistics[2] += rp[1];
            localStatistics[3] += rp[0] * rp[1];
        }

        ippl::Comm->allreduce(localStatistics, 4, std::plus<double>());

        const double numParticles = localStatistics[0];

        if (numParticles == 0.0) {
            return 0.0;
        }

        const double perParticle = 1.0 / numParticles;
        const double meanR       = localStatistics[1] * perParticle;
        const double meanP       = localStatistics[2] * perParticle;
        const double RP          = localStatistics[3] * perParticle;

        double varianceStatistics[] = {0.0, 0.0};

        for (OneDIterator_t it = begin; it < end; ++it) {
            const OneDPhaseSpace_t& rp = *it;

            varianceStatistics[0] += std::pow(rp[0] - meanR, 2);
            varianceStatistics[1] += std::pow(rp[1] - meanP, 2);
        }

        ippl::Comm->allreduce(varianceStatistics, 2, std::plus<double>());

        const double stdR = std::sqrt(varianceStatistics[0] / numParticles);
        const double stdP = std::sqrt(varianceStatistics[1] / numParticles);

        const double sumRP      = RP - meanR * meanP;
        const double squaredEps = std::pow(stdR * stdP, 2) - std::pow(sumRP, 2);

        return std::sqrt(std::max(squaredEps, 0.0));
    }

    void computeSetPercentiles(
        const std::vector<OpalParticle>& particles,
        size_t startIdx,
        size_t nLoc,
        SetStatistics& stat) {
        if (!Options::computePercentiles || stat.nTotal_m < 100) {
            return;
        }

        unsigned int numBins = static_cast<unsigned int>(
            3.5 * std::pow(3.0, std::log10(static_cast<double>(stat.nTotal_m))));

        numBins = std::max(1u, numBins);

        std::vector<gsl_histogram*> histograms(3, nullptr);

        Vector_t<double, 3> maxR(0.0);

        for (unsigned int d = 0; d < 3; ++d) {
            maxR(d) = 1.0000001
                * std::max(
                    stat.rmax_m(d) - stat.rmean_m(d),
                    stat.rmean_m(d) - stat.rmin_m(d));

            if (maxR(d) <= 0.0) {
                maxR(d) = 1.0e-30;
            }

            histograms[d] = gsl_histogram_alloc(numBins);
            gsl_histogram_set_ranges_uniform(histograms[d], 0.0, maxR(d));
        }

        for (size_t k = 0; k < nLoc; ++k) {
            const OpalParticle& particle = particles[startIdx + k];

            for (unsigned int d = 0; d < 3; ++d) {
                gsl_histogram_increment(
                    histograms[d],
                    std::abs(particle[2 * d] - stat.rmean_m(d)));
            }
        }

        std::vector<int> localHistogramValues(3 * (numBins + 1), 0);
        std::vector<int> globalHistogramValues(3 * (numBins + 1), 0);

        for (unsigned int d = 0; d < 3; ++d) {
            int accumulated = 0;

            for (unsigned int j = 0; j < numBins; ++j) {
                accumulated += static_cast<int>(gsl_histogram_get(histograms[d], j));
                localHistogramValues[d * (numBins + 1) + j + 1] = accumulated;
            }

            gsl_histogram_free(histograms[d]);
        }

        ippl::Comm->allreduce(
            localHistogramValues.data(),
            globalHistogramValues.data(),
            3 * (numBins + 1),
            std::plus<int>());

        const int numParticles68 =
            percentileParticleCount(stat.nTotal_m, percentileOneSigmaNormalDist);

        const int numParticles95 =
            percentileParticleCount(stat.nTotal_m, percentileTwoSigmasNormalDist);

        const int numParticles99 =
            percentileParticleCount(stat.nTotal_m, percentileThreeSigmasNormalDist);

        const int numParticles99_99 =
            percentileParticleCount(stat.nTotal_m, percentileFourSigmasNormalDist);

        for (unsigned int d = 0; d < 3; ++d) {
            std::vector<OneDPhaseSpace_t> oneDPhaseSpace(nLoc);

            for (size_t k = 0; k < nLoc; ++k) {
                const OpalParticle& particle = particles[startIdx + k];

                oneDPhaseSpace[k][0] = particle[2 * d];
                oneDPhaseSpace[k][1] = particle[2 * d + 1];
            }

            std::sort(
                oneDPhaseSpace.begin(),
                oneDPhaseSpace.end(),
                [d, &stat](const OneDPhaseSpace_t& left, const OneDPhaseSpace_t& right) {
                    return std::abs(left[0] - stat.rmean_m(d))
                         < std::abs(right[0] - stat.rmean_m(d));
                });

            OneDIterator_t endSixtyEight = oneDPhaseSpace.end();
            OneDIterator_t endNinetyFive = oneDPhaseSpace.end();
            OneDIterator_t endNinetyNine = oneDPhaseSpace.end();
            OneDIterator_t endNinetyNine_NinetyNine = oneDPhaseSpace.end();

            std::tie(stat.sixtyEightPercentile_m(d), endSixtyEight) =
                determinePercentileRange(
                    oneDPhaseSpace.begin(),
                    oneDPhaseSpace.end(),
                    globalHistogramValues,
                    localHistogramValues,
                    d,
                    numParticles68,
                    stat.rmean_m(d));

            std::tie(stat.ninetyFivePercentile_m(d), endNinetyFive) =
                determinePercentileRange(
                    oneDPhaseSpace.begin(),
                    oneDPhaseSpace.end(),
                    globalHistogramValues,
                    localHistogramValues,
                    d,
                    numParticles95,
                    stat.rmean_m(d));

            std::tie(stat.ninetyNinePercentile_m(d), endNinetyNine) =
                determinePercentileRange(
                    oneDPhaseSpace.begin(),
                    oneDPhaseSpace.end(),
                    globalHistogramValues,
                    localHistogramValues,
                    d,
                    numParticles99,
                    stat.rmean_m(d));

            std::tie(stat.ninetyNine_NinetyNinePercentile_m(d), endNinetyNine_NinetyNine) =
                determinePercentileRange(
                    oneDPhaseSpace.begin(),
                    oneDPhaseSpace.end(),
                    globalHistogramValues,
                    localHistogramValues,
                    d,
                    numParticles99_99,
                    stat.rmean_m(d));

            stat.normalizedEps68Percentile_m(d) =
                computeNormalizedEmittance(oneDPhaseSpace.begin(), endSixtyEight);

            stat.normalizedEps95Percentile_m(d) =
                computeNormalizedEmittance(oneDPhaseSpace.begin(), endNinetyFive);

            stat.normalizedEps99Percentile_m(d) =
                computeNormalizedEmittance(oneDPhaseSpace.begin(), endNinetyNine);

            stat.normalizedEps99_99Percentile_m(d) =
                computeNormalizedEmittance(oneDPhaseSpace.begin(), endNinetyNine_NinetyNine);
        }
    }

} // namespace


SetStatistics::SetStatistics()
    : outputName_m(""),
      spos_m(0.0),
      refTime_m(0.0),
      tmean_m(0.0),
      trms_m(0.0),
      totalCharge_m(0.0),
      totalMass_m(0.0),
      meanKineticEnergy_m(0.0),
      stdKineticEnergy_m(0.0),
      meanGamma_m(0.0),
      nTotal_m(0),
      RefPartR_m(0.0),
      RefPartP_m(0.0),
      rmin_m(0.0),
      rmax_m(0.0),
      rmean_m(0.0),
      pmean_m(0.0),
      rrms_m(0.0),
      prms_m(0.0),
      rprms_m(0.0),
      normEmit_m(0.0),
      geomEmit_m(0.0),
      maxR_m(0.0),
      rsqsum_m(0.0),
      psqsum_m(0.0),
      rpsum_m(0.0),
      eps2_m(0.0),
      eps_norm_m(0.0),
      fac_m(0.0),
      sixtyEightPercentile_m(0.0),
      ninetyFivePercentile_m(0.0),
      ninetyNinePercentile_m(0.0),
      ninetyNine_NinetyNinePercentile_m(0.0),
      normalizedEps68Percentile_m(0.0),
      normalizedEps95Percentile_m(0.0),
      normalizedEps99Percentile_m(0.0),
      normalizedEps99_99Percentile_m(0.0) {
}

LossDataSink::LossDataSink(std::string outfn, bool hdf5Save, CollectionType collectionType)
    : h5hut_mode_m(hdf5Save),
      H5file_m(0),
      outputName_m(outfn),
      H5call_m(0),
      collectionType_m(collectionType) {
    particles_m.clear();
    turnNumber_m.clear();
    bunchNumber_m.clear();

    if (h5hut_mode_m && !Options::enableHDF5) {
        throw GeneralOpalException(
                "LossDataSink::LossDataSink",
                "You must select an OPTION to save Loss data files\n"
                "Please, choose 'ENABLEHDF5=TRUE' or 'ASCIIDUMP=TRUE'");
    }

    OpalData::getInstance()->checkAndAddOutputFileName(outputName_m);

    if (OpalData::getInstance()->hasBunchAllocated()) {
        OpalData::getInstance()->setOpenMode(OpalData::OpenMode::APPEND);
    }
}

LossDataSink::LossDataSink(const LossDataSink& rhs)
    : h5hut_mode_m(rhs.h5hut_mode_m),
      H5file_m(0),
      outputName_m(rhs.outputName_m),
      H5call_m(rhs.H5call_m),
      RefPartR_m(rhs.RefPartR_m),
      RefPartP_m(rhs.RefPartP_m),
      globalTrackStep_m(rhs.globalTrackStep_m),
      refTime_m(rhs.refTime_m),
      spos_m(rhs.spos_m),
      collectionType_m(rhs.collectionType_m) {
    particles_m.clear();
    turnNumber_m.clear();
    bunchNumber_m.clear();
}

LossDataSink::~LossDataSink() noexcept(false) {
    if (H5file_m) {
        closeFile();
    }
    ippl::Comm->barrier();
}

void LossDataSink::openH5(h5_int32_t mode) {
    h5_prop_t props = H5CreateFileProp();
    MPI_Comm comm   = ippl::Comm->getCommunicator();
    H5SetPropFileMPIOCollective(props, &comm);
    openFile(fileName_m.c_str(), mode, props);
    H5CloseProp(props);
}

void LossDataSink::writeHeaderH5() {
    // Write file attributes to describe phase space to H5 file.
    std::stringstream OPAL_version;
    OPAL_version << buildinfo::project_name << " " << buildinfo::project_version << " # git rev. "
                 << Util::getGitRevision();
    writeFileAttribString("OPAL_version", OPAL_version.str().c_str());

    writeFileAttribString("SPOSUnit", "m");
    writeFileAttribString("TIMEUnit", "s");
    writeFileAttribString("RefPartRUnit", "m");
    writeFileAttribString("RefPartPUnit", "#beta#gamma");
    writeFileAttribString("GlobalTrackStepUnit", "1");

    writeFileAttribString("centroidUnit", "m");
    writeFileAttribString("RMSXUnit", "m");
    writeFileAttribString("MEANPUnit", "#beta#gamma");
    writeFileAttribString("RMSPUnit", "#beta#gamma");
    writeFileAttribString("#varepsilonUnit", "m rad");
    writeFileAttribString("#varepsilon-geomUnit", "m rad");
    writeFileAttribString("maxRUnit", "m");
    writeFileAttribString("ENERGYUnit", "MeV");
    writeFileAttribString("dEUnit", "MeV");
    writeFileAttribString("TotalChargeUnit", "C");
    writeFileAttribString("TotalMassUnit", "MeV");

    writeFileAttribString("idUnit", "1");
    writeFileAttribString("xUnit", "m");
    writeFileAttribString("yUnit", "m");
    writeFileAttribString("zUnit", "m");
    writeFileAttribString("pxUnit", "#beta#gamma");
    writeFileAttribString("pyUnit", "#beta#gamma");
    writeFileAttribString("pzUnit", "#beta#gamma");
    writeFileAttribString("qUnit", "C");
    writeFileAttribString("mUnit", "MeV");

    writeFileAttribString("turnUnit", "1");
    writeFileAttribString("bunchNumUnit", "1");

    writeFileAttribString("timeUnit", "s");
    writeFileAttribString("meanTimeUnit", "s");
    writeFileAttribString("rmsTimeUnit", "s");

    if (Options::computePercentiles) {
        writeFileAttribString("68-percentileUnit", "m");
        writeFileAttribString("95-percentileUnit", "m");
        writeFileAttribString("99-percentileUnit", "m");
        writeFileAttribString("99_99-percentileUnit", "m");
        writeFileAttribString("normalizedEps68PercentileUnit", "m rad");
        writeFileAttribString("normalizedEps95PercentileUnit", "m rad");
        writeFileAttribString("normalizedEps99PercentileUnit", "m rad");
        writeFileAttribString("normalizedEps99_99PercentileUnit", "m rad");
    }

    if (collectionType_m == CollectionType::TEMPORAL) {
        writeFileAttribString("type", "temporal");
    } else {
        writeFileAttribString("type", "spatial");
    }
}

void LossDataSink::writeHeaderASCII() {
    bool hasTurn = hasTurnInformations();
    if (ippl::Comm->rank() == 0) {
        os_m << "# x (m),  y (m),  z (m),  px ( ),  py ( ),  pz ( ), id";
        if (hasTurn) {
            os_m << ",  turn ( ), bunchNumber ( )";
        }
        os_m << ", time (s)" << std::endl;
    }
}

void LossDataSink::addReferenceParticle(
        const Vector_t<double, 3>& x, const Vector_t<double, 3>& p, double time, double spos,
        long long globalTrackStep) {
    RefPartR_m.push_back(x);
    RefPartP_m.push_back(p);
    globalTrackStep_m.push_back((h5_int64_t)globalTrackStep);
    spos_m.push_back(spos);
    refTime_m.push_back(time);
}

void LossDataSink::addParticle(
    const OpalParticle& particle,
    const std::optional<std::pair<int, short>>& turnBunchNumPair) {

    if (turnBunchNumPair) {
        if (turnNumber_m.size() != particles_m.size()) {
            throw GeneralOpalException(
                "LossDataSink::addParticle",
                "Either all particles or no particles must have turn number and bunch number");
        }

        turnNumber_m.push_back(turnBunchNumPair.value().first);
        bunchNumber_m.push_back(turnBunchNumPair.value().second);

    } else {
        if (!turnNumber_m.empty()) {
            throw GeneralOpalException(
                "LossDataSink::addParticle",
                "Either all particles or no particles must have turn number and bunch number");
        }
    }

    particles_m.push_back(particle);
}

void LossDataSink::save(unsigned int numSets, OpalData::OpenMode openMode) {
    if (outputName_m.empty()) return;
    if (hasNoParticlesToDump()) return;

    if (openMode == OpalData::OpenMode::UNDEFINED) {
        openMode = OpalData::getInstance()->getOpenMode();
    }

    namespace fs = std::filesystem;
    if (h5hut_mode_m) {
        fileName_m = outputName_m + std::string(".h5");
        if (openMode == OpalData::OpenMode::WRITE || !fs::exists(fileName_m)) {
            openH5();
            writeHeaderH5();
        } else {
            openH5(H5_O_APPENDONLY);
            getNumSteps();
        }

        splitSets(numSets);

        for (unsigned int i = 0; i < numSets; ++i) {
            saveH5(i);
        }
        closeFile();

    } else {
        fileName_m = outputName_m + std::string(".loss");
        if (openMode == OpalData::OpenMode::WRITE || !fs::exists(fileName_m)) {
            openASCII();
            writeHeaderASCII();
        } else {
            appendASCII();
        }
        saveASCII();
        closeASCII();
    }
    *gmsg << level2 << "Save '" << fileName_m << "'" << endl;

    ippl::Comm->barrier();

    particles_m.clear();
    turnNumber_m.clear();
    bunchNumber_m.clear();
    spos_m.clear();
    refTime_m.clear();
    RefPartR_m.clear();
    RefPartP_m.clear();
    globalTrackStep_m.clear();
    startSet_m.clear();
}

// Note: This was changed to calculate the global number of dumped particles
// because there are two cases to be considered:
// 1. ALL nodes have 0 lost particles -> nothing to be done.
// 2. Some nodes have 0 lost particles, some not -> H5 can handle that but all
// nodes HAVE to participate, otherwise H5 waits endlessly for a response from
// the nodes that didn't enter the saveH5 function. -DW
bool LossDataSink::hasNoParticlesToDump() const {
    const unsigned long long nLocal =
        static_cast<unsigned long long>(particles_m.size());

    unsigned long long nGlobal = 0;
    ippl::Comm->allreduce(&nLocal, &nGlobal, 1, std::plus<unsigned long long>());

    return nGlobal == 0;
}

bool LossDataSink::hasTurnInformations() const {
    bool hasTurnInformation = !turnNumber_m.empty();

    ippl::Comm->allreduce(hasTurnInformation, 1, std::logical_or<bool>());

    return hasTurnInformation > 0;
}

void LossDataSink::saveH5(unsigned int setIdx) {
    size_t nLoc     = particles_m.size();
    size_t startIdx = 0;
    size_t endIdx   = nLoc;
    if (setIdx + 1 < startSet_m.size()) {
        startIdx = startSet_m[setIdx];
        endIdx   = startSet_m[setIdx + 1];
        nLoc     = endIdx - startIdx;
    }

    SetStatistics stat = computeSetStatistics(setIdx);

    if (stat.nTotal_m == 0) {
        return;
    }

    const bool hasTurn = hasTurnInformations();

    if (hasTurn
        && (turnNumber_m.size() < startIdx + nLoc
            || bunchNumber_m.size() < startIdx + nLoc)) {
        throw GeneralOpalException(
            "LossDataSink::saveH5",
            "Turn/bunch information is globally present, but this rank does not have "
            "turn/bunch data for all particles in this loss set.");
    }

    /// Set current record/time step.
    setStep();
    setNumParticles(static_cast<h5_int64_t>(nLoc));

    if (setIdx < spos_m.size()) {
        writeStepAttribFloat64("SPOS", &(spos_m[setIdx]), 1);
    }

    if (setIdx < refTime_m.size()) {
        writeStepAttribFloat64("TIME", &(refTime_m[setIdx]), 1);
    }

    if (setIdx < RefPartR_m.size()) {
        Vector_t<double, 3> refPartR = RefPartR_m[setIdx];
        writeStepAttribFloat64("RefPartR", &refPartR[0], 3);
    }

    if (setIdx < RefPartP_m.size()) {
        Vector_t<double, 3> refPartP = RefPartP_m[setIdx];
        writeStepAttribFloat64("RefPartP", &refPartP[0], 3);
    }

    if (setIdx < globalTrackStep_m.size()) {
        writeStepAttribInt64("GlobalTrackStep", &(globalTrackStep_m[setIdx]), 1);
    }

    Vector_t<double, 3> tmpVector;
    double tmpDouble;

    tmpVector = stat.rmean_m;
    writeStepAttribFloat64("centroid", &tmpVector[0], 3);

    tmpVector = stat.rrms_m;
    writeStepAttribFloat64("RMSX", &tmpVector[0], 3);

    tmpVector = stat.pmean_m;
    writeStepAttribFloat64("MEANP", &tmpVector[0], 3);

    tmpVector = stat.prms_m;
    writeStepAttribFloat64("RMSP", &tmpVector[0], 3);

    tmpVector = stat.eps_norm_m;
    writeStepAttribFloat64("#varepsilon", &tmpVector[0], 3);

    tmpVector = stat.geomEmit_m;
    writeStepAttribFloat64("#varepsilon-geom", &tmpVector[0], 3);

    tmpDouble = stat.meanKineticEnergy_m;
    writeStepAttribFloat64("ENERGY", &tmpDouble, 1);

    tmpDouble = stat.stdKineticEnergy_m;
    writeStepAttribFloat64("dE", &tmpDouble, 1);

    tmpDouble = stat.totalCharge_m;
    writeStepAttribFloat64("TotalCharge", &tmpDouble, 1);

    tmpDouble = stat.totalMass_m;
    writeStepAttribFloat64("TotalMass", &tmpDouble, 1);

    tmpDouble = stat.tmean_m;
    writeStepAttribFloat64("meanTime", &tmpDouble, 1);

    tmpDouble = stat.trms_m;
    writeStepAttribFloat64("rmsTime", &tmpDouble, 1);

    if (Options::computePercentiles) {
        tmpVector = stat.sixtyEightPercentile_m;
        writeStepAttribFloat64("68-percentile", &tmpVector[0], 3);

        tmpVector = stat.ninetyFivePercentile_m;
        writeStepAttribFloat64("95-percentile", &tmpVector[0], 3);

        tmpVector = stat.ninetyNinePercentile_m;
        writeStepAttribFloat64("99-percentile", &tmpVector[0], 3);

        tmpVector = stat.ninetyNine_NinetyNinePercentile_m;
        writeStepAttribFloat64("99_99-percentile", &tmpVector[0], 3);

        tmpVector = stat.normalizedEps68Percentile_m;
        writeStepAttribFloat64("normalizedEps68Percentile", &tmpVector[0], 3);

        tmpVector = stat.normalizedEps95Percentile_m;
        writeStepAttribFloat64("normalizedEps95Percentile", &tmpVector[0], 3);

        tmpVector = stat.normalizedEps99Percentile_m;
        writeStepAttribFloat64("normalizedEps99Percentile", &tmpVector[0], 3);

        tmpVector = stat.normalizedEps99_99Percentile_m;
        writeStepAttribFloat64("normalizedEps99_99Percentile", &tmpVector[0], 3);
    }

    tmpVector = stat.maxR_m;
    writeStepAttribFloat64("maxR", &tmpVector[0], 3);

    // Write all data
    std::vector<h5_float64_t> f64bufferStorage(std::max<size_t>(nLoc, 1));
    std::vector<h5_int64_t> i64bufferStorage(std::max<size_t>(nLoc, 1));

    h5_float64_t* f64buffer = f64bufferStorage.data();
    h5_int64_t* i64buffer   = i64bufferStorage.data();

    ::i64transform(particles_m, startIdx, nLoc, i64buffer, [](const OpalParticle& particle) {
        return particle.getId();
    });
    writeDataInt64("id", i64buffer);

    ::f64transform(particles_m, startIdx, nLoc, f64buffer, [](const OpalParticle& particle) {
        return particle.getX();
    });
    writeDataFloat64("x", f64buffer);

    ::f64transform(particles_m, startIdx, nLoc, f64buffer, [](const OpalParticle& particle) {
        return particle.getY();
    });
    writeDataFloat64("y", f64buffer);

    ::f64transform(particles_m, startIdx, nLoc, f64buffer, [](const OpalParticle& particle) {
        return particle.getZ();
    });
    writeDataFloat64("z", f64buffer);

    ::f64transform(particles_m, startIdx, nLoc, f64buffer, [](const OpalParticle& particle) {
        return particle.getPx();
    });
    writeDataFloat64("px", f64buffer);

    ::f64transform(particles_m, startIdx, nLoc, f64buffer, [](const OpalParticle& particle) {
        return particle.getPy();
    });
    writeDataFloat64("py", f64buffer);

    ::f64transform(particles_m, startIdx, nLoc, f64buffer, [](const OpalParticle& particle) {
        return particle.getPz();
    });
    writeDataFloat64("pz", f64buffer);

    ::f64transform(particles_m, startIdx, nLoc, f64buffer, [](const OpalParticle& particle) {
        return particle.getCharge();
    });
    writeDataFloat64("q", f64buffer);

    ::f64transform(particles_m, startIdx, nLoc, f64buffer, [](const OpalParticle& particle) {
        return particle.getMass();
    });
    writeDataFloat64("m", f64buffer);

    if (hasTurn) {
        std::copy(
            turnNumber_m.begin() + startIdx, 
            turnNumber_m.begin() + startIdx + nLoc, 
            i64buffer);
        writeDataInt64("turn", i64buffer);

        std::copy(
            bunchNumber_m.begin() + startIdx, 
            bunchNumber_m.begin() + startIdx + nLoc, 
            i64buffer);
        writeDataInt64("bunchNumber", i64buffer);
    }

    ::f64transform(particles_m, startIdx, nLoc, f64buffer, [](const OpalParticle& particle) {
        return particle.getTime();
    });
    writeDataFloat64("time", f64buffer);

    ++H5call_m;
}

void LossDataSink::saveASCII() {
    const bool hasTurn = hasTurnInformations();

    if (hasTurn
        && (turnNumber_m.size() != particles_m.size()
            || bunchNumber_m.size() != particles_m.size())) {
        throw GeneralOpalException(
            "LossDataSink::saveASCII",
            "Turn/bunch information is globally present, but this rank does not have "
            "turn/bunch data for all particles.");
    }

    constexpr int tagCount  = 43001;
    constexpr int tagDouble = 43002;
    constexpr int tagInt    = 43003;

    MPI_Comm comm = ippl::Comm->getCommunicator();

    const int rank = ippl::Comm->rank();
    const int size = ippl::Comm->size();

    const std::size_t nLoc = particles_m.size();
    const std::size_t nIntColumns = hasTurn ? 3 : 1;

    auto packDoubleData = [&]() {
        std::vector<double> data(7 * nLoc);

        for (std::size_t i = 0; i < nLoc; ++i) {
            const OpalParticle& particle = particles_m[i];

            data[7 * i + 0] = particle.getX();
            data[7 * i + 1] = particle.getY();
            data[7 * i + 2] = particle.getZ();
            data[7 * i + 3] = particle.getPx();
            data[7 * i + 4] = particle.getPy();
            data[7 * i + 5] = particle.getPz();
            data[7 * i + 6] = particle.getTime();
        }

        return data;
    };

    auto packIntegerData = [&]() {
        std::vector<long long> data(nIntColumns * nLoc);

        for (std::size_t i = 0; i < nLoc; ++i) {
            const OpalParticle& particle = particles_m[i];

            data[nIntColumns * i + 0] =
                static_cast<long long>(particle.getId());

            if (hasTurn) {
                data[nIntColumns * i + 1] =
                    static_cast<long long>(turnNumber_m[i]);

                data[nIntColumns * i + 2] =
                    static_cast<long long>(bunchNumber_m[i]);
            }
        }

        return data;
    };

    auto writeRecords =
        [&](const std::vector<double>& doubleData,
            const std::vector<long long>& integerData,
            std::size_t count) {
            for (std::size_t i = 0; i < count; ++i) {
                os_m << doubleData[7 * i + 0] << "   ";
                os_m << doubleData[7 * i + 1] << "   ";
                os_m << doubleData[7 * i + 2] << "   ";
                os_m << doubleData[7 * i + 3] << "   ";
                os_m << doubleData[7 * i + 4] << "   ";
                os_m << doubleData[7 * i + 5] << "   ";
                os_m << integerData[nIntColumns * i + 0] << "   ";

                if (hasTurn) {
                    os_m << integerData[nIntColumns * i + 1] << "   ";
                    os_m << integerData[nIntColumns * i + 2] << "   ";
                }

                os_m << doubleData[7 * i + 6] << std::endl;
            }
        };

    if (rank == 0) {
        std::vector<double> localDoubleData = packDoubleData();
        std::vector<long long> localIntegerData = packIntegerData();

        writeRecords(localDoubleData, localIntegerData, nLoc);

        for (int src = 1; src < size; ++src) {
            unsigned long long remoteCount = 0;

            MPI_Recv(
                &remoteCount,
                1,
                MPI_UNSIGNED_LONG_LONG,
                src,
                tagCount,
                comm,
                MPI_STATUS_IGNORE);

            std::vector<double> remoteDoubleData(7 * remoteCount);
            std::vector<long long> remoteIntegerData(nIntColumns * remoteCount);

            if (remoteCount > 0) {
                MPI_Recv(
                    remoteDoubleData.data(),
                    static_cast<int>(remoteDoubleData.size()),
                    MPI_DOUBLE,
                    src,
                    tagDouble,
                    comm,
                    MPI_STATUS_IGNORE);

                MPI_Recv(
                    remoteIntegerData.data(),
                    static_cast<int>(remoteIntegerData.size()),
                    MPI_LONG_LONG,
                    src,
                    tagInt,
                    comm,
                    MPI_STATUS_IGNORE);
            }

            writeRecords(remoteDoubleData, remoteIntegerData, remoteCount);
        }
    } else {
        const unsigned long long localCount =
            static_cast<unsigned long long>(nLoc);

        MPI_Send(
            &localCount,
            1,
            MPI_UNSIGNED_LONG_LONG,
            0,
            tagCount,
            comm);

        if (localCount > 0) {
            std::vector<double> localDoubleData = packDoubleData();
            std::vector<long long> localIntegerData = packIntegerData();

            MPI_Send(
                localDoubleData.data(),
                static_cast<int>(localDoubleData.size()),
                MPI_DOUBLE,
                0,
                tagDouble,
                comm);

            MPI_Send(
                localIntegerData.data(),
                static_cast<int>(localIntegerData.size()),
                MPI_LONG_LONG,
                0,
                tagInt,
                comm);
        }
    }

    ippl::Comm->barrier();
}

/**
 * In Opal-T monitors can be traversed several times. We know how
 * many times the bunch passes because we register the passage of
 * the reference particle. This code tries to determine to which
 * bunch (same bunch but different times) a particle belongs. For
 * this we could use algorithms from data science such as k-means
 * or dbscan. But they are an overkill for this application because
 * the bunches are well separated.
 *
 * In a first step we a assign to each bunch the same number of
 * particles and compute the mean time of passage and with it a time
 * range. Of course this is only an approximation. So we reassign
 * the particles to the bunches using the time ranges compute a
 * better approximation. Two iterations should be sufficient for
 * Opal-T where the temporal separation is large.
 *
 * @param numSets number of passes of the reference particle
 *
 */
void LossDataSink::splitSets(unsigned int numSets) {
    startSet_m.clear();

    if (numSets <= 1 || particles_m.empty()) {
        return;
    }

    const size_t nLoc   = particles_m.size();
    size_t avgNumPerSet = nLoc / numSets;

    std::vector<size_t> numPartsInSet(numSets, avgNumPerSet);
    for (unsigned int j = 0; j < (nLoc - numSets * avgNumPerSet); ++j) {
        ++numPartsInSet[j];
    }

    startSet_m.resize(numSets + 1, 0);

    std::vector<double> data(2 * numSets, 0.0);
    double* meanT        = data.data();
    double* numParticles = data.data() + numSets;

    std::vector<double> timeRange(numSets, 0.0);

    for (unsigned int iteration = 0; iteration < 2; ++iteration) {
        std::fill(data.begin(), data.end(), 0.0);
        std::fill(timeRange.begin(), timeRange.end(), 0.0);

        double maxT = particles_m[0].getTime();

        size_t partIdx = 0;
        for (unsigned int j = 0; j < numSets; ++j) {
            const size_t numThisSet = numPartsInSet[j];

            for (size_t k = 0; k < numThisSet && partIdx < nLoc; ++k, ++partIdx) {
                meanT[j] += particles_m[partIdx].getTime();
                maxT = std::max(maxT, particles_m[partIdx].getTime());
            }

            numParticles[j] = static_cast<double>(numThisSet);
        }

        ippl::Comm->allreduce(data.data(), 2 * numSets, std::plus<double>());

        for (unsigned int j = 0; j < numSets; ++j) {
            if (numParticles[j] > 0.0) {
                meanT[j] /= numParticles[j];
            }
        }

        for (unsigned int j = 0; j < numSets - 1; ++j) {
            timeRange[j] = 0.5 * (meanT[j] + meanT[j + 1]);
        }
        timeRange[numSets - 1] = maxT;

        std::fill(numPartsInSet.begin(), numPartsInSet.end(), 0);

        size_t setNum   = 0;
        size_t idxPrior = 0;

        for (size_t idx = 0; idx < nLoc && setNum + 1 < numSets; ++idx) {
            if (particles_m[idx].getTime() > timeRange[setNum]) {
                numPartsInSet[setNum] = idx - idxPrior;
                idxPrior              = idx;
                ++setNum;
            }
        }

        numPartsInSet[numSets - 1] = nLoc - idxPrior;
    }

    for (unsigned int i = 0; i < numSets; ++i) {
        startSet_m[i + 1] = startSet_m[i] + numPartsInSet[i];
    }
}

SetStatistics LossDataSink::computeSetStatistics(unsigned int setIdx) {
    SetStatistics stat;
    double part[6];

    // Layout:
    // 0      : number of particles
    // 1..6   : sums of x, px, y, py, z, pz
    // 7..42  : 6x6 second moments
    // 43     : time sum
    // 44     : time^2 sum
    // 45     : total charge
    // 46     : total mass
    // 47     : kinetic energy sum
    // 48     : kinetic energy^2 sum
    // 49     : gamma sum
    const unsigned int totalSize = 50;

    double plainData[totalSize] = {};
    double rMinMax[6];

    for (unsigned int i = 0; i < 6; ++i) {
        rMinMax[i] = std::numeric_limits<double>::lowest();
    }

    Util::KahanAccumulation data[totalSize] = {};
    Util::KahanAccumulation* localCentroid = data + 1;
    Util::KahanAccumulation* localMoments  = data + 7;
    Util::KahanAccumulation* localOthers   = data + 43;

    size_t startIdx = 0;
    size_t nLoc     = particles_m.size();

    if (setIdx + 1 < startSet_m.size()) {
        startIdx = startSet_m[setIdx];
        nLoc     = startSet_m[setIdx + 1] - startSet_m[setIdx];
    }

    data[0].sum = static_cast<double>(nLoc);

    size_t idx = startIdx;

    for (size_t k = 0; k < nLoc; ++k, ++idx) {
        const OpalParticle& particle = particles_m[idx];

        part[0] = particle.getX();
        part[1] = particle.getPx();
        part[2] = particle.getY();
        part[3] = particle.getPy();
        part[4] = particle.getZ();
        part[5] = particle.getPz();

        for (int i = 0; i < 6; ++i) {
            localCentroid[i] += part[i];
            for (int j = 0; j <= i; ++j) {
                localMoments[i * 6 + j] += part[i] * part[j];
            }
        }

        for (unsigned int d = 0; d < 3; ++d) {
            const double r = part[2 * d];

            rMinMax[2 * d]     = std::max(rMinMax[2 * d], -r);
            rMinMax[2 * d + 1] = std::max(rMinMax[2 * d + 1], r);
        }

        const double time = particle.getTime();

        const double gamma = Util::getGamma(particle.getP());
        const double eKin  = Util::getKineticEnergy(particle.getP(), particle.getMass());

        localOthers[0] += time;
        localOthers[1] += time * time;
        localOthers[2] += particle.getCharge();
        localOthers[3] += particle.getMass();
        localOthers[4] += eKin;
        localOthers[5] += eKin * eKin;
        localOthers[6] += gamma;
    }

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < i; ++j) {
            localMoments[j * 6 + i] = localMoments[i * 6 + j];
        }
    }

    for (unsigned int i = 0; i < totalSize; ++i) {
        plainData[i] = data[i].sum;
    }

    ippl::Comm->allreduce(plainData, totalSize, std::plus<double>());
    ippl::Comm->allreduce(rMinMax, 6, std::greater<double>());

    if (plainData[0] == 0.0) {
        return stat;
    }

    double* centroid = plainData + 1;
    double* moments  = plainData + 7;
    double* others   = plainData + 43;

    stat.outputName_m = outputName_m;

    if (setIdx < spos_m.size()) {
        stat.spos_m = spos_m[setIdx];
    }

    if (setIdx < refTime_m.size()) {
        stat.refTime_m = refTime_m[setIdx];
    }

    if (setIdx < RefPartR_m.size()) {
        stat.RefPartR_m = RefPartR_m[setIdx];
    }

    if (setIdx < RefPartP_m.size()) {
        stat.RefPartP_m = RefPartP_m[setIdx];
    }

    stat.nTotal_m = static_cast<unsigned long>(std::round(plainData[0]));

    for (unsigned int i = 0; i < 3u; ++i) {
        stat.rmean_m(i) = centroid[2 * i] / stat.nTotal_m;
        stat.pmean_m(i) = centroid[2 * i + 1] / stat.nTotal_m;

        stat.rsqsum_m(i) =
            moments[(2 * i) * 6 + 2 * i]
            - stat.nTotal_m * std::pow(stat.rmean_m(i), 2);

        stat.psqsum_m(i) =
            moments[(2 * i + 1) * 6 + 2 * i + 1]
            - stat.nTotal_m * std::pow(stat.pmean_m(i), 2);

        stat.psqsum_m(i) = std::max(0.0, stat.psqsum_m(i));

        stat.rpsum_m(i) =
            moments[(2 * i) * 6 + 2 * i + 1]
            - stat.nTotal_m * stat.rmean_m(i) * stat.pmean_m(i);
    }

    stat.tmean_m = others[0] / stat.nTotal_m;
    stat.trms_m  = std::sqrt(
        std::max(0.0, others[1] / stat.nTotal_m - std::pow(stat.tmean_m, 2)));

    stat.totalCharge_m = others[2];
    stat.totalMass_m   = others[3];

    stat.meanKineticEnergy_m = others[4] / stat.nTotal_m;

    Util::KahanAccumulation localKineticEnergyVariance = {};

    size_t varianceIdx = startIdx;
    for (size_t k = 0; k < nLoc; ++k, ++varianceIdx) {
        const OpalParticle& particle = particles_m[varianceIdx];

        const double eKin =
            Util::getKineticEnergy(particle.getP(), particle.getMass());

        const double diff = eKin - stat.meanKineticEnergy_m;
        localKineticEnergyVariance += diff * diff;
    }

    double kineticEnergyVariance[] = {
        static_cast<double>(localKineticEnergyVariance.sum)
    };

    ippl::Comm->allreduce(kineticEnergyVariance, 1, std::plus<double>());

    stat.stdKineticEnergy_m =
        std::sqrt(std::max(0.0, kineticEnergyVariance[0] / stat.nTotal_m));

    stat.meanGamma_m = others[6] / stat.nTotal_m;

    stat.eps2_m =
        (stat.rsqsum_m * stat.psqsum_m - stat.rpsum_m * stat.rpsum_m)
        / (1.0 * stat.nTotal_m * stat.nTotal_m);

    stat.rpsum_m /= stat.nTotal_m;

    const double betaGamma =
        std::sqrt(std::max(0.0, stat.meanGamma_m * stat.meanGamma_m - 1.0));

    for (unsigned int i = 0; i < 3u; ++i) {
        stat.rrms_m(i) = std::sqrt(std::max(0.0, stat.rsqsum_m(i)) / stat.nTotal_m);
        stat.prms_m(i) = std::sqrt(std::max(0.0, stat.psqsum_m(i)) / stat.nTotal_m);

        stat.eps_norm_m(i) = std::sqrt(std::max(0.0, stat.eps2_m(i)));

        stat.geomEmit_m(i) =
            betaGamma > 0.0 ? stat.eps_norm_m(i) / betaGamma : 0.0;

        const double tmp = stat.rrms_m(i) * stat.prms_m(i);
        stat.fac_m(i) = tmp == 0.0 ? 0.0 : 1.0 / tmp;

        stat.rmin_m(i) = -rMinMax[2 * i];
        stat.rmax_m(i) =  rMinMax[2 * i + 1];
        stat.maxR_m(i) =  std::max(std::abs(stat.rmin_m(i)), std::abs(stat.rmax_m(i))); //stat.rmax_m(i);
    }

    stat.rprms_m = stat.rpsum_m * stat.fac_m;

    computeSetPercentiles(particles_m, startIdx, nLoc, stat);

    return stat;
}
