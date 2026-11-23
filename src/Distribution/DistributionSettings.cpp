//
// Class DistributionSettings
//   The values of one DISTRIBUTION that the samplers read, worked out from its description.
//
// Copyright (c) 2008 - 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Distribution/DistributionSettings.h"

#include <cmath>
#include <limits>
#include <map>
#include <utility>
#include <vector>

#include "Utilities/OpalException.h"
#include "Utility/Inform.h"

extern Inform* gmsg;

DistributionSettings::DistributionSettings(
        const ObjectValues& values, double avrgpz, int seed, std::string inputDirectory)
    : seed_m(seed), inputDirectory_m(std::move(inputDirectory)) {
    static const std::map<std::string, DistributionType> typeStringToDistType_s = {
            {"NODIST", DistributionType::NODIST},
            {"UNIFORM", DistributionType::UNIFORM},
            {"GAUSS", DistributionType::GAUSS},
            {"MULTIVARIATEGAUSS", DistributionType::MULTIVARIATEGAUSS},
            {"FLATTOP", DistributionType::FLATTOP},
            {"OPALFLATTOP", DistributionType::OPALFLATTOP},
            {"FROMFILE", DistributionType::FROMFILE},
            {"EMITTEDFROMFILE", DistributionType::EMITTEDFROMFILE}};

    distT_m = values.getString("TYPE");
    if (distT_m.empty()) {
        throw OpalException(
                "Distribution::setDistType",
                "The attribute \"TYPE\" isn't set for the \"DISTRIBUTION\"!");
    }
    distrTypeT_m = typeStringToDistType_s.at(distT_m);

    // This flag is used for every distribution type during checkpoint restart
    // validation. Initialize it from the input before dispatching to the
    // type-specific setup.
    emitting_m = values.getBool("EMITTED");

    switch (distrTypeT_m) {
        case DistributionType::UNIFORM:
            setDistParametersUniform(values);
            break;
        case DistributionType::GAUSS:
            setDistParametersGauss(values);
            break;
        case DistributionType::MULTIVARIATEGAUSS:
            setDistParametersMultiVariateGauss(values);
            break;
        case DistributionType::FLATTOP:
        case DistributionType::OPALFLATTOP:
            setDistParametersFlatTop(values);
            break;
        case DistributionType::FROMFILE:
            // File-based distributions read their records in the sampler class.
            break;
        case DistributionType::EMITTEDFROMFILE:
            // The file stores sampled particle birth times, but not necessarily the full old-OPAL
            // emission window. Parse the flat-top timing attributes when present so the sampler can
            // reconstruct the original phase shift.
            setDistParametersFlatTop(values);
            break;
        default:
            throw OpalException("Distribution Param", "Unknown \"TYPE\" of \"DISTRIBUTION\"");
    }

    totalNumberParticles_m = numParticles(values);
    filename_m             = values.getString("FNAME");
    const double steps     = values.getReal("EMISSIONSTEPS");
    emissionSteps_m        = steps <= 0.0 ? 100 : static_cast<std::size_t>(std::ceil(steps));
    avrgpz_m               = avrgpz;
}

std::size_t DistributionSettings::numParticles(const ObjectValues& values) {
    // A value <= 0 means "unspecified"; the BEAM NALLOC is used then.
    const double npartDist = values.getReal("NPARTDIST");
    return npartDist > 0.0 ? static_cast<std::size_t>(npartDist) : static_cast<std::size_t>(0);
}

void DistributionSettings::setDistParametersUniform(const ObjectValues& values) {
    setSigmaR_m(values);
    for (unsigned d = 0; d < 3; ++d) {
        if (sigmaR_m[d] <= 0.0) {
            throw OpalException(
                    "Distribution::setDistParametersUniform",
                    "UNIFORM requires SIGMAX, SIGMAY and SIGMAZ to be greater than zero.");
        }
    }
}

void DistributionSettings::setDistParametersGauss(const ObjectValues& values) {
    /*
     * Set distribution parameters. Do all the necessary checks depending
     * on the input attributes.
     * In case of DistributionType::MATCHEDGAUSS we only need to set the cutoff parameters
     */

    // Cutoffs in units of sigma. As in OPAL, the transverse cutoffs define an
    // ellipse in (x, y) and in (px, py) while the longitudinal cutoffs are
    // applied per axis (see Gaussian sampler).
    cutoffP_m = ippl::Vector<double, 3>(
            std::abs(values.getReal("CUTOFFPX")), std::abs(values.getReal("CUTOFFPY")),
            std::abs(values.getReal("CUTOFFPZ")));

    cutoffR_m = ippl::Vector<double, 3>(
            std::abs(values.getReal("CUTOFFX")), std::abs(values.getReal("CUTOFFY")),
            std::abs(values.getReal("CUTOFFLONG")));

    // As in OPAL, a (numerically) zero cutoff means "no cut" for CUTOFFX, CUTOFFY,
    // CUTOFFPX, CUTOFFPY and CUTOFFPZ. CUTOFFLONG is excluded: |z| <= 0 would reject
    // every sample (OPAL loops forever in that case), so it must be positive.
    constexpr double smallestCutoff = 1e-12;
    for (unsigned d = 0; d < 3; ++d) {
        if (d != 2 && cutoffR_m[d] < smallestCutoff) {
            cutoffR_m[d] = std::numeric_limits<double>::max();
        }
        if (cutoffP_m[d] < smallestCutoff) {
            cutoffP_m[d] = std::numeric_limits<double>::max();
        }
    }
    if (cutoffR_m[2] < smallestCutoff) {
        throw OpalException(
                "Distribution::setDistParametersGauss",
                "CUTOFFLONG must be greater than zero (a zero cutoff would reject all samples).");
    }

    setSigmaR_m(values);
    setSigmaP_m(values);
}

void DistributionSettings::setDistParametersMultiVariateGauss(const ObjectValues& values) {
    cutoffR_m = 3.;
    cutoffP_m = 3.;

    for (unsigned int i = 0; i < 6; ++i) {
        for (unsigned int j = 0; j < 6; ++j) {
            if (i == j)
                correlationMatrix_m[i][j] = 1.0;
            else
                correlationMatrix_m[i][j] = 0.0;
        }
    }

    setSigmaR_m(values);
    setSigmaP_m(values);

    std::vector<double> cr = values.getRealArray("CORR");

    if (!cr.empty()) {
        // read off-diagonal correlation coefficients from input file
        if (cr.size() == 15) {
            *gmsg << "* Using CORR to specify off-diagonal correlation coefficients" << endl;
            unsigned int k = 0;
            for (unsigned int i = 0; i < 5; ++i) {
                for (unsigned int j = i + 1; j < 6; ++j, ++k) {
                    correlationMatrix_m[j][i] = cr.at(k);
                    correlationMatrix_m[i][j] = cr.at(k);  // impose symmetry
                }
            }
        } else {
            throw OpalException(
                    "Distribution::SetDistParametersGauss",
                    "Inconsistent set of correlations specified, check manual");
        }
    }
}

void DistributionSettings::setDistParametersFlatTop(const ObjectValues& values) {
    cutoffR_m = 3.;
    cutoffP_m = 3.;

    // set diagonal elements first
    setSigmaR_m(values);
    setSigmaP_m(values);

    // initialize the covariance matrix to identity
    for (unsigned int i = 0; i < 6; ++i) {
        for (unsigned int j = 0; j < 6; ++j) {
            if (i == j)
                correlationMatrix_m[i][j] = 1.0;
            else
                correlationMatrix_m[i][j] = 0.0;
        }
    }

    cutoffR_m = ippl::Vector<double, 3>(
            std::abs(values.getReal("CUTOFFX")), std::abs(values.getReal("CUTOFFY")),
            std::abs(values.getReal("CUTOFFLONG")));

    correlationMatrix_m[1][0] = values.getReal("CORRX");
    correlationMatrix_m[3][2] = values.getReal("CORRY");
    correlationMatrix_m[5][4] = values.getReal("CORRT");

    if (values.getReal("CORRZ") != 0.0) correlationMatrix_m[5][4] = values.getReal("CORRZ");

    emitting_m = values.getBool("EMITTED");

    if (emitting_m) {
        sigmaR_m[2] = 0.0;

        sigmaTRise_m       = std::abs(values.getReal("SIGMAT"));
        sigmaTFall_m       = std::abs(values.getReal("SIGMAT"));
        tPulseLengthFWHM_m = std::abs(values.getReal("TPULSEFWHM"));

        FTOSCAmplitude_m = std::abs(values.getReal("FTOSCAMPLITUDE"));
        FTOSCPeriods_m   = std::abs(values.getReal("FTOSCPERIODS"));

        // If TRISE and TFALL are defined > 0.0 then these attributes
        // override SIGMAT.
        //
        if (std::abs(values.getReal("TRISE")) > 0.0 || std::abs(values.getReal("TFALL")) > 0.0) {
            double timeRatio =
                    std::sqrt(2.0 * std::log(10.0)) - std::sqrt(2.0 * std::log(10.0 / 9.0));

            if (std::abs(values.getReal("TRISE")) > 0.0)
                sigmaTRise_m = std::abs(values.getReal("TRISE")) / timeRatio;

            if (std::abs(values.getReal("TFALL")) > 0.0)
                sigmaTFall_m = std::abs(values.getReal("TFALL")) / timeRatio;
        }

        // For an emitted beam, the longitudinal cutoff >= 0.
        cutoffR_m[2] = std::abs(cutoffR_m[2]);
    }
}

void DistributionSettings::setSigmaR_m(const ObjectValues& values) {
    sigmaR_m = ippl::Vector<double, 3>(
            std::abs(values.getReal("SIGMAX")), std::abs(values.getReal("SIGMAY")),
            std::abs(values.getReal("SIGMAZ")));
}

void DistributionSettings::setSigmaP_m(const ObjectValues& values) {
    sigmaP_m = ippl::Vector<double, 3>(
            std::abs(values.getReal("SIGMAPX")), std::abs(values.getReal("SIGMAPY")),
            std::abs(values.getReal("SIGMAPZ")));
}
