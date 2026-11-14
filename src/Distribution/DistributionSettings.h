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
#ifndef OPALX_DISTRIBUTION_SETTINGS_H
#define OPALX_DISTRIBUTION_SETTINGS_H

#include <cstddef>
#include <string>

#include "Ippl.h"
#include "Track/ObjectValues.h"

enum class DistributionType : short {
    NODIST = -1,
    UNIFORM,  ///< Cold distribution uniform by volume inside an ellipsoid.
    GAUSS,
    MULTIVARIATEGAUSS,
    FLATTOP,
    OPALFLATTOP,
    FROMFILE,
    EMITTEDFROMFILE
};

/**
 * @brief The values of one DISTRIBUTION that the samplers read.
 *
 * Worked out from the attribute values of the DISTRIBUTION when a RUN builds its samplers, and
 * once when the DISTRIBUTION statement is read, to check them.
 */
class DistributionSettings {
public:
    using Matrix_t = ippl::Vector<ippl::Vector<double, 6>, 6>;

    /**
     * @brief Work out the sampler values of a DISTRIBUTION.
     * @param values The attribute values of the DISTRIBUTION.
     * @param avrgpz Mean longitudinal momentum of the beam (beta*gamma).
     * @throws OpalException if TYPE is missing, UNIFORM has a SIGMAX, SIGMAY or SIGMAZ that is
     *         not positive, GAUSS has a zero CUTOFFLONG, or CORR does not have 15 values.
     */
    DistributionSettings(const ObjectValues& values, double avrgpz);

    /// NPARTDIST of the DISTRIBUTION; 0 if it does not set its own number of particles.
    static std::size_t numParticles(const ObjectValues& values);

    std::string getTypeofDistribution() const { return distT_m; }
    DistributionType getType() const { return distrTypeT_m; }

    /// Return the number of macroparticles assigned to this distribution.
    /// If 0, the distribution did not specify its own NPARTDIST.
    std::size_t getNumParticles() const { return totalNumberParticles_m; }

    ippl::Vector<double, 3> getSigmaR() const { return sigmaR_m; }
    ippl::Vector<double, 3> getSigmaP() const { return sigmaP_m; }
    ippl::Vector<double, 3> getCutoffR() const { return cutoffR_m; }
    ippl::Vector<double, 3> getCutoffP() const { return cutoffP_m; }
    double getSigmaTRise() const { return sigmaTRise_m; }
    double getSigmaTFall() const { return sigmaTFall_m; }
    double getTPulseLengthFWHM() const { return tPulseLengthFWHM_m; }
    double getFTOSCAmplitude() const { return FTOSCAmplitude_m; }
    double getFTOSCPeriods() const { return FTOSCPeriods_m; }
    double getAvrgpz() const { return avrgpz_m; }
    /// EMISSIONSTEPS rounded up; 100 if not positive.
    std::size_t getEmissionSteps() const { return emissionSteps_m; }
    std::string getFilename() const { return filename_m; }

    Matrix_t correlationMatrix_m;
    /// True if particles are emitted over time; false if injected at initialization.
    bool emitting_m = false;

private:
    void setDistParametersUniform(const ObjectValues& values);
    void setDistParametersGauss(const ObjectValues& values);
    void setDistParametersMultiVariateGauss(const ObjectValues& values);
    void setDistParametersFlatTop(const ObjectValues& values);
    void setSigmaR_m(const ObjectValues& values);
    void setSigmaP_m(const ObjectValues& values);

    std::string distT_m;  ///< TYPE as given.
    DistributionType distrTypeT_m      = DistributionType::NODIST;
    std::size_t totalNumberParticles_m = 0;
    ippl::Vector<double, 3> sigmaR_m = 0.0, sigmaP_m = 0.0, cutoffR_m = 0.0, cutoffP_m = 0.0;
    double sigmaTRise_m         = 0;
    double sigmaTFall_m         = 0;
    double tPulseLengthFWHM_m   = 0;
    double avrgpz_m             = 0;
    double FTOSCAmplitude_m     = 0;
    double FTOSCPeriods_m       = 0;
    std::size_t emissionSteps_m = 100;
    std::string filename_m;
};

#endif  // OPALX_DISTRIBUTION_SETTINGS_H
