//
// Class DistributionMoments
//   Computes the statistics of particle distributions.
//
// Copyright (c) 2021, Christof Metzger-Kraus
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
#ifndef DISTRIBUTIONMOMENTS_H
#define DISTRIBUTIONMOMENTS_H

#include "AbsBeamline/Component.h"
#include "PartBunch/PartBunch.hpp"

// #include "Matrix.h"

#include "Physics/Physics.h"
#include "Physics/Units.h"

#include <vector>

class OpalParticle;

class DistributionMoments {
public:
    DistributionMoments();

    void compute(
        const std::vector<OpalParticle>::const_iterator&,
        const std::vector<OpalParticle>::const_iterator&);
    void compute(PartBunch_t const&);
    void computeMeanKineticEnergy(PartBunch_t const&);
    void computeDebyeLength(PartBunch_t const&, double density);
    void computePlasmaParameter(double);

    Vector_t<double, 3> getMeanPosition() const;
    Vector_t<double, 3> getStandardDeviationPosition() const;
    Vector_t<double, 3> getMeanMomentum() const;
    Vector_t<double, 3> getStandardDeviationMomentum() const;
    Vector_t<double, 3> getNormalizedEmittance() const;
    Vector_t<double, 3> getGeometricEmittance() const;
    Vector_t<double, 3> getStandardDeviationRP() const;
    Vector_t<double, 3> getHalo() const;
    Vector_t<double, 3> getMaxR() const;

    Vector_t<double, 3> get68Percentile() const;
    Vector_t<double, 3> get95Percentile() const;
    Vector_t<double, 3> get99Percentile() const;
    Vector_t<double, 3> get99_99Percentile() const;
    Vector_t<double, 3> getNormalizedEmittance68Percentile() const;
    Vector_t<double, 3> getNormalizedEmittance95Percentile() const;
    Vector_t<double, 3> getNormalizedEmittance99Percentile() const;
    Vector_t<double, 3> getNormalizedEmittance99_99Percentile() const;

    double getMeanTime() const;
    double getStdTime() const;
    double getMeanGamma() const;
    double getMeanKineticEnergy() const;
    double getTemperature() const;
    double getDebyeLength() const;
    double getPlasmaParameter() const;
    double getStdKineticEnergy() const;
    double getDx() const;
    double getDDx() const;
    double getDy() const;
    double getDDy() const;
    matrix_t getMoments6x6() const;
    double getTotalCharge() const;
    double getTotalMass() const;
    double getTotalNumParticles() const;

private:
    bool isParticleExcluded(const OpalParticle&) const;
    template <class InputIt>
    void computeMeans(const InputIt&, const InputIt&);
    // template <class InputIt>
    // void computeStatistics(const InputIt&, const InputIt&);
    template <class InputIt>
    void computePercentiles(const InputIt&, const InputIt&);
    using iterator_t = std::vector<Vector_t<double, 2>>::const_iterator;
    std::pair<double, iterator_t> determinePercentilesDetail(
        const iterator_t& begin, const iterator_t& end,
        const std::vector<int>& globalAccumulatedHistogram,
        const std::vector<int>& localAccumulatedHistogram, unsigned int dimension,
        int numRequiredParticles) const;
    double computeNormalizedEmittance(const iterator_t& begin, const iterator_t& end) const;

    void fillMembers(std::vector<double>&);

    void reset();

    void resetPlasmaParameters();

    Vector_t<double, 3> meanR_m;
    Vector_t<double, 3> meanP_m;
    Vector_t<double, 3> stdR_m;
    Vector_t<double, 3> stdP_m;
    Vector_t<double, 3> stdRP_m;
    Vector_t<double, 3> normalizedEps_m;
    Vector_t<double, 3> geometricEps_m;
    Vector_t<double, 3> halo_m;
    Vector_t<double, 3> maxR_m;
    Vector_t<double, 3> minR_m;
    Vector_t<double, 3> sixtyEightPercentile_m;
    Vector_t<double, 3> normalizedEps68Percentile_m;
    Vector_t<double, 3> ninetyFivePercentile_m;
    Vector_t<double, 3> normalizedEps95Percentile_m;
    Vector_t<double, 3> ninetyNinePercentile_m;
    Vector_t<double, 3> normalizedEps99Percentile_m;
    Vector_t<double, 3> ninetyNine_NinetyNinePercentile_m;
    Vector_t<double, 3> normalizedEps99_99Percentile_m;

    double meanTime_m;
    double stdTime_m;
    double meanKineticEnergy_m;
    double temperature_m;
    double debyeLength_m;
    double plasmaParameter_m;
    double stdKineticEnergy_m;
    double meanGamma_m;

    double centroid_m[6];
    matrix_t moments_m;

    double totalCharge_m;
    double totalMass_m;
    unsigned int totalNumParticles_m;

    static const double percentileOneSigmaNormalDist_m;
    static const double percentileTwoSigmasNormalDist_m;
    static const double percentileThreeSigmasNormalDist_m;
    static const double percentileFourSigmasNormalDist_m;
};

inline Vector_t<double, 3> DistributionMoments::getMeanPosition() const {
    return meanR_m;
}

inline Vector_t<double, 3> DistributionMoments::getStandardDeviationPosition() const {
    return stdR_m;
}

inline Vector_t<double, 3> DistributionMoments::getMeanMomentum() const {
    return meanP_m;
}

inline Vector_t<double, 3> DistributionMoments::getStandardDeviationMomentum() const {
    return stdP_m;
}

inline Vector_t<double, 3> DistributionMoments::getNormalizedEmittance() const {
    return normalizedEps_m;
}

inline Vector_t<double, 3> DistributionMoments::getGeometricEmittance() const {
    return geometricEps_m;
}

inline Vector_t<double, 3> DistributionMoments::getStandardDeviationRP() const {
    return stdRP_m;
}

inline Vector_t<double, 3> DistributionMoments::getHalo() const {
    return halo_m;
}

inline double DistributionMoments::getMeanTime() const {
    return meanTime_m;
}

inline double DistributionMoments::getStdTime() const {
    return stdTime_m;
}

inline double DistributionMoments::getMeanGamma() const {
    return meanGamma_m;
}

inline double DistributionMoments::getMeanKineticEnergy() const {
    return meanKineticEnergy_m;
}

// Compute and return the value of temperature in K
inline double DistributionMoments::getTemperature() const {
    return (temperature_m / (Physics::kB * Units::eV2kg * Physics::c * Physics::c));
}
inline double DistributionMoments::getDebyeLength() const {
    return debyeLength_m;
}
inline double DistributionMoments::getPlasmaParameter() const {
    return plasmaParameter_m;
}

inline double DistributionMoments::getStdKineticEnergy() const {
    return stdKineticEnergy_m;
}

inline double DistributionMoments::getDx() const {
    return moments_m(0, 5);
}

inline double DistributionMoments::getDDx() const {
    return moments_m(1, 5);
}

inline double DistributionMoments::getDy() const {
    return moments_m(2, 5);
}

inline double DistributionMoments::getDDy() const {
    return moments_m(3, 5);
}

inline matrix_t DistributionMoments::getMoments6x6() const {
    return moments_m;
}

inline double DistributionMoments::getTotalCharge() const {
    return totalCharge_m;
}

inline double DistributionMoments::getTotalMass() const {
    return totalMass_m;
}

inline double DistributionMoments::getTotalNumParticles() const {
    return totalNumParticles_m;
}

inline Vector_t<double, 3> DistributionMoments::get68Percentile() const {
    return sixtyEightPercentile_m;
}

inline Vector_t<double, 3> DistributionMoments::getNormalizedEmittance68Percentile() const {
    return normalizedEps68Percentile_m;
}

inline Vector_t<double, 3> DistributionMoments::get95Percentile() const {
    return ninetyFivePercentile_m;
}

inline Vector_t<double, 3> DistributionMoments::getNormalizedEmittance95Percentile() const {
    return normalizedEps95Percentile_m;
}

inline Vector_t<double, 3> DistributionMoments::get99Percentile() const {
    return ninetyNinePercentile_m;
}

inline Vector_t<double, 3> DistributionMoments::getNormalizedEmittance99Percentile() const {
    return normalizedEps99Percentile_m;
}

inline Vector_t<double, 3> DistributionMoments::get99_99Percentile() const {
    return ninetyNine_NinetyNinePercentile_m;
}

inline Vector_t<double, 3> DistributionMoments::getNormalizedEmittance99_99Percentile() const {
    return normalizedEps99_99Percentile_m;
}

inline Vector_t<double, 3> DistributionMoments::getMaxR() const {
    Vector_t<double, 3> maxDistance;
    for (unsigned int i = 0; i < 3; ++i) {
        maxDistance[i] = std::max(std::abs(maxR_m[i]), std::abs(minR_m[i]));
    }
    return maxDistance;
}

#endif
