
// Class DistributionMoments
//   Computes the statistics of particle distributions.
//
// Copyright (c) 2025, Paul Scherrer Institut, Villigen PSI, Switzerland
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

#include "Algorithms/DistributionMoments.h"

#include <limits>

#include "AbstractObjects/OpalParticle.h"
#include "Utilities/GSLHistogram.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

const double DistributionMoments::percentileOneSigmaNormalDist_m    = std::erf(1 / sqrt(2));
const double DistributionMoments::percentileTwoSigmasNormalDist_m   = std::erf(2 / sqrt(2));
const double DistributionMoments::percentileThreeSigmasNormalDist_m = std::erf(3 / sqrt(2));
const double DistributionMoments::percentileFourSigmasNormalDist_m  = std::erf(4 / sqrt(2));

DistributionMoments::DistributionMoments() {
    reset();
    resetPlasmaParameters();
    moments_m        = matrix6x6_t(0.0);
    notCentMoments_m = matrix6x6_t(0.0);
}

// ---------------------------------------------------------------------------
// computeMeans -- single Kokkos reduce for 6 centroids + Ekin + gamma
// ---------------------------------------------------------------------------
void DistributionMoments::computeMeans(
        ippl::ParticleAttrib<Vector_t<double, 3>>::view_type Rview,
        ippl::ParticleAttrib<Vector_t<double, 3>>::view_type Pview,
        ippl::ParticleAttrib<double>::view_type Mview, size_t Np, size_t Nlocal) {
    /*
     Np is the total number of particles (reduced over ranks). In this function, it is only used to
     average over the number of total particles. For an empty simulation, this leads to divisions by
     0. Since, however, some of the computed moments might be needed regardless, we compute them but
     need to make sure that we do not divide by 0.
     Solution: Set Np to 1 if it is 0.
     */
    Np = (Np == 0) ? 1 : Np;

    const bool rescaleToReference = rescaleToReference_m;
    const double referenceMassGeV = referenceMassGeV_m;

    /*
    If the storage mode is SingleValue, then the mass is a scalar value that is the same for all
    particles. If the storage mode is Attributes, then the mass is a per-particle attribute.
    */
    const bool massIsScalarView = (Mview.extent(0) == 1);

    // Reduce 8 quantities in a single kernel:
    //   [0..5] centroids  (x, px, y, py, z, pz)
    //   [6]    kinetic energy sum
    //   [7]    gamma sum
    SumArray<8> loc;
    Kokkos::parallel_reduce(
            "computeMeans", Nlocal,
            KOKKOS_LAMBDA(const size_t k, SumArray<8>& upd) {
                upd.data[0] += Rview(k)[0];
                upd.data[1] += Pview(k)[0];
                upd.data[2] += Rview(k)[1];
                upd.data[3] += Pview(k)[1];
                upd.data[4] += Rview(k)[2];
                upd.data[5] += Pview(k)[2];

                double gamma0 = Kokkos::sqrt(dot(Pview(k)) + 1.0);

                const double massGeV = rescaleToReference
                                               ? referenceMassGeV
                                               : (massIsScalarView ? Mview(0) : Mview(k));
                upd.data[6] += (gamma0 - 1.0) * massGeV * Units::GeV2MeV;
                upd.data[7] += gamma0;
            },
            Kokkos::Sum<SumArray<8>>(loc));
    Kokkos::fence();

    double glob[8];
    ippl::Comm->allreduce(&loc.data[0], &glob[0], 8, std::plus<double>());

    for (size_t i = 0; i < 2 * Dim; ++i) {
        centroid_m(i) = glob[i];
        means_m(i)    = centroid_m[i] / Np;
    }

    meanKineticEnergy_m = glob[6] / (1. * Np);
    meanGamma_m         = glob[7] / (1. * Np);

    // gammaz reuses the pz centroid (index 5): gammaz = sqrt( (<pz>/Np)^2 + 1 )
    double gammaz = glob[5] / (1. * Np);
    gammaz        = std::sqrt(gammaz * gammaz + 1.0);
    meanGammaZ_m  = gammaz;

    for (size_t i = 0; i < Dim; ++i) {
        meanR_m(i) = means_m[2 * i];
        meanP_m(i) = means_m[2 * i + 1];
    }
}

// ---------------------------------------------------------------------------
// computeMoments -- three Kokkos reducers in a single parallel_reduce:
//   SumMatrix6x6  central moments,  SumMatrix6x6  non-central moments,
//   double        std kinetic energy
// ---------------------------------------------------------------------------
void DistributionMoments::computeMoments(
        ippl::ParticleAttrib<Vector_t<double, 3>>::view_type Rview,
        ippl::ParticleAttrib<Vector_t<double, 3>>::view_type Pview,
        ippl::ParticleAttrib<double>::view_type Mview, size_t Np, size_t Nlocal) {
    // Check that the BunchStateHandler is set. Should always be set, so throw an
    // exception if not.
    if (!bunchStateHandler_m) {
        throw OpalException(
                "DistributionMoments::computeMoments",
                "BunchStateHandler not set, cannot use "
                "DistributionMoments instance correctly.");
    }

    Np = (Np == 0) ? 1 : Np;  // Explanation: see DistributionMoments::computeMeans
                              // implementation

    if (!bunchStateHandler_m->isMomentsDirty()) {
        return;
    }

    IpplTimings::TimerRef momentsTimer = IpplTimings::getTimer("computeMoments");
    IpplTimings::startTimer(momentsTimer);

    reset();
    computeMeans(Rview, Pview, Mview, Np, Nlocal);

    double meanR_loc[Dim];
    double meanP_loc[Dim];
    for (size_t i = 0; i < Dim; ++i) {
        meanR_loc[i] = meanR_m[i];
        meanP_loc[i] = meanP_m[i];
    }

    const double mekin               = meanKineticEnergy_m;
    const bool rescaleToReferenceStd = rescaleToReference_m;
    const double referenceMassGeVStd = referenceMassGeV_m;
    const bool massIsScalarView      = (Mview.extent(0) == 1);

    SumMatrix6x6 loc_cent;
    SumMatrix6x6 loc_ncent;
    double loc_std_ekin = 0.0;

    Kokkos::parallel_reduce(
            "computeMoments", Nlocal,
            KOKKOS_LAMBDA(
                    const size_t k, SumMatrix6x6& cent_upd, SumMatrix6x6& ncent_upd,
                    double& ekin_upd) {
                double cent[6];
                cent[0] = Rview(k)[0] - meanR_loc[0];
                cent[1] = Pview(k)[0] - meanP_loc[0];
                cent[2] = Rview(k)[1] - meanR_loc[1];
                cent[3] = Pview(k)[1] - meanP_loc[1];
                cent[4] = Rview(k)[2] - meanR_loc[2];
                cent[5] = Pview(k)[2] - meanP_loc[2];

                double raw[6];
                raw[0] = Rview(k)[0];
                raw[1] = Pview(k)[0];
                raw[2] = Rview(k)[1];
                raw[3] = Pview(k)[1];
                raw[4] = Rview(k)[2];
                raw[5] = Pview(k)[2];

                for (size_t i = 0; i < 6; ++i) {
                    for (size_t j = 0; j < 6; ++j) {
                        cent_upd.data[i][j] += cent[i] * cent[j];
                        ncent_upd.data[i][j] += raw[i] * raw[j];
                    }
                }

                double gamma0        = Kokkos::sqrt(dot(Pview(k)) + 1.0);
                const double massGeV = rescaleToReferenceStd
                                               ? referenceMassGeVStd
                                               : (massIsScalarView ? Mview(0) : Mview(k));
                double ekin0         = (gamma0 - 1.0) * massGeV * Units::GeV2MeV;
                ekin_upd += (ekin0 - mekin) * (ekin0 - mekin);
            },
            Kokkos::Sum<SumMatrix6x6>(loc_cent), Kokkos::Sum<SumMatrix6x6>(loc_ncent),
            Kokkos::Sum<double>(loc_std_ekin));

    SumMatrix6x6 glob_cent, glob_ncent;
    ippl::Comm->allreduce(&loc_cent.data[0][0], &glob_cent.data[0][0], 36, std::plus<double>());
    ippl::Comm->allreduce(&loc_ncent.data[0][0], &glob_ncent.data[0][0], 36, std::plus<double>());

    double glob_std_ekin;
    ippl::Comm->allreduce(&loc_std_ekin, &glob_std_ekin, 1, std::plus<double>());

    for (size_t i = 0; i < 2 * Dim; ++i) {
        for (size_t j = 0; j < 2 * Dim; ++j) {
            moments_m(i, j)        = glob_cent.data[i][j] / Np;
            notCentMoments_m(i, j) = glob_ncent.data[i][j];
        }
    }

    for (size_t i = 0; i < Dim; ++i) {
        stdR_m(i) = std::sqrt(moments_m(2 * i, 2 * i));
        stdP_m(i) = std::sqrt(moments_m(2 * i + 1, 2 * i + 1));
    }

    stdKineticEnergy_m = std::sqrt(glob_std_ekin / Np);

    double perParticle = 1. / (1. * Np);
    Vector_t<double, 3> squaredEps, sumRP;
    size_t l = 0;
    for (size_t i = 0; i < 3; ++i, l += 2) {
        double w1  = centroid_m[2 * i] * perParticle;
        double w2  = moments_m(2 * i, 2 * i);
        double w3  = notCentMoments_m(l, l) * perParticle;
        double w4  = notCentMoments_m(l + 1, l + 1) * perParticle;
        double tmp = w2 - std::pow(w1, 2);

        halo_m(i) = (w4 + w1 * (-4 * w3 + 3 * w1 * (tmp + w2))) / tmp;
        halo_m(i) -= Options::haloShift;
    }

    for (size_t i = 0; i < 3; ++i) {
        sumRP(i)      = notCentMoments_m(2 * i, 2 * i + 1) * perParticle - meanR_m(i) * meanP_m(i);
        stdRP_m(i)    = sumRP(i) / (stdR_m(i) * stdP_m(i));
        squaredEps(i) = std::pow(stdR_m(i) * stdP_m(i), 2) - std::pow(sumRP(i), 2);
        normalizedEps_m(i) = std::sqrt(std::max(squaredEps(i), 0.0));
    }

    double betaGamma = std::sqrt(std::pow(meanGamma_m, 2) - 1.0);
    geometricEps_m   = normalizedEps_m / Vector_t<double, 3>(betaGamma);

    IpplTimings::stopTimer(momentsTimer);
}

// ---------------------------------------------------------------------------
// computeMinMaxPosition -- MaxArray<3> + MinArray<3> in a single reduce
// ---------------------------------------------------------------------------
void DistributionMoments::computeMinMaxPosition(
        ippl::ParticleAttrib<Vector_t<double, 3>>::view_type Rview, size_t Nlocal) {
    // Identity values (lowest / max) are set by the default constructors, so
    // ranks with Nlocal == 0 naturally contribute neutral elements to the allreduce.
    Inform m("DistributionMoments::computeMinMaxPosition");
    MaxArray<3> loc_max;
    MinArray<3> loc_min;

    Kokkos::parallel_reduce(
            "computeMinMaxPosition", Nlocal,
            KOKKOS_LAMBDA(const size_t k, MaxArray<3>& rmax, MinArray<3>& rmin) {
                for (size_t d = 0; d < 3; ++d) {
                    const double r = Rview(k)[d];
                    if (r > rmax.data[d]) rmax.data[d] = r;
                    if (r < rmin.data[d]) rmin.data[d] = r;
                }
            },
            Kokkos::Sum<MaxArray<3>>(loc_max), Kokkos::Sum<MinArray<3>>(loc_min));

    double rmax[Dim], rmin[Dim];
    ippl::Comm->allreduce(&loc_max.data[0], &rmax[0], Dim, std::greater<double>());
    ippl::Comm->allreduce(&loc_min.data[0], &rmin[0], Dim, std::less<double>());

    // Empty bunch on all ranks: min stays id_min, max stays id_max -> min > max.
    // Use a tiny valid box so mesh/spacing remain valid (e.g. for first step before emission). The
    // same problem arises when one dimension is 0 (e.g. flattop first timestep), then rmin>=rmax
    // would trigger "=0". Therefore, we HAVE to check rmin>rmax instead.
    for (size_t i = 0; i < Dim; ++i) {
        if (rmin[i] > rmax[i]) {
            rmin[i] = rmax[i] = 0.0;
            m << level5 << "Min > max in dimension " << i << ". Setting to 0 (degenerate case)."
              << endl;
        }
    }

    for (size_t i = 0; i < Dim; ++i) {
        minR_m(i) = rmin[i];
        maxR_m(i) = rmax[i];
    }
    m << level5 << "Min/max computation done." << endl;
}

// ---------------------------------------------------------------------------

void DistributionMoments::compute(
        const std::vector<OpalParticle>::const_iterator& /*first*/,
        const std::vector<OpalParticle>::const_iterator& /*last*/) {
    throw OpalException(
            "DistributionMoments::compute", "OpalParticle-iterator interface is not implemented.");
}

template <class InputIt>
void DistributionMoments::computePercentiles(const InputIt& first, const InputIt& last) {
    if (!Options::computePercentiles || totalNumParticles_m < 100) {
        return;
    }

    std::vector<gsl_histogram*> histograms(3);
    // For a normal distribution the number of exchanged data between the cores is minimized
    // if the number of histogram bins follows the following formula. Since we can't know
    // how many particles are in each bin for the real distribution we use this formula.
    unsigned int numBins = 3.5 * std::pow(3, std::log10(totalNumParticles_m));

    Vector_t<double, 3> maxR;
    for (size_t d = 0; d < 3; ++d) {
        maxR(d)       = 1.0000001 * std::max(maxR_m[d] - meanR_m[d], meanR_m[d] - minR_m[d]);
        histograms[d] = gsl_histogram_alloc(numBins);
        gsl_histogram_set_ranges_uniform(histograms[d], 0.0, maxR(d));
    }
    for (InputIt it = first; it != last; ++it) {
        OpalParticle const& particle = *it;
        for (size_t d = 0; d < 3; ++d) {
            gsl_histogram_increment(histograms[d], std::abs(particle[2 * d] - meanR_m[d]));
        }
    }

    std::vector<int> localHistogramValues(3 * (numBins + 1)),
            globalHistogramValues(3 * (numBins + 1));
    for (size_t d = 0; d < 3; ++d) {
        int j              = 0;
        size_t accumulated = 0;
        std::generate(
                localHistogramValues.begin() + d * (numBins + 1) + 1,
                localHistogramValues.begin() + (d + 1) * (numBins + 1),
                [&histograms, &d, &j, &accumulated]() {
                    accumulated += gsl_histogram_get(histograms[d], j++);
                    return accumulated;
                });

        gsl_histogram_free(histograms[d]);
    }

    ippl::Comm->allreduce(
            localHistogramValues.data(), globalHistogramValues.data(), 3 * (numBins + 1),
            std::plus<int>());

    int numParticles68 = static_cast<int>(
            std::floor(totalNumParticles_m * percentileOneSigmaNormalDist_m + 0.5));
    int numParticles95 = static_cast<int>(
            std::floor(totalNumParticles_m * percentileTwoSigmasNormalDist_m + 0.5));
    int numParticles99 = static_cast<int>(
            std::floor(totalNumParticles_m * percentileThreeSigmasNormalDist_m + 0.5));
    int numParticles99_99 = static_cast<int>(
            std::floor(totalNumParticles_m * percentileFourSigmasNormalDist_m + 0.5));

    for (int d = 0; d < 3; ++d) {
        unsigned int localNum = last - first, current = 0;
        std::vector<Vector_t<double, 2>> oneDPhaseSpace(localNum);
        for (InputIt it = first; it != last; ++it, ++current) {
            OpalParticle const& particle = *it;
            oneDPhaseSpace[current](0)   = particle[2 * d];
            oneDPhaseSpace[current](1)   = particle[2 * d + 1];
        }
        std::sort(
                oneDPhaseSpace.begin(), oneDPhaseSpace.end(),
                [d, this](Vector_t<double, 2>& left, Vector_t<double, 2>& right) {
                    return std::abs(left[0] - meanR_m[d]) < std::abs(right[0] - meanR_m[d]);
                });

        iterator_t endSixtyEight, endNinetyFive, endNinetyNine, endNinetyNine_NinetyNine;
        endSixtyEight = endNinetyFive = endNinetyNine = endNinetyNine_NinetyNine =
                oneDPhaseSpace.end();

        std::tie(sixtyEightPercentile_m[d], endSixtyEight) = determinePercentilesDetail(
                oneDPhaseSpace.begin(), oneDPhaseSpace.end(), globalHistogramValues,
                localHistogramValues, d, numParticles68);

        std::tie(ninetyFivePercentile_m[d], endNinetyFive) = determinePercentilesDetail(
                oneDPhaseSpace.begin(), oneDPhaseSpace.end(), globalHistogramValues,
                localHistogramValues, d, numParticles95);

        std::tie(ninetyNinePercentile_m[d], endNinetyNine) = determinePercentilesDetail(
                oneDPhaseSpace.begin(), oneDPhaseSpace.end(), globalHistogramValues,
                localHistogramValues, d, numParticles99);

        std::tie(ninetyNine_NinetyNinePercentile_m[d], endNinetyNine_NinetyNine) =
                determinePercentilesDetail(
                        oneDPhaseSpace.begin(), oneDPhaseSpace.end(), globalHistogramValues,
                        localHistogramValues, d, numParticles99_99);

        normalizedEps68Percentile_m[d] =
                computeNormalizedEmittance(oneDPhaseSpace.begin(), endSixtyEight);
        normalizedEps95Percentile_m[d] =
                computeNormalizedEmittance(oneDPhaseSpace.begin(), endNinetyFive);
        normalizedEps99Percentile_m[d] =
                computeNormalizedEmittance(oneDPhaseSpace.begin(), endNinetyNine);
        normalizedEps99_99Percentile_m[d] =
                computeNormalizedEmittance(oneDPhaseSpace.begin(), endNinetyNine_NinetyNine);
    }
}

/** Computes the percentile and the range of all local particles that are contained therein.
 *  In a first step the container globalAccumulatedHistogram is looped through until accumulated
 *  histogram value exceeds the required number of particles. The percentile then is between the
 *  boundaries of the last histogram bin before the loop stopped. Then all particle coordinates
 *  that are between the boundaries of this bin are communicated acros all nodes and sorted.
 *  The exact percentile is then determined by counting the n smallest coordinates such that
 *  the total number of partiles results is equal to 'numRequiredParticles'. In accordance with
 *  matlab (?) the percentile is the midpoint between the last particle within the percentile and
 *  tje first particle outside. Finally each node determines which of its particles are contained
 *  in the percentile.
 *
 *  To determine the histogram, the coordinates should not be used directly. Instead, the
 *  absolute value of the difference between a coordinate and the mean, |x - <x>|, should be used
 *  so that the percentile values are similar to the standard deviation.

 * @param begin: begin of a container containing the one dimensional phase space of all local
 *               particles.
 * @param end:   end of the container.
 * @param globalAccumulatedHistogram: container with partial sum of histogram values of position
 *                                    coordinates summed up across all nodes. The first value
 *                                    should be 0.
 * @param localAccumulatedHistogram: container with partial sum of histogram values of position
 *                                   coordinates of all local particles. The first value should be
 0.
 * @param dimension: dimension of the one dimensional phase space.
 * @param numRequiredParticles: number of particles that are contained in the requested percentile.
 *                              Is determined by the total number of particles and the percentile.
 * @return: pair of percentile and iterator pointing to the element after the range contained in the
 *          percentile.
 */
std::pair<double, DistributionMoments::iterator_t> DistributionMoments::determinePercentilesDetail(
        const DistributionMoments::iterator_t& begin, const DistributionMoments::iterator_t& end,
        const std::vector<int>& globalAccumulatedHistogram,
        const std::vector<int>& localAccumulatedHistogram, unsigned int dimension,
        int numRequiredParticles) const {
    unsigned int numBins     = globalAccumulatedHistogram.size() / 3;
    double percentile        = 0.0;
    iterator_t endPercentile = end;
    for (unsigned int i = 1; i < numBins; ++i) {
        unsigned int idx = dimension * numBins + i;
        if (globalAccumulatedHistogram[idx] > numRequiredParticles) {
            iterator_t beginBin = begin + localAccumulatedHistogram[idx - 1];
            iterator_t endBin   = begin + localAccumulatedHistogram[idx];
            unsigned int numMissingParticles =
                    numRequiredParticles - globalAccumulatedHistogram[idx - 1];
            unsigned int shift = 2;
            while (numMissingParticles == 0) {
                beginBin = begin + localAccumulatedHistogram[idx - shift];
                numMissingParticles =
                        numRequiredParticles - globalAccumulatedHistogram[idx - shift];
                ++shift;
            }

            std::vector<unsigned int> numParticlesInBin(ippl::Comm->size() + 1);
            numParticlesInBin[ippl::Comm->rank() + 1] = endBin - beginBin;

            ippl::Comm->allreduce(
                    &(numParticlesInBin[1]), ippl::Comm->size(), std::plus<unsigned int>());

            std::partial_sum(
                    numParticlesInBin.begin(), numParticlesInBin.end(), numParticlesInBin.begin());

            std::vector<double> positions(numParticlesInBin.back());
            std::transform(
                    beginBin, endBin, positions.begin() + numParticlesInBin[ippl::Comm->rank()],
                    [&dimension, this](Vector_t<double, 2> const& particle) {
                        return std::abs(particle[0] - meanR_m[dimension]);
                    });
            ippl::Comm->allreduce(&(positions[0]), positions.size(), std::plus<double>());
            std::sort(positions.begin(), positions.end());

            percentile = (*(positions.begin() + numMissingParticles - 1)
                          + *(positions.begin() + numMissingParticles))
                         / 2;
            for (iterator_t it = beginBin; it != endBin; ++it) {
                if (std::abs((*it)[0] - meanR_m[dimension]) > percentile) {
                    return std::make_pair(percentile, it);
                }
            }
            endPercentile = endBin;
            break;
        }
    }
    return std::make_pair(percentile, endPercentile);
}

double DistributionMoments::computeNormalizedEmittance(
        const DistributionMoments::iterator_t& begin,
        const DistributionMoments::iterator_t& end) const {
    double localStatistics[] = {0.0, 0.0, 0.0, 0.0};
    localStatistics[0]       = end - begin;
    for (iterator_t it = begin; it < end; ++it) {
        const Vector_t<double, 2>& rp = *it;
        localStatistics[1] += rp(0);
        localStatistics[2] += rp(1);
        localStatistics[3] += rp(0) * rp(1);
    }
    ippl::Comm->allreduce(&(localStatistics[0]), 4, std::plus<double>());

    double numParticles = localStatistics[0];
    double perParticle  = 1 / localStatistics[0];
    double meanR        = localStatistics[1] * perParticle;
    double meanP        = localStatistics[2] * perParticle;
    double RP           = localStatistics[3] * perParticle;

    localStatistics[0] = 0.0;
    localStatistics[1] = 0.0;
    for (iterator_t it = begin; it < end; ++it) {
        const Vector_t<double, 2>& rp = *it;
        localStatistics[0] += std::pow(rp(0) - meanR, 2);
        localStatistics[1] += std::pow(rp(1) - meanP, 2);
    }
    ippl::Comm->allreduce(&(localStatistics[0]), 2, std::plus<double>());

    double stdR          = std::sqrt(localStatistics[0] / numParticles);
    double stdP          = std::sqrt(localStatistics[1] / numParticles);
    double sumRP         = RP - meanR * meanP;
    double squaredEps    = std::pow(stdR * stdP, 2) - std::pow(sumRP, 2);
    double normalizedEps = std::sqrt(std::max(squaredEps, 0.0));

    return normalizedEps;
}

// ---------------------------------------------------------------------------

void DistributionMoments::fillMembers(std::vector<double>& /*localMoments*/) {
    throw OpalException(
            "DistributionMoments::fillMembers", "Not implemented -- must not be called.");
}

void DistributionMoments::computeMeanKineticEnergy() {
    throw OpalException(
            "DistributionMoments::computeMeanKineticEnergy",
            "Not implemented -- must not be called.");
}

// ---------------------------------------------------------------------------
// computeDebyeLength -- SumArray<3> for avg velocity, scalar for temperature
// ---------------------------------------------------------------------------
void DistributionMoments::computeDebyeLength(
        ippl::ParticleAttrib<Vector_t<double, 3>>::view_type Pview, size_t Np, size_t Nlocal,
        double density) {
    Np = (Np == 0) ? 1 : Np;

    resetPlasmaParameters();
    const double c = Physics::c;

    // v = (P * c) / gamma  -->  average velocity per component
    SumArray<3> loc_vel;
    Kokkos::parallel_reduce(
            "computeDebyeLength avgVel", Nlocal,
            KOKKOS_LAMBDA(const size_t k, SumArray<3>& upd) {
                double gamma0 = Kokkos::sqrt(dot(Pview(k)) + 1.0);

                upd.data[0] += Pview(k)[0] * c / gamma0;
                upd.data[1] += Pview(k)[1] * c / gamma0;
                upd.data[2] += Pview(k)[2] * c / gamma0;
            },
            Kokkos::Sum<SumArray<3>>(loc_vel));
    Kokkos::fence();

    double avgVel[3];
    ippl::Comm->allreduce(&loc_vel.data[0], &avgVel[0], Dim, std::plus<double>());

    for (unsigned i = 0; i < 3; ++i)
        avgVel[i] /= Np;

    double loc_tempAvg = 0.0;
    Kokkos::parallel_reduce(
            "computeDebyeLength temperature", Nlocal,
            KOKKOS_LAMBDA(const size_t k, double& mom0) {
                double gamma0 = Kokkos::sqrt(dot(Pview(k)) + 1.0);

                mom0 += Kokkos::pow(Pview(k)[0] * c / gamma0 - avgVel[0], 2);
                mom0 += Kokkos::pow(Pview(k)[1] * c / gamma0 - avgVel[1], 2);
                mom0 += Kokkos::pow(Pview(k)[2] * c / gamma0 - avgVel[2], 2);
            },
            Kokkos::Sum<double>(loc_tempAvg));
    Kokkos::fence();

    double tempAvg;
    ippl::Comm->allreduce(&loc_tempAvg, &tempAvg, 1, std::plus<double>());

    // k_B T in units of kg m^2/s^2
    temperature_m = (1.0 / 3.0) * Units::eV2kg * Units::GeV2eV * Physics::m_e * (tempAvg / Np);

    debyeLength_m =
            std::sqrt((temperature_m * Physics::epsilon_0) / (density * std::pow(Physics::q_e, 2)));

    computePlasmaParameter(density);
}

void DistributionMoments::computePlasmaParameter(double density) {
    plasmaParameter_m = (4.0 / 3.0) * Physics::pi * std::pow(debyeLength_m, 3) * density;
}

// ---------------------------------------------------------------------------

void DistributionMoments::reset() {
    std::fill(std::begin(centroid_m), std::end(centroid_m), 0.0);
    meanR_m         = 0.0;
    meanP_m         = 0.0;
    stdR_m          = 0.0;
    stdP_m          = 0.0;
    stdRP_m         = 0.0;
    normalizedEps_m = 0.0;
    geometricEps_m  = 0.0;
    halo_m          = 0.0;

    meanKineticEnergy_m = 0.0;
    stdKineticEnergy_m  = 0.0;

    totalCharge_m       = 0.0;
    totalMass_m         = 0.0;
    totalNumParticles_m = 0;

    sixtyEightPercentile_m            = 0.0;
    normalizedEps68Percentile_m       = 0.0;
    ninetyFivePercentile_m            = 0.0;
    normalizedEps95Percentile_m       = 0.0;
    ninetyNinePercentile_m            = 0.0;
    normalizedEps99Percentile_m       = 0.0;
    ninetyNine_NinetyNinePercentile_m = 0.0;
    normalizedEps99_99Percentile_m    = 0.0;
}

void DistributionMoments::resetPlasmaParameters() {
    temperature_m     = 0.0;
    debyeLength_m     = 0.0;
    plasmaParameter_m = 0.0;
}

bool DistributionMoments::isParticleExcluded(const OpalParticle& /*particle*/) const {
    throw OpalException(
            "DistributionMoments::isParticleExcluded", "Not implemented -- must not be called.");
}
