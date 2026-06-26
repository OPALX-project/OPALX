#ifndef OPALX_TEST_LINEAR_BREIT_WHEELER_BENCHMARK_COMMON_H
#define OPALX_TEST_LINEAR_BREIT_WHEELER_BENCHMARK_COMMON_H

#include "Physics/LinearBreitWheeler.h"
#include "Physics/Physics.h"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace LinearBreitWheelerBenchmark {

    enum class FinalState {
        Electron,
        Positron,
    };

    enum class Observable {
        Energy,
        Theta,
    };

    /**
     * @brief Configuration for folding the linear Breit-Wheeler kernel over a finite incoming
     * photon beam.
     *
     * This benchmark keeps the laser photon fixed and samples only the incoming
     * high-energy photon beam. The present implementation models the beam by a
     * Gaussian spread in the transverse photon slopes around the reference beam
     * axis and, optionally, by a Gaussian relative energy spread.
     *
     * The benchmark intentionally stays momentum-space only. It does not model
     * photon-beam position spread or laser-overlap weighting yet.
     */
    struct FinitePhotonBeamConfig {
        double centralHighEnergyPhotonEnergyGeV          = 0.5;
        double wavelength_m                              = 1.0e-9;
        Vector_t<double, 3> referenceHighEnergyDirection = [] {
            Vector_t<double, 3> value(0.0);
            value(2) = 1.0;
            return value;
        }();
        Vector_t<double, 3> laserDirection = [] {
            Vector_t<double, 3> value(0.0);
            value(2) = -1.0;
            return value;
        }();
        double sigmaThetaXRad       = 1.0e-3;
        double sigmaThetaYRad       = 1.0e-3;
        double relativeEnergySpread = 0.0;
        double sigmaX_m             = 0.0;
        double sigmaY_m             = 0.0;
        double sigmaS_m             = 0.0;
        double laserRayleighX_m     = 1.0e-6;
        double laserRayleighY_m     = 1.0e-6;
        double laserSigmaT_m        = 5.0e-12 * Physics::c;
        bool overlapWeighting       = false;
        std::size_t bins            = 80;
        double minValue             = 0.0;
        double maxValue             = 0.5;
    };

    struct HistogramConfig {
        double highEnergyPhotonEnergyGeV        = 0.5;
        double wavelength_m                     = 1.0e-9;
        Vector_t<double, 3> highEnergyDirection = [] {
            Vector_t<double, 3> value(0.0);
            value(2) = 1.0;
            return value;
        }();
        Vector_t<double, 3> laserDirection = [] {
            Vector_t<double, 3> value(0.0);
            value(2) = -1.0;
            return value;
        }();
        std::size_t bins = 80;
        double minValue  = 0.0;
        double maxValue  = 0.5;
    };

    struct Histogram {
        std::vector<double> centers;
        std::vector<double> density;
        std::vector<double> counts;
        double binWidth    = 0.0;
        double totalWeight = 0.0;
    };

    struct JointHistogramConfig {
        double highEnergyPhotonEnergyGeV        = 0.5;
        double wavelength_m                     = 1.0e-9;
        Vector_t<double, 3> highEnergyDirection = [] {
            Vector_t<double, 3> value(0.0);
            value(2) = 1.0;
            return value;
        }();
        Vector_t<double, 3> laserDirection = [] {
            Vector_t<double, 3> value(0.0);
            value(2) = -1.0;
            return value;
        }();
        std::size_t energyBins = 80;
        std::size_t thetaBins  = 80;
        double energyMinGeV    = 0.0;
        double energyMaxGeV    = 0.5;
        double thetaMinRad     = 0.0;
        double thetaMaxRad     = 0.0045;
    };

    struct FinitePhotonBeamJointConfig {
        double centralHighEnergyPhotonEnergyGeV          = 0.5;
        double wavelength_m                              = 1.0e-9;
        Vector_t<double, 3> referenceHighEnergyDirection = [] {
            Vector_t<double, 3> value(0.0);
            value(2) = 1.0;
            return value;
        }();
        Vector_t<double, 3> laserDirection = [] {
            Vector_t<double, 3> value(0.0);
            value(2) = -1.0;
            return value;
        }();
        double sigmaThetaXRad       = 1.0e-3;
        double sigmaThetaYRad       = 1.0e-3;
        double relativeEnergySpread = 0.0;
        double sigmaX_m             = 0.0;
        double sigmaY_m             = 0.0;
        double sigmaS_m             = 0.0;
        double laserRayleighX_m     = 1.0e-6;
        double laserRayleighY_m     = 1.0e-6;
        double laserSigmaT_m        = 5.0e-12 * Physics::c;
        bool overlapWeighting       = false;
        std::size_t energyBins      = 80;
        std::size_t thetaBins       = 80;
        double energyMinGeV         = 0.0;
        double energyMaxGeV         = 0.5;
        double thetaMinRad          = 0.0;
        double thetaMaxRad          = 0.0060;
    };

    struct JointHistogram {
        std::vector<double> energyCentersGeV;
        std::vector<double> thetaCentersRad;
        std::vector<std::vector<double>> densityPerGeVRad;
        std::vector<std::vector<double>> counts;
        double energyBinWidthGeV = 0.0;
        double thetaBinWidthRad  = 0.0;
        double totalWeight       = 0.0;
    };

    inline double polarAngleRad(const Vector_t<double, 3>& momentum) {
        const double transverse = std::hypot(momentum(0), momentum(1));
        return std::atan2(transverse, momentum(2));
    }

    inline double sampledObservable(
            const Physics::LinearBreitWheeler::SampledEvent& event, FinalState state,
            Observable observable) {
        if (observable == Observable::Theta) {
            return state == FinalState::Positron ? polarAngleRad(event.positronMomentumLabGeV)
                                                 : polarAngleRad(event.electronMomentumLabGeV);
        }

        return state == FinalState::Positron ? event.positronEnergyLabGeV
                                             : event.electronEnergyLabGeV;
    }

    struct SampledPhotonBeamState {
        double energyGeV              = 0.0;
        double slopeXRad              = 0.0;
        double slopeYRad              = 0.0;
        double x_m                    = 0.0;
        double y_m                    = 0.0;
        double s_m                    = 0.0;
        Vector_t<double, 3> direction = [] {
            Vector_t<double, 3> value(0.0);
            value(2) = 1.0;
            return value;
        }();
    };

    inline double sampleHighEnergyPhotonEnergyGeV(
            double centralEnergyGeV, double relativeEnergySpread, std::mt19937_64& engine);

    inline void buildTransverseBasis(
            const Vector_t<double, 3>& referenceDirection, Vector_t<double, 3>& axis1,
            Vector_t<double, 3>& axis2) {
        axis1    = Vector_t<double, 3>(0.0);
        axis1(0) = 1.0;
        if (std::abs(referenceDirection(0)) > 0.9) {
            axis1(0) = 0.0;
            axis1(1) = 1.0;
        }
        axis1                  = cross(referenceDirection, axis1);
        const double axis1Norm = std::sqrt(dot(axis1, axis1));
        axis1 /= axis1Norm;
        axis2                  = cross(referenceDirection, axis1);
        const double axis2Norm = std::sqrt(dot(axis2, axis2));
        axis2 /= axis2Norm;
    }

    inline SampledPhotonBeamState samplePhotonBeamState(
            const Vector_t<double, 3>& referenceDirection, double sigmaThetaXRad,
            double sigmaThetaYRad, double sigmaX_m, double sigmaY_m, double sigmaS_m,
            double centralEnergyGeV, double relativeEnergySpread, std::mt19937_64& engine) {
        std::normal_distribution<double> unitNormal(0.0, 1.0);
        SampledPhotonBeamState state;
        state.slopeXRad = sigmaThetaXRad * unitNormal(engine);
        state.slopeYRad = sigmaThetaYRad * unitNormal(engine);
        state.x_m       = sigmaX_m > 0.0 ? sigmaX_m * unitNormal(engine) : 0.0;
        state.y_m       = sigmaY_m > 0.0 ? sigmaY_m * unitNormal(engine) : 0.0;
        state.s_m       = sigmaS_m > 0.0 ? sigmaS_m * unitNormal(engine) : 0.0;
        state.energyGeV =
                sampleHighEnergyPhotonEnergyGeV(centralEnergyGeV, relativeEnergySpread, engine);

        Vector_t<double, 3> axis1(0.0);
        Vector_t<double, 3> axis2(0.0);
        buildTransverseBasis(referenceDirection, axis1, axis2);
        state.direction = referenceDirection + state.slopeXRad * axis1 + state.slopeYRad * axis2;
        const double directionNorm = std::sqrt(dot(state.direction, state.direction));
        state.direction /= directionNorm;
        return state;
    }

    inline Vector_t<double, 3> samplePhotonBeamDirection(
            const Vector_t<double, 3>& referenceDirection, double sigmaThetaXRad,
            double sigmaThetaYRad, std::mt19937_64& engine) {
        return samplePhotonBeamState(
                       referenceDirection, sigmaThetaXRad, sigmaThetaYRad, 0.0, 0.0, 0.0, 1.0, 0.0,
                       engine)
                .direction;
    }

    inline double effectiveSpatialOverlapSigma(double beamSigma_m, double waist_m) {
        if (beamSigma_m <= 0.0) {
            return 0.0;
        }
        if (waist_m <= 0.0) {
            return beamSigma_m;
        }

        const double inverseVariance =
                1.0 / (beamSigma_m * beamSigma_m) + 4.0 / (waist_m * waist_m);
        return std::sqrt(1.0 / inverseVariance);
    }

    inline double effectiveTemporalOverlapSigma(double beamSigmaS_m, double laserSigmaT_m) {
        if (beamSigmaS_m <= 0.0) {
            return 0.0;
        }
        if (laserSigmaT_m <= 0.0) {
            return beamSigmaS_m;
        }

        const double inverseVariance =
                1.0 / (beamSigmaS_m * beamSigmaS_m) + 1.0 / (laserSigmaT_m * laserSigmaT_m);
        return std::sqrt(1.0 / inverseVariance);
    }

    inline SampledPhotonBeamState sampleOverlapPhotonBeamState(
            const FinitePhotonBeamConfig& config, std::mt19937_64& engine) {
        constexpr double pi = 3.14159265358979323846;
        const double waistX_m =
                config.laserRayleighX_m > 0.0
                        ? std::sqrt(config.wavelength_m * config.laserRayleighX_m / pi)
                        : 0.0;
        const double waistY_m =
                config.laserRayleighY_m > 0.0
                        ? std::sqrt(config.wavelength_m * config.laserRayleighY_m / pi)
                        : 0.0;
        return samplePhotonBeamState(
                config.referenceHighEnergyDirection, config.sigmaThetaXRad, config.sigmaThetaYRad,
                effectiveSpatialOverlapSigma(config.sigmaX_m, waistX_m),
                effectiveSpatialOverlapSigma(config.sigmaY_m, waistY_m),
                effectiveTemporalOverlapSigma(config.sigmaS_m, config.laserSigmaT_m),
                config.centralHighEnergyPhotonEnergyGeV, config.relativeEnergySpread, engine);
    }

    inline SampledPhotonBeamState sampleOverlapPhotonBeamState(
            const FinitePhotonBeamJointConfig& config, std::mt19937_64& engine) {
        constexpr double pi = 3.14159265358979323846;
        const double waistX_m =
                config.laserRayleighX_m > 0.0
                        ? std::sqrt(config.wavelength_m * config.laserRayleighX_m / pi)
                        : 0.0;
        const double waistY_m =
                config.laserRayleighY_m > 0.0
                        ? std::sqrt(config.wavelength_m * config.laserRayleighY_m / pi)
                        : 0.0;
        return samplePhotonBeamState(
                config.referenceHighEnergyDirection, config.sigmaThetaXRad, config.sigmaThetaYRad,
                effectiveSpatialOverlapSigma(config.sigmaX_m, waistX_m),
                effectiveSpatialOverlapSigma(config.sigmaY_m, waistY_m),
                effectiveTemporalOverlapSigma(config.sigmaS_m, config.laserSigmaT_m),
                config.centralHighEnergyPhotonEnergyGeV, config.relativeEnergySpread, engine);
    }

    inline double sampleHighEnergyPhotonEnergyGeV(
            double centralEnergyGeV, double relativeEnergySpread, std::mt19937_64& engine) {
        if (relativeEnergySpread <= 0.0) {
            return centralEnergyGeV;
        }

        std::normal_distribution<double> unitNormal(0.0, 1.0);
        double energyGeV = centralEnergyGeV;
        do {
            energyGeV = centralEnergyGeV * (1.0 + relativeEnergySpread * unitNormal(engine));
        } while (energyGeV <= 0.0);
        return energyGeV;
    }

    /**
     * @brief Sample a one-dimensional Breit-Wheeler histogram for a finite incoming photon beam.
     *
     * The laser photon stays fixed. For each sampled event the incoming high-energy
     * photon direction is drawn from a Gaussian angular spread around the reference
     * beam axis, and the energy is optionally smeared by a Gaussian relative energy
     * spread. The resulting per-event kernel is then sampled with the same
     * host-side Breit-Wheeler event generator used in the fixed-geometry benchmark.
     */
    inline Histogram sampleFinitePhotonBeamHistogram(
            const FinitePhotonBeamConfig& config, FinalState state, Observable observable,
            std::size_t sampleCount, std::uint64_t streamIndex = 0) {
        if (config.bins == 0) {
            throw std::runtime_error(
                    "LinearBreitWheelerBenchmark: number of bins must be positive.");
        }
        if (config.maxValue <= config.minValue) {
            throw std::runtime_error(
                    "LinearBreitWheelerBenchmark: histogram range must satisfy max > min.");
        }
        if (sampleCount == 0) {
            throw std::runtime_error("LinearBreitWheelerBenchmark: sample count must be positive.");
        }

        Histogram histogram;
        histogram.binWidth = (config.maxValue - config.minValue) / static_cast<double>(config.bins);
        histogram.centers.resize(config.bins);
        histogram.density.assign(config.bins, 0.0);
        histogram.counts.assign(config.bins, 0.0);
        for (std::size_t i = 0; i < config.bins; ++i) {
            histogram.centers[i] =
                    config.minValue + (static_cast<double>(i) + 0.5) * histogram.binWidth;
        }

        const double laserPhotonEnergyGeV =
                Physics::LinearBreitWheeler::photonEnergyFromWavelengthGeV(config.wavelength_m);
        auto engine           = Physics::LinearBreitWheeler::makeHostRandomEngine(streamIndex);
        histogram.totalWeight = 0.0;

        for (std::size_t i = 0; i < sampleCount; ++i) {
            const auto beamState =
                    config.overlapWeighting
                            ? sampleOverlapPhotonBeamState(config, engine)
                            : samplePhotonBeamState(
                                      config.referenceHighEnergyDirection, config.sigmaThetaXRad,
                                      config.sigmaThetaYRad, config.sigmaX_m, config.sigmaY_m,
                                      config.sigmaS_m, config.centralHighEnergyPhotonEnergyGeV,
                                      config.relativeEnergySpread, engine);
            constexpr double overlapWeight = 1.0;
            histogram.totalWeight += overlapWeight;
            const auto kernel = Physics::LinearBreitWheeler::makeSamplingKernel(
                    beamState.energyGeV, laserPhotonEnergyGeV, beamState.direction,
                    config.laserDirection);
            const auto event   = Physics::LinearBreitWheeler::sampleEvent(kernel, engine);
            const double value = sampledObservable(event, state, observable);
            if (value < config.minValue || value >= config.maxValue) {
                continue;
            }
            const std::size_t bin =
                    static_cast<std::size_t>((value - config.minValue) / histogram.binWidth);
            if (bin < histogram.counts.size()) {
                histogram.counts[bin] += overlapWeight;
            }
        }

        if (histogram.totalWeight <= 0.0) {
            throw std::runtime_error(
                    "LinearBreitWheelerBenchmark: zero accepted finite-photon-beam events.");
        }
        for (std::size_t i = 0; i < histogram.counts.size(); ++i) {
            histogram.density[i] =
                    histogram.counts[i] / (histogram.totalWeight * histogram.binWidth);
        }
        return histogram;
    }

    /**
     * @brief Sample a one-dimensional Breit-Wheeler benchmark histogram.
     *
     * The first CAIN-backed OPALX validation path is intentionally host-only and
     * fixed-geometry. This helper mirrors the standalone benchmark executable: it
     * builds one cached `LinearBreitWheeler::SamplingKernel`, draws a fixed number
     * of sampled events, and bins either the outgoing electron or positron energy
     * or polar angle in the lab frame.
     */
    inline Histogram sampleHistogram(
            const HistogramConfig& config, FinalState state, Observable observable,
            std::size_t sampleCount, std::uint64_t streamIndex = 0) {
        if (config.bins == 0) {
            throw std::runtime_error(
                    "LinearBreitWheelerBenchmark: number of bins must be positive.");
        }
        if (config.maxValue <= config.minValue) {
            throw std::runtime_error(
                    "LinearBreitWheelerBenchmark: histogram range must satisfy max > min.");
        }
        if (sampleCount == 0) {
            throw std::runtime_error("LinearBreitWheelerBenchmark: sample count must be positive.");
        }

        Histogram histogram;
        histogram.binWidth = (config.maxValue - config.minValue) / static_cast<double>(config.bins);
        histogram.centers.resize(config.bins);
        histogram.density.assign(config.bins, 0.0);
        histogram.counts.assign(config.bins, 0.0);
        for (std::size_t i = 0; i < config.bins; ++i) {
            histogram.centers[i] =
                    config.minValue + (static_cast<double>(i) + 0.5) * histogram.binWidth;
        }

        const double laserPhotonEnergyGeV =
                Physics::LinearBreitWheeler::photonEnergyFromWavelengthGeV(config.wavelength_m);
        const auto kernel = Physics::LinearBreitWheeler::makeSamplingKernel(
                config.highEnergyPhotonEnergyGeV, laserPhotonEnergyGeV, config.highEnergyDirection,
                config.laserDirection);
        auto engine = Physics::LinearBreitWheeler::makeHostRandomEngine(streamIndex);

        histogram.totalWeight = static_cast<double>(sampleCount);
        for (std::size_t i = 0; i < sampleCount; ++i) {
            const auto event   = Physics::LinearBreitWheeler::sampleEvent(kernel, engine);
            const double value = sampledObservable(event, state, observable);
            if (value < config.minValue || value >= config.maxValue) {
                continue;
            }
            const std::size_t bin =
                    static_cast<std::size_t>((value - config.minValue) / histogram.binWidth);
            if (bin < histogram.counts.size()) {
                histogram.counts[bin] += 1.0;
            }
        }

        for (std::size_t i = 0; i < histogram.counts.size(); ++i) {
            histogram.density[i] =
                    histogram.counts[i] / (histogram.totalWeight * histogram.binWidth);
        }
        return histogram;
    }

    /**
     * @brief Sample the joint Breit-Wheeler laboratory distribution in energy and polar angle.
     *
     * Each sampled event contributes to a two-dimensional histogram in
     * @f$(E,\theta)@f$ for either the outgoing electron or positron. The resulting
     * density is normalized so that integrating over the full histogram cell area
     * yields unity.
     */

    inline JointHistogram sampleFinitePhotonBeamJointHistogram(
            const FinitePhotonBeamJointConfig& config, FinalState state, std::size_t sampleCount,
            std::uint64_t streamIndex = 0) {
        if (config.energyBins == 0 || config.thetaBins == 0) {
            throw std::runtime_error(
                    "LinearBreitWheelerBenchmark: joint histogram bin counts must be positive.");
        }
        if (config.energyMaxGeV <= config.energyMinGeV
            || config.thetaMaxRad <= config.thetaMinRad) {
            throw std::runtime_error(
                    "LinearBreitWheelerBenchmark: joint histogram ranges must satisfy max > min.");
        }
        if (sampleCount == 0) {
            throw std::runtime_error("LinearBreitWheelerBenchmark: sample count must be positive.");
        }

        JointHistogram histogram;
        histogram.energyBinWidthGeV = (config.energyMaxGeV - config.energyMinGeV)
                                      / static_cast<double>(config.energyBins);
        histogram.thetaBinWidthRad =
                (config.thetaMaxRad - config.thetaMinRad) / static_cast<double>(config.thetaBins);
        histogram.energyCentersGeV.resize(config.energyBins);
        histogram.thetaCentersRad.resize(config.thetaBins);
        histogram.densityPerGeVRad.assign(
                config.energyBins, std::vector<double>(config.thetaBins, 0.0));
        histogram.counts.assign(config.energyBins, std::vector<double>(config.thetaBins, 0.0));

        for (std::size_t i = 0; i < config.energyBins; ++i) {
            histogram.energyCentersGeV[i] =
                    config.energyMinGeV
                    + (static_cast<double>(i) + 0.5) * histogram.energyBinWidthGeV;
        }
        for (std::size_t j = 0; j < config.thetaBins; ++j) {
            histogram.thetaCentersRad[j] =
                    config.thetaMinRad
                    + (static_cast<double>(j) + 0.5) * histogram.thetaBinWidthRad;
        }

        const double laserPhotonEnergyGeV =
                Physics::LinearBreitWheeler::photonEnergyFromWavelengthGeV(config.wavelength_m);
        auto engine           = Physics::LinearBreitWheeler::makeHostRandomEngine(streamIndex);
        histogram.totalWeight = 0.0;

        for (std::size_t i = 0; i < sampleCount; ++i) {
            const auto beamState =
                    config.overlapWeighting
                            ? sampleOverlapPhotonBeamState(config, engine)
                            : samplePhotonBeamState(
                                      config.referenceHighEnergyDirection, config.sigmaThetaXRad,
                                      config.sigmaThetaYRad, config.sigmaX_m, config.sigmaY_m,
                                      config.sigmaS_m, config.centralHighEnergyPhotonEnergyGeV,
                                      config.relativeEnergySpread, engine);
            constexpr double overlapWeight = 1.0;
            histogram.totalWeight += overlapWeight;
            const auto kernel = Physics::LinearBreitWheeler::makeSamplingKernel(
                    beamState.energyGeV, laserPhotonEnergyGeV, beamState.direction,
                    config.laserDirection);
            const auto event       = Physics::LinearBreitWheeler::sampleEvent(kernel, engine);
            const double energyGeV = sampledObservable(event, state, Observable::Energy);
            const double thetaRad  = sampledObservable(event, state, Observable::Theta);
            if (energyGeV < config.energyMinGeV || energyGeV >= config.energyMaxGeV
                || thetaRad < config.thetaMinRad || thetaRad >= config.thetaMaxRad) {
                continue;
            }

            const std::size_t ie = static_cast<std::size_t>(
                    (energyGeV - config.energyMinGeV) / histogram.energyBinWidthGeV);
            const std::size_t jt = static_cast<std::size_t>(
                    (thetaRad - config.thetaMinRad) / histogram.thetaBinWidthRad);
            if (ie < histogram.counts.size() && jt < histogram.counts[ie].size()) {
                histogram.counts[ie][jt] += overlapWeight;
            }
        }

        if (histogram.totalWeight <= 0.0) {
            throw std::runtime_error(
                    "LinearBreitWheelerBenchmark: zero accepted finite-photon-beam events.");
        }
        const double cellArea = histogram.energyBinWidthGeV * histogram.thetaBinWidthRad;
        for (std::size_t i = 0; i < histogram.counts.size(); ++i) {
            for (std::size_t j = 0; j < histogram.counts[i].size(); ++j) {
                histogram.densityPerGeVRad[i][j] =
                        histogram.counts[i][j] / (histogram.totalWeight * cellArea);
            }
        }

        return histogram;
    }

    inline JointHistogram sampleJointHistogram(
            const JointHistogramConfig& config, FinalState state, std::size_t sampleCount,
            std::uint64_t streamIndex = 0) {
        if (config.energyBins == 0 || config.thetaBins == 0) {
            throw std::runtime_error(
                    "LinearBreitWheelerBenchmark: joint histogram bin counts must be positive.");
        }
        if (config.energyMaxGeV <= config.energyMinGeV
            || config.thetaMaxRad <= config.thetaMinRad) {
            throw std::runtime_error(
                    "LinearBreitWheelerBenchmark: joint histogram ranges must satisfy max > min.");
        }
        if (sampleCount == 0) {
            throw std::runtime_error("LinearBreitWheelerBenchmark: sample count must be positive.");
        }

        JointHistogram histogram;
        histogram.energyBinWidthGeV = (config.energyMaxGeV - config.energyMinGeV)
                                      / static_cast<double>(config.energyBins);
        histogram.thetaBinWidthRad =
                (config.thetaMaxRad - config.thetaMinRad) / static_cast<double>(config.thetaBins);
        histogram.energyCentersGeV.resize(config.energyBins);
        histogram.thetaCentersRad.resize(config.thetaBins);
        histogram.densityPerGeVRad.assign(
                config.energyBins, std::vector<double>(config.thetaBins, 0.0));
        histogram.counts.assign(config.energyBins, std::vector<double>(config.thetaBins, 0.0));

        for (std::size_t i = 0; i < config.energyBins; ++i) {
            histogram.energyCentersGeV[i] =
                    config.energyMinGeV
                    + (static_cast<double>(i) + 0.5) * histogram.energyBinWidthGeV;
        }
        for (std::size_t j = 0; j < config.thetaBins; ++j) {
            histogram.thetaCentersRad[j] =
                    config.thetaMinRad
                    + (static_cast<double>(j) + 0.5) * histogram.thetaBinWidthRad;
        }

        const double laserPhotonEnergyGeV =
                Physics::LinearBreitWheeler::photonEnergyFromWavelengthGeV(config.wavelength_m);
        const auto kernel = Physics::LinearBreitWheeler::makeSamplingKernel(
                config.highEnergyPhotonEnergyGeV, laserPhotonEnergyGeV, config.highEnergyDirection,
                config.laserDirection);
        auto engine = Physics::LinearBreitWheeler::makeHostRandomEngine(streamIndex);

        histogram.totalWeight = static_cast<double>(sampleCount);
        for (std::size_t i = 0; i < sampleCount; ++i) {
            const auto event       = Physics::LinearBreitWheeler::sampleEvent(kernel, engine);
            const double energyGeV = sampledObservable(event, state, Observable::Energy);
            const double thetaRad  = sampledObservable(event, state, Observable::Theta);
            if (energyGeV < config.energyMinGeV || energyGeV >= config.energyMaxGeV
                || thetaRad < config.thetaMinRad || thetaRad >= config.thetaMaxRad) {
                continue;
            }

            const std::size_t ie = static_cast<std::size_t>(
                    (energyGeV - config.energyMinGeV) / histogram.energyBinWidthGeV);
            const std::size_t jt = static_cast<std::size_t>(
                    (thetaRad - config.thetaMinRad) / histogram.thetaBinWidthRad);
            if (ie < histogram.counts.size() && jt < histogram.counts[ie].size()) {
                histogram.counts[ie][jt] += 1.0;
            }
        }

        const double cellArea = histogram.energyBinWidthGeV * histogram.thetaBinWidthRad;
        for (std::size_t i = 0; i < histogram.counts.size(); ++i) {
            for (std::size_t j = 0; j < histogram.counts[i].size(); ++j) {
                histogram.densityPerGeVRad[i][j] =
                        histogram.counts[i][j] / (histogram.totalWeight * cellArea);
            }
        }

        return histogram;
    }

    inline void writeHistogramCSV(
            const Histogram& histogram, const std::filesystem::path& outputPath) {
        std::ofstream output(outputPath);
        if (!output) {
            throw std::runtime_error(
                    "LinearBreitWheelerBenchmark: unable to open output file "
                    + outputPath.string());
        }
        output << "# center,density,count\n";
        output.precision(17);
        for (std::size_t i = 0; i < histogram.centers.size(); ++i) {
            output << histogram.centers[i] << ',' << histogram.density[i] << ','
                   << histogram.counts[i] << '\n';
        }
    }

    inline void writeJointHistogramCSV(
            const JointHistogram& histogram, const std::filesystem::path& outputPath) {
        std::ofstream output(outputPath);
        if (!output) {
            throw std::runtime_error(
                    "LinearBreitWheelerBenchmark: unable to open output file "
                    + outputPath.string());
        }
        output << "# energy_center_GeV,theta_center_rad,density_per_GeV_rad,count\n";
        output.precision(17);
        for (std::size_t i = 0; i < histogram.energyCentersGeV.size(); ++i) {
            for (std::size_t j = 0; j < histogram.thetaCentersRad.size(); ++j) {
                output << histogram.energyCentersGeV[i] << ',' << histogram.thetaCentersRad[j]
                       << ',' << histogram.densityPerGeVRad[i][j] << ',' << histogram.counts[i][j]
                       << '\n';
            }
        }
    }

    inline Histogram readHistogramCSV(const std::filesystem::path& inputPath) {
        Histogram histogram;
        std::ifstream input(inputPath);
        if (!input) {
            throw std::runtime_error(
                    "LinearBreitWheelerBenchmark: unable to open input file " + inputPath.string());
        }

        std::string line;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }
            std::stringstream parser(line);
            std::string field;
            std::vector<double> values;
            while (std::getline(parser, field, ',')) {
                values.push_back(std::stod(field));
            }
            if (values.size() < 2) {
                continue;
            }
            histogram.centers.push_back(values[0]);
            histogram.density.push_back(values[1]);
            histogram.counts.push_back(values.size() >= 3 ? values[2] : 0.0);
        }

        if (histogram.centers.size() >= 2) {
            histogram.binWidth = histogram.centers[1] - histogram.centers[0];
        }
        histogram.totalWeight = 0.0;
        for (double value : histogram.density) {
            histogram.totalWeight += value * histogram.binWidth;
        }
        return histogram;
    }

    inline JointHistogram readJointHistogramCSV(const std::filesystem::path& inputPath) {
        JointHistogram histogram;
        std::ifstream input(inputPath);
        if (!input) {
            throw std::runtime_error(
                    "LinearBreitWheelerBenchmark: unable to open input file " + inputPath.string());
        }

        struct Row {
            double energyCenter = 0.0;
            double thetaCenter  = 0.0;
            double density      = 0.0;
            double count        = 0.0;
        };
        std::vector<Row> rows;
        std::string line;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }
            std::stringstream parser(line);
            std::string field;
            std::vector<double> values;
            while (std::getline(parser, field, ',')) {
                values.push_back(std::stod(field));
            }
            if (values.size() < 3) {
                continue;
            }
            rows.push_back(
                    Row{values[0], values[1], values[2], values.size() >= 4 ? values[3] : 0.0});
            if (histogram.energyCentersGeV.empty()
                || histogram.energyCentersGeV.back() != values[0]) {
                histogram.energyCentersGeV.push_back(values[0]);
            }
            bool seenTheta = false;
            for (double existing : histogram.thetaCentersRad) {
                if (existing == values[1]) {
                    seenTheta = true;
                    break;
                }
            }
            if (!seenTheta) {
                histogram.thetaCentersRad.push_back(values[1]);
            }
        }

        histogram.densityPerGeVRad.assign(
                histogram.energyCentersGeV.size(),
                std::vector<double>(histogram.thetaCentersRad.size(), 0.0));
        histogram.counts.assign(
                histogram.energyCentersGeV.size(),
                std::vector<double>(histogram.thetaCentersRad.size(), 0.0));
        for (const auto& row : rows) {
            std::size_t i = 0;
            while (i < histogram.energyCentersGeV.size()
                   && histogram.energyCentersGeV[i] != row.energyCenter) {
                ++i;
            }
            std::size_t j = 0;
            while (j < histogram.thetaCentersRad.size()
                   && histogram.thetaCentersRad[j] != row.thetaCenter) {
                ++j;
            }
            if (i < histogram.energyCentersGeV.size() && j < histogram.thetaCentersRad.size()) {
                histogram.densityPerGeVRad[i][j] = row.density;
                histogram.counts[i][j]           = row.count;
            }
        }

        if (histogram.energyCentersGeV.size() >= 2) {
            histogram.energyBinWidthGeV =
                    histogram.energyCentersGeV[1] - histogram.energyCentersGeV[0];
        }
        if (histogram.thetaCentersRad.size() >= 2) {
            histogram.thetaBinWidthRad =
                    histogram.thetaCentersRad[1] - histogram.thetaCentersRad[0];
        }
        histogram.totalWeight = 0.0;
        for (const auto& row : histogram.densityPerGeVRad) {
            for (double value : row) {
                histogram.totalWeight +=
                        value * histogram.energyBinWidthGeV * histogram.thetaBinWidthRad;
            }
        }
        return histogram;
    }

    inline double histogramArea(const Histogram& histogram) {
        double area = 0.0;
        for (double density : histogram.density) {
            area += density * histogram.binWidth;
        }
        return area;
    }

    inline double histogramMean(const Histogram& histogram) {
        double mean = 0.0;
        for (std::size_t i = 0; i < histogram.centers.size(); ++i) {
            mean += histogram.centers[i] * histogram.density[i] * histogram.binWidth;
        }
        return mean;
    }

    inline double histogramL1Distance(const Histogram& lhs, const Histogram& rhs) {
        if (lhs.centers.size() != rhs.centers.size()) {
            throw std::runtime_error("LinearBreitWheelerBenchmark: histogram sizes do not match.");
        }
        double distance = 0.0;
        for (std::size_t i = 0; i < lhs.centers.size(); ++i) {
            distance += std::abs(lhs.density[i] - rhs.density[i]) * lhs.binWidth;
        }
        return distance;
    }

    inline double jointHistogramArea(const JointHistogram& histogram) {
        double area = 0.0;
        for (const auto& row : histogram.densityPerGeVRad) {
            for (double density : row) {
                area += density * histogram.energyBinWidthGeV * histogram.thetaBinWidthRad;
            }
        }
        return area;
    }

    inline double jointHistogramMeanEnergyGeV(const JointHistogram& histogram) {
        double mean = 0.0;
        for (std::size_t i = 0; i < histogram.energyCentersGeV.size(); ++i) {
            for (std::size_t j = 0; j < histogram.thetaCentersRad.size(); ++j) {
                mean += histogram.energyCentersGeV[i] * histogram.densityPerGeVRad[i][j]
                        * histogram.energyBinWidthGeV * histogram.thetaBinWidthRad;
            }
        }
        return mean;
    }

    inline double jointHistogramMeanThetaRad(const JointHistogram& histogram) {
        double mean = 0.0;
        for (std::size_t i = 0; i < histogram.energyCentersGeV.size(); ++i) {
            for (std::size_t j = 0; j < histogram.thetaCentersRad.size(); ++j) {
                mean += histogram.thetaCentersRad[j] * histogram.densityPerGeVRad[i][j]
                        * histogram.energyBinWidthGeV * histogram.thetaBinWidthRad;
            }
        }
        return mean;
    }

    inline double jointHistogramL1Distance(const JointHistogram& lhs, const JointHistogram& rhs) {
        if (lhs.energyCentersGeV.size() != rhs.energyCentersGeV.size()
            || lhs.thetaCentersRad.size() != rhs.thetaCentersRad.size()) {
            throw std::runtime_error(
                    "LinearBreitWheelerBenchmark: joint histogram sizes do not match.");
        }
        double distance = 0.0;
        for (std::size_t i = 0; i < lhs.energyCentersGeV.size(); ++i) {
            for (std::size_t j = 0; j < lhs.thetaCentersRad.size(); ++j) {
                distance += std::abs(lhs.densityPerGeVRad[i][j] - rhs.densityPerGeVRad[i][j])
                            * lhs.energyBinWidthGeV * lhs.thetaBinWidthRad;
            }
        }
        return distance;
    }

}  // namespace LinearBreitWheelerBenchmark

#endif
