#ifndef OPALX_TEST_LINEAR_COMPTON_SPECTRUM_COMMON_H
#define OPALX_TEST_LINEAR_COMPTON_SPECTRUM_COMMON_H

#include "Physics/LinearCompton.h"
#include "Physics/Physics.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace LinearComptonBenchmark {

    struct SpectrumConfig {
        double electronTotalEnergyGeV     = 1.0;
        double wavelength_m               = 1.03e-6;
        Vector_t<double, 3> beamDirection = [] {
            Vector_t<double, 3> value(0.0);
            value(2) = 1.0;
            return value;
        }();
        Vector_t<double, 3> laserDirection = [] {
            Vector_t<double, 3> value(0.0);
            value(0) = 1.0;
            return value;
        }();
        std::size_t bins    = 80;
        double energyMinGeV = 0.0;
        double energyMaxGeV = 0.01;
        int cosinePanels    = 720;
        int azimuthPanels   = 720;
    };

    struct SpectrumHistogram {
        std::vector<double> centersGeV;
        std::vector<double> densityPerGeV;
        std::vector<double> counts;
        double binWidthGeV = 0.0;
        double totalWeight = 0.0;
    };

    inline SpectrumHistogram integrateLabSpectrum(const SpectrumConfig& config) {
        SpectrumHistogram histogram;
        histogram.binWidthGeV =
                (config.energyMaxGeV - config.energyMinGeV) / static_cast<double>(config.bins);
        histogram.centersGeV.resize(config.bins);
        histogram.densityPerGeV.assign(config.bins, 0.0);
        histogram.counts.assign(config.bins, 0.0);

        for (std::size_t i = 0; i < config.bins; ++i) {
            histogram.centersGeV[i] =
                    config.energyMinGeV + (static_cast<double>(i) + 0.5) * histogram.binWidthGeV;
        }

        const double laserPhotonEnergyGeV =
                Physics::LinearCompton::photonEnergyFromWavelengthGeV(config.wavelength_m);
        const double incomingPhotonEnergyERFGeV =
                Physics::LinearCompton::restFrameIncomingPhotonEnergyGeV(
                        config.electronTotalEnergyGeV, laserPhotonEnergyGeV, config.beamDirection,
                        config.laserDirection);

        const double dMu  = 2.0 / static_cast<double>(config.cosinePanels);
        const double dPhi = Physics::two_pi / static_cast<double>(config.azimuthPanels);

        for (int i = 0; i < config.cosinePanels; ++i) {
            const double mu     = -1.0 + (static_cast<double>(i) + 0.5) * dMu;
            const double kernel = Physics::LinearCompton::differentialCrossSectionSolidAngleERF(
                    incomingPhotonEnergyERFGeV, mu);

            for (int j = 0; j < config.azimuthPanels; ++j) {
                const double phi       = (static_cast<double>(j) + 0.5) * dPhi;
                const double energyGeV = Physics::LinearCompton::labPhotonEnergyGeV(
                        config.electronTotalEnergyGeV, laserPhotonEnergyGeV, config.beamDirection,
                        config.laserDirection, mu, phi);
                const double weight = kernel * dMu * dPhi;
                histogram.totalWeight += weight;

                if (energyGeV < config.energyMinGeV || energyGeV >= config.energyMaxGeV) {
                    continue;
                }

                const std::size_t bin = static_cast<std::size_t>(
                        (energyGeV - config.energyMinGeV) / histogram.binWidthGeV);
                if (bin < histogram.counts.size()) {
                    histogram.counts[bin] += weight;
                }
            }
        }

        if (histogram.totalWeight <= 0.0) {
            throw std::runtime_error(
                    "LinearComptonBenchmark: total histogram weight is not positive.");
        }

        for (std::size_t i = 0; i < histogram.densityPerGeV.size(); ++i) {
            histogram.densityPerGeV[i] =
                    histogram.counts[i] / (histogram.totalWeight * histogram.binWidthGeV);
        }

        return histogram;
    }

    inline SpectrumHistogram sampleLabSpectrum(
            const SpectrumConfig& config, std::size_t sampleCount, std::uint64_t streamIndex = 0) {
        SpectrumHistogram histogram;
        histogram.binWidthGeV =
                (config.energyMaxGeV - config.energyMinGeV) / static_cast<double>(config.bins);
        histogram.centersGeV.resize(config.bins);
        histogram.densityPerGeV.assign(config.bins, 0.0);
        histogram.counts.assign(config.bins, 0.0);

        for (std::size_t i = 0; i < config.bins; ++i) {
            histogram.centersGeV[i] =
                    config.energyMinGeV + (static_cast<double>(i) + 0.5) * histogram.binWidthGeV;
        }

        const double laserPhotonEnergyGeV =
                Physics::LinearCompton::photonEnergyFromWavelengthGeV(config.wavelength_m);
        const auto kernel = Physics::LinearCompton::makeSamplingKernel(
                config.electronTotalEnergyGeV, laserPhotonEnergyGeV, config.beamDirection,
                config.laserDirection);
        auto engine = Physics::LinearCompton::makeHostRandomEngine(streamIndex);

        histogram.totalWeight = static_cast<double>(sampleCount);
        for (std::size_t i = 0; i < sampleCount; ++i) {
            const auto event       = Physics::LinearCompton::sampleEvent(kernel, engine);
            const double energyGeV = event.scatteredPhotonEnergyLabGeV;
            if (energyGeV < config.energyMinGeV || energyGeV >= config.energyMaxGeV) {
                continue;
            }

            const std::size_t bin = static_cast<std::size_t>(
                    (energyGeV - config.energyMinGeV) / histogram.binWidthGeV);
            if (bin < histogram.counts.size()) {
                histogram.counts[bin] += 1.0;
            }
        }

        if (histogram.totalWeight <= 0.0) {
            throw std::runtime_error("LinearComptonBenchmark: sample count must be positive.");
        }

        for (std::size_t i = 0; i < histogram.densityPerGeV.size(); ++i) {
            histogram.densityPerGeV[i] =
                    histogram.counts[i] / (histogram.totalWeight * histogram.binWidthGeV);
        }

        return histogram;
    }

    inline void writeSpectrumCSV(
            const SpectrumHistogram& histogram, const std::filesystem::path& outputPath) {
        std::ofstream output(outputPath);
        if (!output) {
            throw std::runtime_error(
                    "LinearComptonBenchmark: unable to open output file " + outputPath.string());
        }

        output << "# center_GeV,density_per_GeV,count\n";
        output << std::setprecision(17);
        for (std::size_t i = 0; i < histogram.centersGeV.size(); ++i) {
            output << histogram.centersGeV[i] << ',' << histogram.densityPerGeV[i] << ','
                   << histogram.counts[i] << '\n';
        }
    }

    inline SpectrumHistogram readSpectrumCSV(const std::filesystem::path& inputPath) {
        SpectrumHistogram histogram;
        std::ifstream input(inputPath);
        if (!input) {
            throw std::runtime_error(
                    "LinearComptonBenchmark: unable to open input file " + inputPath.string());
        }

        std::string line;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }

            std::stringstream parser(line);
            std::string field;
            std::vector<double> fields;
            while (std::getline(parser, field, ',')) {
                fields.push_back(std::stod(field));
            }
            if (fields.size() < 2) {
                continue;
            }

            histogram.centersGeV.push_back(fields[0]);
            histogram.densityPerGeV.push_back(fields[1]);
            histogram.counts.push_back(fields.size() >= 3 ? fields[2] : 0.0);
        }

        if (histogram.centersGeV.size() >= 2) {
            histogram.binWidthGeV = histogram.centersGeV[1] - histogram.centersGeV[0];
        }
        histogram.totalWeight = 0.0;
        for (double density : histogram.densityPerGeV) {
            histogram.totalWeight += density * histogram.binWidthGeV;
        }

        return histogram;
    }

    inline double histogramMeanEnergyGeV(const SpectrumHistogram& histogram) {
        double mean = 0.0;
        for (std::size_t i = 0; i < histogram.centersGeV.size(); ++i) {
            mean += histogram.centersGeV[i] * histogram.densityPerGeV[i] * histogram.binWidthGeV;
        }
        return mean;
    }

    inline double histogramArea(const SpectrumHistogram& histogram) {
        double area = 0.0;
        for (double density : histogram.densityPerGeV) {
            area += density * histogram.binWidthGeV;
        }
        return area;
    }

    inline double histogramL1Distance(const SpectrumHistogram& lhs, const SpectrumHistogram& rhs) {
        if (lhs.centersGeV.size() != rhs.centersGeV.size()) {
            throw std::runtime_error("LinearComptonBenchmark: histogram sizes do not match.");
        }

        double distance = 0.0;
        for (std::size_t i = 0; i < lhs.centersGeV.size(); ++i) {
            distance += std::abs(lhs.densityPerGeV[i] - rhs.densityPerGeV[i]) * lhs.binWidthGeV;
        }
        return distance;
    }

    struct AngleConfig {
        double electronTotalEnergyGeV     = 1.0;
        double wavelength_m               = 1.03e-6;
        Vector_t<double, 3> beamDirection = [] {
            Vector_t<double, 3> value(0.0);
            value(2) = 1.0;
            return value;
        }();
        Vector_t<double, 3> laserDirection = [] {
            Vector_t<double, 3> value(0.0);
            value(0) = 1.0;
            return value;
        }();
        std::size_t bins   = 80;
        double thetaMinRad = 0.0;
        double thetaMaxRad = 0.02;
        int cosinePanels   = 720;
        int azimuthPanels  = 720;
    };

    struct AngleHistogram {
        std::vector<double> centersRad;
        std::vector<double> densityPerRad;
        std::vector<double> counts;
        double binWidthRad = 0.0;
        double totalWeight = 0.0;
    };

    inline double photonPolarAngleRad(
            const Vector_t<double, 3>& photonDirection, const Vector_t<double, 3>& beamDirection) {
        const double beamNorm   = std::sqrt(dot(beamDirection, beamDirection));
        const double photonNorm = std::sqrt(dot(photonDirection, photonDirection));
        if (beamNorm <= 0.0 || photonNorm <= 0.0) {
            throw std::runtime_error(
                    "LinearComptonBenchmark: beam/photon direction norm must be positive.");
        }
        const double cosine = std::max(
                -1.0, std::min(1.0, dot(photonDirection, beamDirection) / (beamNorm * photonNorm)));
        return std::acos(cosine);
    }

    inline AngleHistogram integrateLabAngularSpectrum(const AngleConfig& config) {
        AngleHistogram histogram;
        histogram.binWidthRad =
                (config.thetaMaxRad - config.thetaMinRad) / static_cast<double>(config.bins);
        histogram.centersRad.resize(config.bins);
        histogram.densityPerRad.assign(config.bins, 0.0);
        histogram.counts.assign(config.bins, 0.0);

        for (std::size_t i = 0; i < config.bins; ++i) {
            histogram.centersRad[i] =
                    config.thetaMinRad + (static_cast<double>(i) + 0.5) * histogram.binWidthRad;
        }

        const double laserPhotonEnergyGeV =
                Physics::LinearCompton::photonEnergyFromWavelengthGeV(config.wavelength_m);
        const double incomingPhotonEnergyERFGeV =
                Physics::LinearCompton::restFrameIncomingPhotonEnergyGeV(
                        config.electronTotalEnergyGeV, laserPhotonEnergyGeV, config.beamDirection,
                        config.laserDirection);

        const double dMu  = 2.0 / static_cast<double>(config.cosinePanels);
        const double dPhi = Physics::two_pi / static_cast<double>(config.azimuthPanels);

        for (int i = 0; i < config.cosinePanels; ++i) {
            const double mu     = -1.0 + (static_cast<double>(i) + 0.5) * dMu;
            const double kernel = Physics::LinearCompton::differentialCrossSectionSolidAngleERF(
                    incomingPhotonEnergyERFGeV, mu);

            for (int j = 0; j < config.azimuthPanels; ++j) {
                const double phi                    = (static_cast<double>(j) + 0.5) * dPhi;
                const Vector_t<double, 3> direction = Physics::LinearCompton::labPhotonDirection(
                        config.electronTotalEnergyGeV, config.beamDirection, config.laserDirection,
                        mu, phi);
                const double thetaRad = photonPolarAngleRad(direction, config.beamDirection);
                const double weight   = kernel * dMu * dPhi;
                histogram.totalWeight += weight;

                if (thetaRad < config.thetaMinRad || thetaRad >= config.thetaMaxRad) {
                    continue;
                }

                const std::size_t bin = static_cast<std::size_t>(
                        (thetaRad - config.thetaMinRad) / histogram.binWidthRad);
                if (bin < histogram.counts.size()) {
                    histogram.counts[bin] += weight;
                }
            }
        }

        if (histogram.totalWeight <= 0.0) {
            throw std::runtime_error(
                    "LinearComptonBenchmark: angular histogram total weight is not positive.");
        }

        for (std::size_t i = 0; i < histogram.densityPerRad.size(); ++i) {
            histogram.densityPerRad[i] =
                    histogram.counts[i] / (histogram.totalWeight * histogram.binWidthRad);
        }

        return histogram;
    }

    inline AngleHistogram sampleLabAngularSpectrum(
            const AngleConfig& config, std::size_t sampleCount, std::uint64_t streamIndex = 0) {
        AngleHistogram histogram;
        histogram.binWidthRad =
                (config.thetaMaxRad - config.thetaMinRad) / static_cast<double>(config.bins);
        histogram.centersRad.resize(config.bins);
        histogram.densityPerRad.assign(config.bins, 0.0);
        histogram.counts.assign(config.bins, 0.0);

        for (std::size_t i = 0; i < config.bins; ++i) {
            histogram.centersRad[i] =
                    config.thetaMinRad + (static_cast<double>(i) + 0.5) * histogram.binWidthRad;
        }

        const double laserPhotonEnergyGeV =
                Physics::LinearCompton::photonEnergyFromWavelengthGeV(config.wavelength_m);
        const auto kernel = Physics::LinearCompton::makeSamplingKernel(
                config.electronTotalEnergyGeV, laserPhotonEnergyGeV, config.beamDirection,
                config.laserDirection);
        auto engine = Physics::LinearCompton::makeHostRandomEngine(streamIndex);

        histogram.totalWeight = static_cast<double>(sampleCount);
        for (std::size_t i = 0; i < sampleCount; ++i) {
            const auto event = Physics::LinearCompton::sampleEvent(kernel, engine);
            const double thetaRad =
                    photonPolarAngleRad(event.scatteredPhotonDirectionLab, config.beamDirection);
            if (thetaRad < config.thetaMinRad || thetaRad >= config.thetaMaxRad) {
                continue;
            }

            const std::size_t bin = static_cast<std::size_t>(
                    (thetaRad - config.thetaMinRad) / histogram.binWidthRad);
            if (bin < histogram.counts.size()) {
                histogram.counts[bin] += 1.0;
            }
        }

        if (histogram.totalWeight <= 0.0) {
            throw std::runtime_error(
                    "LinearComptonBenchmark: angular sample count must be positive.");
        }

        for (std::size_t i = 0; i < histogram.densityPerRad.size(); ++i) {
            histogram.densityPerRad[i] =
                    histogram.counts[i] / (histogram.totalWeight * histogram.binWidthRad);
        }

        return histogram;
    }

    inline void writeAngleCSV(
            const AngleHistogram& histogram, const std::filesystem::path& outputPath) {
        std::ofstream output(outputPath);
        if (!output) {
            throw std::runtime_error(
                    "LinearComptonBenchmark: unable to open angular output file "
                    + outputPath.string());
        }

        output << "# center_rad,density_per_rad,count\n";
        output << std::setprecision(17);
        for (std::size_t i = 0; i < histogram.centersRad.size(); ++i) {
            output << histogram.centersRad[i] << ',' << histogram.densityPerRad[i] << ','
                   << histogram.counts[i] << '\n';
        }
    }

    inline AngleHistogram readAngleCSV(const std::filesystem::path& inputPath) {
        AngleHistogram histogram;
        std::ifstream input(inputPath);
        if (!input) {
            throw std::runtime_error(
                    "LinearComptonBenchmark: unable to open angular input file "
                    + inputPath.string());
        }

        std::string line;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }

            std::stringstream parser(line);
            std::string field;
            std::vector<double> fields;
            while (std::getline(parser, field, ',')) {
                fields.push_back(std::stod(field));
            }
            if (fields.size() < 2) {
                continue;
            }

            histogram.centersRad.push_back(fields[0]);
            histogram.densityPerRad.push_back(fields[1]);
            histogram.counts.push_back(fields.size() >= 3 ? fields[2] : 0.0);
        }

        if (histogram.centersRad.size() >= 2) {
            histogram.binWidthRad = histogram.centersRad[1] - histogram.centersRad[0];
        }
        histogram.totalWeight = 0.0;
        for (double density : histogram.densityPerRad) {
            histogram.totalWeight += density * histogram.binWidthRad;
        }

        return histogram;
    }

    inline double angleHistogramMeanRad(const AngleHistogram& histogram) {
        double mean = 0.0;
        for (std::size_t i = 0; i < histogram.centersRad.size(); ++i) {
            mean += histogram.centersRad[i] * histogram.densityPerRad[i] * histogram.binWidthRad;
        }
        return mean;
    }

    inline double angleHistogramArea(const AngleHistogram& histogram) {
        double area = 0.0;
        for (double density : histogram.densityPerRad) {
            area += density * histogram.binWidthRad;
        }
        return area;
    }

    inline double angleHistogramL1Distance(const AngleHistogram& lhs, const AngleHistogram& rhs) {
        if (lhs.centersRad.size() != rhs.centersRad.size()) {
            throw std::runtime_error(
                    "LinearComptonBenchmark: angular histogram sizes do not match.");
        }

        double distance = 0.0;
        for (std::size_t i = 0; i < lhs.centersRad.size(); ++i) {
            distance += std::abs(lhs.densityPerRad[i] - rhs.densityPerRad[i]) * lhs.binWidthRad;
        }
        return distance;
    }

    struct JointConfig {
        double electronTotalEnergyGeV     = 1.0;
        double wavelength_m               = 1.03e-6;
        Vector_t<double, 3> beamDirection = [] {
            Vector_t<double, 3> value(0.0);
            value(2) = 1.0;
            return value;
        }();
        Vector_t<double, 3> laserDirection = [] {
            Vector_t<double, 3> value(0.0);
            value(0) = 1.0;
            return value;
        }();
        std::size_t energyBins = 60;
        std::size_t thetaBins  = 60;
        double energyMinGeV    = 0.0;
        double energyMaxGeV    = 0.01;
        double thetaMinRad     = 0.0;
        double thetaMaxRad     = 0.02;
        int cosinePanels       = 720;
        int azimuthPanels      = 720;
    };

    struct JointHistogram {
        std::vector<double> energyCentersGeV;
        std::vector<double> thetaCentersRad;
        std::vector<double> densityPerGeVRad;
        std::vector<double> counts;
        double energyBinWidthGeV = 0.0;
        double thetaBinWidthRad  = 0.0;
        double totalWeight       = 0.0;
    };

    inline std::size_t jointHistogramIndex(
            const JointHistogram& histogram, std::size_t energyBin, std::size_t thetaBin) {
        return energyBin * histogram.thetaCentersRad.size() + thetaBin;
    }

    inline JointHistogram integrateLabJointSpectrum(const JointConfig& config) {
        JointHistogram histogram;
        histogram.energyBinWidthGeV = (config.energyMaxGeV - config.energyMinGeV)
                                      / static_cast<double>(config.energyBins);
        histogram.thetaBinWidthRad =
                (config.thetaMaxRad - config.thetaMinRad) / static_cast<double>(config.thetaBins);
        histogram.energyCentersGeV.resize(config.energyBins);
        histogram.thetaCentersRad.resize(config.thetaBins);
        histogram.densityPerGeVRad.assign(config.energyBins * config.thetaBins, 0.0);
        histogram.counts.assign(config.energyBins * config.thetaBins, 0.0);

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
                Physics::LinearCompton::photonEnergyFromWavelengthGeV(config.wavelength_m);
        const double incomingPhotonEnergyERFGeV =
                Physics::LinearCompton::restFrameIncomingPhotonEnergyGeV(
                        config.electronTotalEnergyGeV, laserPhotonEnergyGeV, config.beamDirection,
                        config.laserDirection);

        const double dMu  = 2.0 / static_cast<double>(config.cosinePanels);
        const double dPhi = Physics::two_pi / static_cast<double>(config.azimuthPanels);

        for (int i = 0; i < config.cosinePanels; ++i) {
            const double mu     = -1.0 + (static_cast<double>(i) + 0.5) * dMu;
            const double kernel = Physics::LinearCompton::differentialCrossSectionSolidAngleERF(
                    incomingPhotonEnergyERFGeV, mu);

            for (int j = 0; j < config.azimuthPanels; ++j) {
                const double phi       = (static_cast<double>(j) + 0.5) * dPhi;
                const double energyGeV = Physics::LinearCompton::labPhotonEnergyGeV(
                        config.electronTotalEnergyGeV, laserPhotonEnergyGeV, config.beamDirection,
                        config.laserDirection, mu, phi);
                const Vector_t<double, 3> direction = Physics::LinearCompton::labPhotonDirection(
                        config.electronTotalEnergyGeV, config.beamDirection, config.laserDirection,
                        mu, phi);
                const double thetaRad = photonPolarAngleRad(direction, config.beamDirection);
                const double weight   = kernel * dMu * dPhi;
                histogram.totalWeight += weight;

                if (energyGeV < config.energyMinGeV || energyGeV >= config.energyMaxGeV
                    || thetaRad < config.thetaMinRad || thetaRad >= config.thetaMaxRad) {
                    continue;
                }

                const std::size_t energyBin = static_cast<std::size_t>(
                        (energyGeV - config.energyMinGeV) / histogram.energyBinWidthGeV);
                const std::size_t thetaBin = static_cast<std::size_t>(
                        (thetaRad - config.thetaMinRad) / histogram.thetaBinWidthRad);
                if (energyBin < histogram.energyCentersGeV.size()
                    && thetaBin < histogram.thetaCentersRad.size()) {
                    histogram.counts[jointHistogramIndex(histogram, energyBin, thetaBin)] += weight;
                }
            }
        }

        if (histogram.totalWeight <= 0.0) {
            throw std::runtime_error(
                    "LinearComptonBenchmark: joint histogram total weight is not positive.");
        }

        const double cellArea = histogram.energyBinWidthGeV * histogram.thetaBinWidthRad;
        for (std::size_t i = 0; i < histogram.densityPerGeVRad.size(); ++i) {
            histogram.densityPerGeVRad[i] =
                    histogram.counts[i] / (histogram.totalWeight * cellArea);
        }

        return histogram;
    }

    inline JointHistogram sampleLabJointSpectrum(
            const JointConfig& config, std::size_t sampleCount, std::uint64_t streamIndex = 0) {
        JointHistogram histogram;
        histogram.energyBinWidthGeV = (config.energyMaxGeV - config.energyMinGeV)
                                      / static_cast<double>(config.energyBins);
        histogram.thetaBinWidthRad =
                (config.thetaMaxRad - config.thetaMinRad) / static_cast<double>(config.thetaBins);
        histogram.energyCentersGeV.resize(config.energyBins);
        histogram.thetaCentersRad.resize(config.thetaBins);
        histogram.densityPerGeVRad.assign(config.energyBins * config.thetaBins, 0.0);
        histogram.counts.assign(config.energyBins * config.thetaBins, 0.0);

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
                Physics::LinearCompton::photonEnergyFromWavelengthGeV(config.wavelength_m);
        const auto kernel = Physics::LinearCompton::makeSamplingKernel(
                config.electronTotalEnergyGeV, laserPhotonEnergyGeV, config.beamDirection,
                config.laserDirection);
        auto engine = Physics::LinearCompton::makeHostRandomEngine(streamIndex);

        histogram.totalWeight = static_cast<double>(sampleCount);
        for (std::size_t i = 0; i < sampleCount; ++i) {
            const auto event       = Physics::LinearCompton::sampleEvent(kernel, engine);
            const double energyGeV = event.scatteredPhotonEnergyLabGeV;
            const double thetaRad =
                    photonPolarAngleRad(event.scatteredPhotonDirectionLab, config.beamDirection);

            if (energyGeV < config.energyMinGeV || energyGeV >= config.energyMaxGeV
                || thetaRad < config.thetaMinRad || thetaRad >= config.thetaMaxRad) {
                continue;
            }

            const std::size_t energyBin = static_cast<std::size_t>(
                    (energyGeV - config.energyMinGeV) / histogram.energyBinWidthGeV);
            const std::size_t thetaBin = static_cast<std::size_t>(
                    (thetaRad - config.thetaMinRad) / histogram.thetaBinWidthRad);
            if (energyBin < histogram.energyCentersGeV.size()
                && thetaBin < histogram.thetaCentersRad.size()) {
                histogram.counts[jointHistogramIndex(histogram, energyBin, thetaBin)] += 1.0;
            }
        }

        if (histogram.totalWeight <= 0.0) {
            throw std::runtime_error(
                    "LinearComptonBenchmark: sampled joint histogram total weight is not "
                    "positive.");
        }

        const double cellArea = histogram.energyBinWidthGeV * histogram.thetaBinWidthRad;
        for (std::size_t i = 0; i < histogram.densityPerGeVRad.size(); ++i) {
            histogram.densityPerGeVRad[i] =
                    histogram.counts[i] / (histogram.totalWeight * cellArea);
        }

        return histogram;
    }

    inline void writeJointCSV(
            const JointHistogram& histogram, const std::filesystem::path& outputPath) {
        std::ofstream output(outputPath);
        if (!output) {
            throw std::runtime_error(
                    "LinearComptonBenchmark: unable to open joint output file "
                    + outputPath.string());
        }

        output << std::setprecision(17);
        output << "# energy_center_GeV,theta_center_rad,density_per_GeV_rad,count\n";
        for (std::size_t i = 0; i < histogram.energyCentersGeV.size(); ++i) {
            for (std::size_t j = 0; j < histogram.thetaCentersRad.size(); ++j) {
                const std::size_t index = jointHistogramIndex(histogram, i, j);
                output << histogram.energyCentersGeV[i] << ',' << histogram.thetaCentersRad[j]
                       << ',' << histogram.densityPerGeVRad[index] << ',' << histogram.counts[index]
                       << '\n';
            }
        }
    }

    inline JointHistogram readJointCSV(const std::filesystem::path& inputPath) {
        struct Row {
            double energyCenter = 0.0;
            double thetaCenter  = 0.0;
            double density      = 0.0;
            double count        = 0.0;
        };

        JointHistogram histogram;
        std::ifstream input(inputPath);
        if (!input) {
            throw std::runtime_error(
                    "LinearComptonBenchmark: unable to open joint input file "
                    + inputPath.string());
        }

        std::vector<Row> rows;
        std::string line;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }

            std::stringstream parser(line);
            std::string field;
            std::vector<double> fields;
            while (std::getline(parser, field, ',')) {
                fields.push_back(std::stod(field));
            }
            if (fields.size() < 3) {
                continue;
            }

            rows.push_back(
                    Row{fields[0], fields[1], fields[2], fields.size() >= 4 ? fields[3] : 0.0});
        }

        auto containsCenter = [](const std::vector<double>& values, double target) {
            for (double value : values) {
                if (std::abs(value - target)
                    <= 1.0e-12 * std::max(1.0, std::max(std::abs(value), std::abs(target)))) {
                    return true;
                }
            }
            return false;
        };

        for (const Row& row : rows) {
            if (!containsCenter(histogram.energyCentersGeV, row.energyCenter)) {
                histogram.energyCentersGeV.push_back(row.energyCenter);
            }
            if (!containsCenter(histogram.thetaCentersRad, row.thetaCenter)) {
                histogram.thetaCentersRad.push_back(row.thetaCenter);
            }
        }

        histogram.densityPerGeVRad.assign(
                histogram.energyCentersGeV.size() * histogram.thetaCentersRad.size(), 0.0);
        histogram.counts.assign(
                histogram.energyCentersGeV.size() * histogram.thetaCentersRad.size(), 0.0);

        auto findIndex = [](const std::vector<double>& values, double target) {
            for (std::size_t i = 0; i < values.size(); ++i) {
                if (std::abs(values[i] - target)
                    <= 1.0e-12 * std::max(1.0, std::max(std::abs(values[i]), std::abs(target)))) {
                    return i;
                }
            }
            throw std::runtime_error(
                    "LinearComptonBenchmark: unable to locate joint histogram center.");
        };

        for (const Row& row : rows) {
            const std::size_t energyBin = findIndex(histogram.energyCentersGeV, row.energyCenter);
            const std::size_t thetaBin  = findIndex(histogram.thetaCentersRad, row.thetaCenter);
            const std::size_t index     = jointHistogramIndex(histogram, energyBin, thetaBin);
            histogram.densityPerGeVRad[index] = row.density;
            histogram.counts[index]           = row.count;
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
        const double cellArea = histogram.energyBinWidthGeV * histogram.thetaBinWidthRad;
        for (double density : histogram.densityPerGeVRad) {
            histogram.totalWeight += density * cellArea;
        }

        return histogram;
    }

    inline double jointHistogramArea(const JointHistogram& histogram) {
        double area           = 0.0;
        const double cellArea = histogram.energyBinWidthGeV * histogram.thetaBinWidthRad;
        for (double density : histogram.densityPerGeVRad) {
            area += density * cellArea;
        }
        return area;
    }

    inline double jointHistogramMeanEnergyGeV(const JointHistogram& histogram) {
        double mean           = 0.0;
        const double cellArea = histogram.energyBinWidthGeV * histogram.thetaBinWidthRad;
        for (std::size_t i = 0; i < histogram.energyCentersGeV.size(); ++i) {
            for (std::size_t j = 0; j < histogram.thetaCentersRad.size(); ++j) {
                mean += histogram.energyCentersGeV[i]
                        * histogram.densityPerGeVRad[jointHistogramIndex(histogram, i, j)]
                        * cellArea;
            }
        }
        return mean;
    }

    inline double jointHistogramMeanThetaRad(const JointHistogram& histogram) {
        double mean           = 0.0;
        const double cellArea = histogram.energyBinWidthGeV * histogram.thetaBinWidthRad;
        for (std::size_t i = 0; i < histogram.energyCentersGeV.size(); ++i) {
            for (std::size_t j = 0; j < histogram.thetaCentersRad.size(); ++j) {
                mean += histogram.thetaCentersRad[j]
                        * histogram.densityPerGeVRad[jointHistogramIndex(histogram, i, j)]
                        * cellArea;
            }
        }
        return mean;
    }

    inline double jointHistogramL1Distance(const JointHistogram& lhs, const JointHistogram& rhs) {
        if (lhs.energyCentersGeV.size() != rhs.energyCentersGeV.size()
            || lhs.thetaCentersRad.size() != rhs.thetaCentersRad.size()) {
            throw std::runtime_error("LinearComptonBenchmark: joint histogram sizes do not match.");
        }

        double distance       = 0.0;
        const double cellArea = lhs.energyBinWidthGeV * lhs.thetaBinWidthRad;
        for (std::size_t i = 0; i < lhs.densityPerGeVRad.size(); ++i) {
            distance += std::abs(lhs.densityPerGeVRad[i] - rhs.densityPerGeVRad[i]) * cellArea;
        }
        return distance;
    }

}  // namespace LinearComptonBenchmark

#endif
