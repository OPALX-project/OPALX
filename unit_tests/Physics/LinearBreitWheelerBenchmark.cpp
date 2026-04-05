#include "LinearBreitWheelerBenchmarkCommon.h"
#include "Utilities/Options.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
/**
 * @brief Command-line options for the standalone linear Breit-Wheeler benchmark driver.
 *
 * The executable supports three validation modes:
 *
 * - raw sampled one-dimensional observables for the fixed-geometry benchmark,
 * - raw sampled one-dimensional observables for a finite incoming photon beam,
 * - sampled joint `E`-`theta` histograms for two-dimensional CAIN comparisons.
 */
struct CliOptions {
    std::filesystem::path output = "opalx-linear-breit-wheeler-spectrum.csv";
    std::string particle = "electron";
    std::string observable = "energy";
    bool joint = false;
    bool finitePhotonBeam = false;
    std::size_t samples = 250000;
    int seed = 13579;
    double wavelength_m = 1.0e-9;
    double highEnergyPhotonEnergyGeV = 0.5;
    double photonSigmaThetaRad = 1.0e-3;
    double photonRelativeEnergySpread = 0.0;
    std::size_t bins = 80;
    double minValue = 0.0;
    double maxValue = 0.5;
    double thetaMaxRad = 0.0045;
};

CliOptions parseArguments(int argc, char** argv) {
    CliOptions options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--particle") {
            options.particle = argv[++i];
        } else if (arg == "--observable") {
            options.observable = argv[++i];
        } else if (arg == "--joint") {
            options.joint = true;
        } else if (arg == "--finite-photon-beam") {
            options.finitePhotonBeam = true;
        } else if (arg == "--samples") {
            options.samples = static_cast<std::size_t>(std::stoull(argv[++i]));
        } else if (arg == "--seed") {
            options.seed = std::stoi(argv[++i]);
        } else if (arg == "--wavelength") {
            options.wavelength_m = std::stod(argv[++i]);
        } else if (arg == "--high-energy-photon") {
            options.highEnergyPhotonEnergyGeV = std::stod(argv[++i]);
        } else if (arg == "--photon-sigma-theta") {
            options.photonSigmaThetaRad = std::stod(argv[++i]);
        } else if (arg == "--photon-relative-energy-spread") {
            options.photonRelativeEnergySpread = std::stod(argv[++i]);
        } else if (arg == "--bins") {
            options.bins = static_cast<std::size_t>(std::stoull(argv[++i]));
        } else if (arg == "--min") {
            options.minValue = std::stod(argv[++i]);
        } else if (arg == "--max") {
            options.maxValue = std::stod(argv[++i]);
        } else if (arg == "--theta-max") {
            options.thetaMaxRad = std::stod(argv[++i]);
        } else if (arg == "-h" || arg == "--help") {
            std::cout
                << "Usage: LinearBreitWheelerBenchmark [output.csv] [--particle electron|positron]"
                << " [--observable energy|theta] [--joint] [--finite-photon-beam] [--samples N] [--seed S] [--bins N]"
                << " [--high-energy-photon E] [--photon-sigma-theta S] [--photon-relative-energy-spread S]"
                << " [--min V] [--max V] [--theta-max V]\n";
            std::exit(0);
        } else if (!arg.empty() && arg[0] == '-') {
            throw std::runtime_error("Unknown option: " + arg);
        } else {
            options.output = arg;
        }
    }
    return options;
}

Vector_t<double, 3> makeDirection(double x, double y, double z) {
    Vector_t<double, 3> value(0.0);
    value(0) = x;
    value(1) = y;
    value(2) = z;
    return value;
}

LinearBreitWheelerBenchmark::FinalState parseFinalState(const std::string& particle) {
    if (particle == "positron") {
        return LinearBreitWheelerBenchmark::FinalState::Positron;
    }
    return LinearBreitWheelerBenchmark::FinalState::Electron;
}

LinearBreitWheelerBenchmark::Observable parseObservable(const std::string& observable) {
    if (observable == "theta") {
        return LinearBreitWheelerBenchmark::Observable::Theta;
    }
    return LinearBreitWheelerBenchmark::Observable::Energy;
}
}  // namespace

int main(int argc, char** argv) {
    const CliOptions options = parseArguments(argc, argv);
    const int previousSeed = Options::seed;
    Options::seed = options.seed;

    const auto state = parseFinalState(options.particle);

    if (options.joint) {
        if (options.finitePhotonBeam) {
            LinearBreitWheelerBenchmark::FinitePhotonBeamJointConfig config;
            config.centralHighEnergyPhotonEnergyGeV = options.highEnergyPhotonEnergyGeV;
            config.wavelength_m = options.wavelength_m;
            config.referenceHighEnergyDirection = makeDirection(0.0, 0.0, 1.0);
            config.laserDirection = makeDirection(0.0, 0.0, -1.0);
            config.sigmaThetaXRad = options.photonSigmaThetaRad;
            config.sigmaThetaYRad = options.photonSigmaThetaRad;
            config.relativeEnergySpread = options.photonRelativeEnergySpread;
            config.energyBins = options.bins;
            config.thetaBins = options.bins;
            config.energyMinGeV = options.minValue;
            config.energyMaxGeV = options.maxValue;
            config.thetaMinRad = 0.0;
            config.thetaMaxRad = options.thetaMaxRad;

            const auto histogram = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamJointHistogram(config,
                                                                                                      state,
                                                                                                      options.samples);
            LinearBreitWheelerBenchmark::writeJointHistogramCSV(histogram, options.output);
        } else {
            LinearBreitWheelerBenchmark::JointHistogramConfig config;
            config.highEnergyPhotonEnergyGeV = options.highEnergyPhotonEnergyGeV;
            config.wavelength_m = options.wavelength_m;
            config.highEnergyDirection = makeDirection(0.0, 0.0, 1.0);
            config.laserDirection = makeDirection(0.0, 0.0, -1.0);
            config.energyBins = options.bins;
            config.thetaBins = options.bins;
            config.energyMinGeV = options.minValue;
            config.energyMaxGeV = options.maxValue;
            config.thetaMinRad = 0.0;
            config.thetaMaxRad = options.thetaMaxRad;

            const auto histogram = LinearBreitWheelerBenchmark::sampleJointHistogram(config,
                                                                                     state,
                                                                                     options.samples);
            LinearBreitWheelerBenchmark::writeJointHistogramCSV(histogram, options.output);
        }
        Options::seed = previousSeed;
        std::cout << "Wrote OPALX linear Breit-Wheeler joint benchmark to " << options.output << '\n';
        return 0;
    }

    const auto observable = parseObservable(options.observable);
    auto engine = Physics::LinearBreitWheeler::makeHostRandomEngine();

    std::ofstream output(options.output);
    output << "# value\n";

    if (options.finitePhotonBeam) {
        const double laserPhotonEnergyGeV = Physics::LinearBreitWheeler::photonEnergyFromWavelengthGeV(options.wavelength_m);
        LinearBreitWheelerBenchmark::FinitePhotonBeamConfig config;
        config.centralHighEnergyPhotonEnergyGeV = options.highEnergyPhotonEnergyGeV;
        config.wavelength_m = options.wavelength_m;
        config.referenceHighEnergyDirection = makeDirection(0.0, 0.0, 1.0);
        config.laserDirection = makeDirection(0.0, 0.0, -1.0);
        config.sigmaThetaXRad = options.photonSigmaThetaRad;
        config.sigmaThetaYRad = options.photonSigmaThetaRad;
        config.relativeEnergySpread = options.photonRelativeEnergySpread;

        for (std::size_t i = 0; i < options.samples; ++i) {
            const double highEnergyPhotonEnergyGeV = LinearBreitWheelerBenchmark::sampleHighEnergyPhotonEnergyGeV(
                config.centralHighEnergyPhotonEnergyGeV,
                config.relativeEnergySpread,
                engine);
            const auto highEnergyDirection = LinearBreitWheelerBenchmark::samplePhotonBeamDirection(
                config.referenceHighEnergyDirection,
                config.sigmaThetaXRad,
                config.sigmaThetaYRad,
                engine);
            const auto kernel = Physics::LinearBreitWheeler::makeSamplingKernel(highEnergyPhotonEnergyGeV,
                                                                                laserPhotonEnergyGeV,
                                                                                highEnergyDirection,
                                                                                config.laserDirection);
            const auto event = Physics::LinearBreitWheeler::sampleEvent(kernel, engine);
            const double value = LinearBreitWheelerBenchmark::sampledObservable(event, state, observable);
            output << value << '\n';
        }
    } else {
        const double laserPhotonEnergyGeV = Physics::LinearBreitWheeler::photonEnergyFromWavelengthGeV(options.wavelength_m);
        const auto kernel = Physics::LinearBreitWheeler::makeSamplingKernel(options.highEnergyPhotonEnergyGeV,
                                                                            laserPhotonEnergyGeV,
                                                                            makeDirection(0.0, 0.0, 1.0),
                                                                            makeDirection(0.0, 0.0, -1.0));
        for (std::size_t i = 0; i < options.samples; ++i) {
            const auto event = Physics::LinearBreitWheeler::sampleEvent(kernel, engine);
            const double value = LinearBreitWheelerBenchmark::sampledObservable(event, state, observable);
            output << value << '\n';
        }
    }

    Options::seed = previousSeed;
    std::cout << "Wrote OPALX linear Breit-Wheeler sampled benchmark to " << options.output << '\n';
    return 0;
}
