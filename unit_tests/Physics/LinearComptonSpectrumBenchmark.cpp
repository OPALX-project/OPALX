#include "LinearComptonSpectrumCommon.h"
#include "Utilities/Options.h"

#include <filesystem>
#include <iostream>
#include <string>

namespace {
struct CliOptions {
    std::filesystem::path output = "opalx-linear-compton-90deg-xi029-spectrum.csv";
    bool sampled = false;
    bool angular = false;
    std::size_t samples = 250000;
    int seed = 13579;
};

CliOptions parseArguments(int argc, char** argv) {
    CliOptions options;

    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--sampled") {
            options.sampled = true;
        } else if (arg == "--angular") {
            options.angular = true;
        } else if (arg == "--samples") {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value after --samples");
            }
            options.samples = static_cast<std::size_t>(std::stoull(argv[++i]));
        } else if (arg == "--seed") {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value after --seed");
            }
            options.seed = std::stoi(argv[++i]);
        } else if (arg == "-h" || arg == "--help") {
            std::cout
                << "Usage: LinearComptonSpectrumBenchmark [output.csv] [--angular] [--sampled] [--samples N] [--seed S]\n"
                << "  default mode : deterministic energy spectrum\n"
                << "  --angular    : write the lab polar-angle histogram instead of the energy histogram\n"
                << "  --sampled    : host-only Monte Carlo benchmark using LinearCompton::sampleEvent\n"
                << "  --samples N  : number of Monte Carlo samples in sampled mode\n"
                << "  --seed S     : Options::seed value used in sampled mode\n";
            std::exit(0);
        } else if (!arg.empty() && arg[0] == '-') {
            throw std::runtime_error("Unknown option: " + arg);
        } else {
            options.output = arg;
        }
    }

    return options;
}
}  // namespace

int main(int argc, char** argv) {
    const CliOptions options = parseArguments(argc, argv);

    if (options.angular) {
        LinearComptonBenchmark::AngleConfig config;
        LinearComptonBenchmark::AngleHistogram histogram;
        if (options.sampled) {
            const int previousSeed = Options::seed;
            Options::seed = options.seed;
            histogram = LinearComptonBenchmark::sampleLabAngularSpectrum(config, options.samples);
            Options::seed = previousSeed;
        } else {
            histogram = LinearComptonBenchmark::integrateLabAngularSpectrum(config);
        }

        LinearComptonBenchmark::writeAngleCSV(histogram, options.output);

        std::cout << "Wrote OPALX linear-Compton angular spectrum to " << options.output << '\n'
                  << "Mode = " << (options.sampled ? "sampled" : "deterministic") << '\n'
                  << "Observable = theta_lab [rad]\n"
                  << "Area = " << LinearComptonBenchmark::angleHistogramArea(histogram) << '\n'
                  << "Mean angle [rad] = "
                  << LinearComptonBenchmark::angleHistogramMeanRad(histogram) << '\n';
        if (options.sampled) {
            std::cout << "Samples = " << options.samples << '\n'
                      << "Seed = " << options.seed << '\n';
        }
        return 0;
    }

    LinearComptonBenchmark::SpectrumConfig config;
    LinearComptonBenchmark::SpectrumHistogram histogram;
    if (options.sampled) {
        const int previousSeed = Options::seed;
        Options::seed = options.seed;
        histogram = LinearComptonBenchmark::sampleLabSpectrum(config, options.samples);
        Options::seed = previousSeed;
    } else {
        histogram = LinearComptonBenchmark::integrateLabSpectrum(config);
    }

    LinearComptonBenchmark::writeSpectrumCSV(histogram, options.output);

    std::cout << "Wrote OPALX linear-Compton spectrum to " << options.output << '\n'
              << "Mode = " << (options.sampled ? "sampled" : "deterministic") << '\n'
              << "Observable = E_gamma [GeV]\n"
              << "Area = " << LinearComptonBenchmark::histogramArea(histogram) << '\n'
              << "Mean energy [GeV] = "
              << LinearComptonBenchmark::histogramMeanEnergyGeV(histogram) << '\n';
    if (options.sampled) {
        std::cout << "Samples = " << options.samples << '\n'
                  << "Seed = " << options.seed << '\n';
    }
    return 0;
}
