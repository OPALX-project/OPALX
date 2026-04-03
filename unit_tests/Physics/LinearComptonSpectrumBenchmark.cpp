#include "LinearComptonSpectrumCommon.h"

#include <filesystem>
#include <iostream>

int main(int argc, char** argv) {
    std::filesystem::path output = argc >= 2
        ? std::filesystem::path(argv[1])
        : std::filesystem::path("opalx-linear-compton-90deg-xi029-spectrum.csv");

    LinearComptonBenchmark::SpectrumConfig config;
    const auto histogram = LinearComptonBenchmark::integrateLabSpectrum(config);
    LinearComptonBenchmark::writeSpectrumCSV(histogram, output);

    std::cout << "Wrote OPALX linear-Compton spectrum to " << output << '\n'
              << "Area = " << LinearComptonBenchmark::histogramArea(histogram) << '\n'
              << "Mean energy [GeV] = "
              << LinearComptonBenchmark::histogramMeanEnergyGeV(histogram) << '\n';
    return 0;
}
