#include "LinearComptonSpectrumCommon.h"
#include "Utilities/Options.h"
#include "gtest/gtest.h"

namespace {
std::filesystem::path referenceSpectrumPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR)
        / "data"
        / "cain_linear_compton_90deg_xi029_histogram.csv";
}
}

TEST(TestLinearComptonSpectrum, WeakFieldSpectrumMatchesCainReference) {
    LinearComptonBenchmark::SpectrumConfig config;
    const auto opalxSpectrum = LinearComptonBenchmark::integrateLabSpectrum(config);
    const auto cainSpectrum = LinearComptonBenchmark::readSpectrumCSV(referenceSpectrumPath());

    ASSERT_EQ(opalxSpectrum.centersGeV.size(), cainSpectrum.centersGeV.size());
    EXPECT_NEAR(LinearComptonBenchmark::histogramArea(opalxSpectrum), 1.0, 5.0e-4);
    EXPECT_NEAR(LinearComptonBenchmark::histogramArea(cainSpectrum), 1.0, 5.0e-4);

    const double opalxMean = LinearComptonBenchmark::histogramMeanEnergyGeV(opalxSpectrum);
    const double cainMean = LinearComptonBenchmark::histogramMeanEnergyGeV(cainSpectrum);
    EXPECT_NEAR(opalxMean, cainMean, cainMean * 3.0e-2);

    const double l1Distance = LinearComptonBenchmark::histogramL1Distance(opalxSpectrum,
                                                                          cainSpectrum);
    EXPECT_LT(l1Distance, 0.1);
}


TEST(TestLinearComptonSpectrum, WeakFieldSampledSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    LinearComptonBenchmark::SpectrumConfig config;
    const auto opalxSpectrum = LinearComptonBenchmark::sampleLabSpectrum(config, 250000);
    const auto cainSpectrum = LinearComptonBenchmark::readSpectrumCSV(referenceSpectrumPath());

    ASSERT_EQ(opalxSpectrum.centersGeV.size(), cainSpectrum.centersGeV.size());
    EXPECT_NEAR(LinearComptonBenchmark::histogramArea(opalxSpectrum), 1.0, 1.5e-2);

    const double opalxMean = LinearComptonBenchmark::histogramMeanEnergyGeV(opalxSpectrum);
    const double cainMean = LinearComptonBenchmark::histogramMeanEnergyGeV(cainSpectrum);
    EXPECT_NEAR(opalxMean, cainMean, cainMean * 3.5e-2);

    const double l1Distance = LinearComptonBenchmark::histogramL1Distance(opalxSpectrum,
                                                                          cainSpectrum);
    EXPECT_LT(l1Distance, 0.14);

    Options::seed = previousSeed;
}
