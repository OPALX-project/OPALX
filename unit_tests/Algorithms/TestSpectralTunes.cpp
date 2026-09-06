#include <cmath>
#include <limits>
#include "Algorithms/SpectralTunes.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

TEST(SpectralTunes, LegacyGridAndOffsetInvariance) {
    const size_t n       = 1440;
    const unsigned turns = 100;
    std::vector<double> v(n), offset(n);
    for (size_t i = 0; i < n; ++i) {
        v[i]      = std::cos(2 * std::acos(-1.0) * 1.155 * turns * i / (n - 1) + 0.31);
        offset[i] = 7 + 0.005 * v[i];
    }
    auto a = SpectralTunes::analyze(v, turns), b = SpectralTunes::analyze(offset, turns);
    EXPECT_DOUBLE_EQ(a.peak.tune, 1.155);
    EXPECT_DOUBLE_EQ(a.peak.spacing, 0.0025);
    EXPECT_EQ(a.tune.size(), 2304u);
    EXPECT_DOUBLE_EQ(a.peak.tune, b.peak.tune);
    EXPECT_NEAR(a.peak.power, b.peak.power, 1e-8);
}
TEST(SpectralTunes, VerticalPeakAndDominantMode) {
    std::vector<double> v(1440);
    for (size_t i = 0; i < v.size(); ++i) {
        const double t = 100.0 * i / (v.size() - 1);
        v[i]           = std::sin(2 * std::acos(-1.0) * 0.89 * t)
               + 0.1 * std::cos(2 * std::acos(-1.0) * 1.155 * t);
    }
    EXPECT_DOUBLE_EQ(SpectralTunes::analyze(v, 100).peak.tune, 0.89);
}
TEST(SpectralTunes, RejectInvalidSignals) {
    EXPECT_THROW(SpectralTunes::analyze(std::vector<double>(5, 1), 100), OpalException);
    EXPECT_THROW(SpectralTunes::analyze(std::vector<double>(20, 1), 100), OpalException);
    EXPECT_THROW(SpectralTunes::analyze(std::vector<double>(20, 1), 0), OpalException);
    EXPECT_THROW(
            SpectralTunes::analyze(
                    std::vector<double>(20, std::numeric_limits<double>::quiet_NaN()), 100),
            OpalException);
}
