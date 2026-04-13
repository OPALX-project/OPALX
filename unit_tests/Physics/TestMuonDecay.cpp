#include "Physics/MuonDecay.h"
#include "Physics/Physics.h"
#include "gtest/gtest.h"

#include <cmath>

namespace {

TEST(MuonDecayTest, MaxElectronEnergyIsHalfMuonMass) {
    EXPECT_DOUBLE_EQ(Physics::MuonDecay::maxElectronEnergy(), Physics::m_mu / 2.0);
}

TEST(MuonDecayTest, MinElectronXIsPositiveAndSmall) {
    const double xMin = Physics::MuonDecay::minElectronX();
    EXPECT_GT(xMin, 0.0);
    EXPECT_LT(xMin, 0.02);  // m_e / (m_mu/2) ~ 0.0097
}

TEST(MuonDecayTest, MichelSpectrumIsNonNegativeInRange) {
    const double xMin = Physics::MuonDecay::minElectronX();
    constexpr int nSamples = 1000;
    for (int i = 0; i <= nSamples; ++i) {
        const double x = xMin + (1.0 - xMin) * static_cast<double>(i) / nSamples;
        EXPECT_GE(Physics::MuonDecay::michelSpectrum(x), 0.0);
    }
}

TEST(MuonDecayTest, MichelSpectrumPeaksAtXEqualsOne) {
    const double fPeak = Physics::MuonDecay::michelSpectrum(1.0);
    EXPECT_NEAR(fPeak, Physics::MuonDecay::michelUpperBound(), 1.0e-14);
}

TEST(MuonDecayTest, MichelSpectrumBelowUpperBound) {
    const double xMin = Physics::MuonDecay::minElectronX();
    const double upper = Physics::MuonDecay::michelUpperBound();
    constexpr int nSamples = 10000;
    for (int i = 0; i <= nSamples; ++i) {
        const double x = xMin + (1.0 - xMin) * static_cast<double>(i) / nSamples;
        EXPECT_LE(Physics::MuonDecay::michelSpectrum(x), upper + 1.0e-14);
    }
}

TEST(MuonDecayTest, MichelSpectrumEndpoints) {
    // f(0) = 0
    EXPECT_DOUBLE_EQ(Physics::MuonDecay::michelSpectrum(0.0), 0.0);
    // f(1) = 2*1*1*(3-2) = 2
    EXPECT_DOUBLE_EQ(Physics::MuonDecay::michelSpectrum(1.0), 2.0);
}

TEST(MuonDecayTest, MichelUpperBoundIs2) {
    EXPECT_DOUBLE_EQ(Physics::MuonDecay::michelUpperBound(), 2.0);
}

}  // namespace
