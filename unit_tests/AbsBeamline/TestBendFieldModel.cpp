//
// Tests for BendFieldModel: the stateless Enge fringe-field shape math shared by
// the analytic SBEND and RBEND bends.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
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

#include "AbsBeamline/BendFieldModel.h"
#include "gtest/gtest.h"

#include <cmath>

namespace {
    constexpr double gap = 0.02;  // full gap used throughout, metres
}

// F(0) ~ 0.3825 at the pole face; F -> 1 deep inside; F -> 0 well outside.
TEST(BendFieldModel, EngeProfileValueAtFaceAndLimits) {
    EXPECT_NEAR(BendFieldModel::engeProfile(0.0, gap).value, 0.38253, 1.0e-4);
    EXPECT_NEAR(BendFieldModel::engeProfile(-10.0 * gap, gap).value, 1.0, 1.0e-6);
    EXPECT_NEAR(BendFieldModel::engeProfile(10.0 * gap, gap).value, 0.0, 1.0e-6);
    // Monotonically decreasing from inside to outside.
    EXPECT_GT(
            BendFieldModel::engeProfile(-gap, gap).value,
            BendFieldModel::engeProfile(gap, gap).value);
}

// A non-positive gap is a hard edge: F = 1, no derivatives.
TEST(BendFieldModel, EngeProfileHardEdgeForZeroGap) {
    const BendFieldModel::EngeValue f = BendFieldModel::engeProfile(0.123, 0.0);
    EXPECT_DOUBLE_EQ(f.value, 1.0);
    EXPECT_DOUBLE_EQ(f.firstDerivative, 0.0);
    EXPECT_DOUBLE_EQ(f.secondDerivative, 0.0);
}

// The analytic derivatives match finite differences of the value.
TEST(BendFieldModel, EngeProfileDerivativesMatchFiniteDifferences) {
    for (const double d : {-0.03, -0.01, 0.0, 0.01, 0.03}) {
        const BendFieldModel::EngeValue f = BendFieldModel::engeProfile(d, gap);

        const double h  = 1.0e-6;
        const double fp = 1.0e-4;
        const double fdFirst =
                (BendFieldModel::engeProfile(d + h, gap).value
                 - BendFieldModel::engeProfile(d - h, gap).value)
                / (2.0 * h);
        const double fdSecond =
                (BendFieldModel::engeProfile(d + fp, gap).value
                 - 2.0 * f.value + BendFieldModel::engeProfile(d - fp, gap).value)
                / (fp * fp);

        EXPECT_NEAR(f.firstDerivative, fdFirst, 1.0e-3 * (1.0 + std::abs(f.firstDerivative)));
        EXPECT_NEAR(f.secondDerivative, fdSecond, 1.0e-2 * (1.0 + std::abs(f.secondDerivative)));
    }
}

// Gap selection and the fringe half width.
TEST(BendFieldModel, GapAndHalfWidthHelpers) {
    EXPECT_DOUBLE_EQ(BendFieldModel::profileGap(0.04, 0.0), 0.04);   // GAP wins
    EXPECT_DOUBLE_EQ(BendFieldModel::profileGap(0.0, 0.02), 0.04);   // fall back to 2*HGAP
    EXPECT_DOUBLE_EQ(BendFieldModel::fringeHalfWidth(0.02), 0.1);    // five gaps
    EXPECT_DOUBLE_EQ(BendFieldModel::safeAbsCos(0.0), 1.0);
    EXPECT_NEAR(BendFieldModel::safeAbsCos(M_PI / 2.0), 1.0e-6, 1.0e-9);  // floored
}

// The combined amplitude is 1 deep inside, ~0.3825 at a face, 0 outside, and
// reports which face limits it.
TEST(BendFieldModel, FringeAmplitudeCombinesTwoFaces) {
    // Deep interior: both distances well inside -> full field.
    const auto interior = BendFieldModel::fringeAmplitude(-10.0 * gap, -10.0 * gap, gap);
    EXPECT_NEAR(interior.value, 1.0, 1.0e-6);

    // At the entrance face, exit far inside -> entrance limits.
    const auto atEntrance = BendFieldModel::fringeAmplitude(0.0, -10.0 * gap, gap);
    EXPECT_NEAR(atEntrance.value, 0.38253, 1.0e-4);
    EXPECT_EQ(atEntrance.activeFace, 0);

    // At the exit face, entrance far inside -> exit limits.
    const auto atExit = BendFieldModel::fringeAmplitude(-10.0 * gap, 0.0, gap);
    EXPECT_NEAR(atExit.value, 0.38253, 1.0e-4);
    EXPECT_EQ(atExit.activeFace, 1);

    // Beyond the entrance fringe -> zero.
    const auto outside = BendFieldModel::fringeAmplitude(10.0 * gap, -10.0 * gap, gap);
    EXPECT_NEAR(outside.value, 0.0, 1.0e-6);
}

// Edge kick: with FINT = 0 the psi correction vanishes, leaving -h tan(E).
TEST(BendFieldModel, EdgeKickReducesToTanWithoutFint) {
    const double h = 0.5, E = 0.2;
    EXPECT_NEAR(
            BendFieldModel::edgeVerticalKickCoefficient(h, 0.02, 0.0, E), -h * std::tan(E), 1.0e-12);

    // A finite FINT shifts the effective edge angle by psi = h*HGAP*FINT*(1+sin^2E)/cosE.
    const double hgap = 0.02, fint = 0.5;
    const double psi  = h * hgap * fint * (1.0 + std::sin(E) * std::sin(E)) / std::cos(E);
    EXPECT_NEAR(
            BendFieldModel::edgeVerticalKickCoefficient(h, hgap, fint, E), -h * std::tan(E - psi),
            1.0e-12);
}
