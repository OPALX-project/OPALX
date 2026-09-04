//
// Tests for MultipoleTFieldModel: the stateless field math of the combined-function
// multipole with a tanh fringe, straight and curved.
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

#include "AbsBeamline/MultipoleTFieldModel.h"
#include "gtest/gtest.h"

#include <cmath>
#include <vector>

namespace mtfm = MultipoleTFieldModel;

namespace {
    /// A plain host table with the two members the model templates need, so these
    /// tests need no Kokkos runtime (see TestBendFieldModel).
    class Table {
    public:
        Table(const unsigned int rows, const unsigned int columns)
            : rows_m(rows), columns_m(columns), data_m(rows * columns, 0.0) {}

        double& operator()(const unsigned int row, const unsigned int column) {
            return data_m[row * columns_m + column];
        }
        double operator()(const unsigned int row, const unsigned int column) const {
            return data_m[row * columns_m + column];
        }
        unsigned int extent(const unsigned int dimension) const {
            return (dimension == 0) ? rows_m : columns_m;
        }

    private:
        unsigned int rows_m;
        unsigned int columns_m;
        std::vector<double> data_m;
    };

    /// The table sized the way the element sizes it: 2 * maxFOrder + 1 derivatives.
    Table makeTable(const unsigned int maxFOrder) {
        const unsigned int numDerivatives = 2 * maxFOrder + 1;
        Table table(numDerivatives + 1, numDerivatives + 2);
        mtfm::tanhDerivativeTable(numDerivatives, table);
        return table;
    }

    mtfm::FieldInputs makeInputs(
            const std::vector<double>& profile, const double length, const double curvature,
            const double fringe, const unsigned int maxFOrder = 5) {
        mtfm::FieldInputs in{};
        for (unsigned int i = 0; i < mtfm::NumPoles; ++i) {
            in.profile[i] = (i < profile.size()) ? profile[i] : 0.0;
        }
        in.maxFOrder   = maxFOrder;
        in.bodyLength  = length;
        in.curvature   = curvature;
        in.fringeLeft  = fringe;
        in.fringeRight = fringe;
        return in;
    }

    /// (radial offset, vertical, arc length) as GeometryHelper::toBendArcCoords gives it.
    Vector_t<double, 3> arc(const double x, const double y, const double s) {
        return Vector_t<double, 3>(x, y, s);
    }
}  // namespace

// ---------------------------------------------------------------------------
// The tanh derivative table
// ---------------------------------------------------------------------------

// Row n holds the polynomial P_n with d^n/du^n tanh u = P_n(tanh u), from
// P_0(t) = t and P_n = (1 - t^2) P_{n-1}'.
TEST(MultipoleTFieldModel, TanhDerivativeTableRows) {
    const unsigned int numDerivatives = 5;
    Table table(numDerivatives + 1, numDerivatives + 2);
    mtfm::tanhDerivativeTable(numDerivatives, table);

    const std::vector<std::vector<double>> expected = {
            {0, 1, 0, 0, 0, 0, 0},
            {1, 0, -1, 0, 0, 0, 0},
            {0, -2, 0, 2, 0, 0, 0},
            {-2, 0, 8, 0, -6, 0, 0},
            {0, 16, 0, -40, 0, 24, 0},
            {16, 0, -136, 0, 240, 0, -120}};
    for (unsigned int n = 0; n < expected.size(); ++n) {
        for (unsigned int k = 0; k < expected[n].size(); ++k) {
            EXPECT_DOUBLE_EQ(table(n, k), expected[n][k]) << "row " << n << " column " << k;
        }
    }
}

// ---------------------------------------------------------------------------
// The transverse profile
// ---------------------------------------------------------------------------

// T(x) = sum_k b_k x^k and its derivatives, at x = 1 for b = {0,1,2,3,4,5}.
TEST(MultipoleTFieldModel, TransverseDerivatives) {
    Kokkos::Array<double, mtfm::NumPoles> profile{{0.0, 1.0, 2.0, 3.0, 4.0, 5.0}};
    Kokkos::Array<double, mtfm::MaxDerivatives> out{};
    mtfm::transverseDerivatives(profile, 7, 1.0, out);

    const std::vector<double> expected = {15.0, 55.0, 170.0, 414.0, 696.0, 600.0, 0.0};
    for (unsigned int i = 0; i < expected.size(); ++i) {
        EXPECT_DOUBLE_EQ(out[i], expected[i]) << "derivative " << i;
    }
}

// A pure dipole profile is constant in x with no derivatives.
TEST(MultipoleTFieldModel, TransverseDerivativesOfADipole) {
    Kokkos::Array<double, mtfm::NumPoles> profile{{0.7, 0.0, 0.0, 0.0, 0.0, 0.0}};
    Kokkos::Array<double, mtfm::MaxDerivatives> out{};
    mtfm::transverseDerivatives(profile, 4, 0.123, out);

    EXPECT_DOUBLE_EQ(out[0], 0.7);
    EXPECT_DOUBLE_EQ(out[1], 0.0);
    EXPECT_DOUBLE_EQ(out[2], 0.0);
    EXPECT_DOUBLE_EQ(out[3], 0.0);
}

// ---------------------------------------------------------------------------
// The tanh fringe
// ---------------------------------------------------------------------------

// S is symmetric about the centre with equal fringe lengths, so even derivatives
// are equal at +/- s and odd ones change sign.
TEST(MultipoleTFieldModel, FringeDerivativeParity) {
    Table table = makeTable(3);
    Kokkos::Array<double, mtfm::MaxDerivatives> plus{};
    Kokkos::Array<double, mtfm::MaxDerivatives> minus{};
    mtfm::fringeDerivatives(2.0, 1.0, 1.0, +2.0, table, plus);
    mtfm::fringeDerivatives(2.0, 1.0, 1.0, -2.0, table, minus);

    for (unsigned int i = 0; i < table.extent(0); ++i) {
        if (i % 2 == 0) {
            EXPECT_NEAR(plus[i], minus[i], 1.0e-10) << "even derivative " << i;
        } else {
            EXPECT_NEAR(plus[i], -minus[i], 1.0e-10) << "odd derivative " << i;
        }
    }
}

// Reference values of S and its derivatives at the exit face of a 4 m body with
// 1 m fringes (the numbers the previous MultipoleTBase test pinned).
TEST(MultipoleTFieldModel, FringeDerivativeValues) {
    Table table = makeTable(3);
    Kokkos::Array<double, mtfm::MaxDerivatives> out{};
    mtfm::fringeDerivatives(2.0, 1.0, 1.0, 2.0, table, out);

    const std::vector<double> expected = {
            0.49966464986953352,
            -0.49932952465848707,
            -0.0013400513070529474,
            1.0026765069198489,
            -0.0053386419156264964,
            -7.9893801387861263,
            -0.021010422121266359};
    for (unsigned int i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(out[i], expected[i], 1.0e-10) << "derivative " << i;
    }
}

// A zero fringe length on both sides is a hard edge: S = 1 inside the body, 0
// outside, no derivatives, and no division by zero.
TEST(MultipoleTFieldModel, HardEdgeFringe) {
    Table table = makeTable(3);
    Kokkos::Array<double, mtfm::MaxDerivatives> inside{};
    Kokkos::Array<double, mtfm::MaxDerivatives> outside{};
    mtfm::fringeDerivatives(2.0, 0.0, 0.0, 0.0, table, inside);
    mtfm::fringeDerivatives(2.0, 0.0, 0.0, 3.0, table, outside);

    EXPECT_DOUBLE_EQ(inside[0], 1.0);
    EXPECT_DOUBLE_EQ(outside[0], 0.0);
    for (unsigned int i = 1; i < table.extent(0); ++i) {
        EXPECT_DOUBLE_EQ(inside[i], 0.0) << "derivative " << i;
        EXPECT_DOUBLE_EQ(outside[i], 0.0) << "derivative " << i;
        EXPECT_FALSE(std::isnan(inside[i]));
    }
}

// One hard side and one soft side: the step acts on its own side only.
TEST(MultipoleTFieldModel, HalfHardEdgeFringe) {
    Table table = makeTable(3);
    Kokkos::Array<double, mtfm::MaxDerivatives> out{};
    // Hard entrance, soft exit: at the exit face S is the tanh half value.
    mtfm::fringeDerivatives(2.0, 0.0, 0.1, 2.0, table, out);
    EXPECT_NEAR(out[0], 0.5, 1.0e-6);
    EXPECT_FALSE(std::isnan(out[1]));
    // Deep inside the body S is 1 from both sides.
    mtfm::fringeDerivatives(2.0, 0.0, 0.1, 0.0, table, out);
    EXPECT_NEAR(out[0], 1.0, 1.0e-6);
}

// ---------------------------------------------------------------------------
// The straight field
// ---------------------------------------------------------------------------

// A dipole profile gives the flat-top field in the body and half of it at the
// face, with nothing in the other components on the mid-plane.
TEST(MultipoleTFieldModel, StraightDipoleAlongTheAxis) {
    Table table          = makeTable(5);
    const double length  = 4.4;
    const mtfm::FieldInputs in = makeInputs({1.0}, length, 0.0, 0.3);

    const Vector_t<double, 3> centre = mtfm::straightField(0.0, 0.0, 0.0, in, table);
    EXPECT_NEAR(centre(1), 1.0, 1.0e-2);
    EXPECT_DOUBLE_EQ(centre(0), 0.0);
    EXPECT_DOUBLE_EQ(centre(2), 0.0);

    const Vector_t<double, 3> face = mtfm::straightField(0.0, 0.0, -0.5 * length, in, table);
    EXPECT_NEAR(face(1), 0.5, 1.0e-2);

    const Vector_t<double, 3> outside =
            mtfm::straightField(0.0, 0.0, -0.5 * length - 6.0 * 0.3, in, table);
    EXPECT_NEAR(outside(1), 0.0, 1.0e-4);
}

// A dipole is constant across the aperture; a gradient is linear in x, with the
// same sign convention as MULTIPOLE and SBEND (positive TP[1] raises By at +x).
TEST(MultipoleTFieldModel, StraightProfileAcrossTheAperture) {
    Table table                = makeTable(5);
    const double length        = 4.4;
    const mtfm::FieldInputs dipole = makeInputs({1.0}, length, 0.0, 0.3);
    const mtfm::FieldInputs quad   = makeInputs({0.0, 1.0}, length, 0.0, 0.3);

    for (const double x : {-1.0, -0.5, 0.0, 0.5, 1.0}) {
        EXPECT_NEAR(mtfm::straightField(x, 0.0, 0.0, dipole, table)(1), 1.0, 1.0e-2);
        EXPECT_NEAR(mtfm::straightField(x, 0.0, 0.0, quad, table)(1), x, 1.0e-2);
        // At the face the whole profile is scaled by S = 1/2.
        EXPECT_NEAR(
                mtfm::straightField(x, 0.0, -0.5 * length, quad, table)(1), 0.5 * x, 1.0e-2);
    }
}

// ---------------------------------------------------------------------------
// The curved field
// ---------------------------------------------------------------------------

// On the design orbit of a curved dipole the field is the flat-top value, purely
// vertical, whichever way the body bends.
TEST(MultipoleTFieldModel, CurvedDipoleOnTheDesignOrbit) {
    Table table         = makeTable(5);
    const double length = 1.0;
    for (const double curvature : {+0.5, -0.5}) {
        const mtfm::FieldInputs in = makeInputs({1.0}, length, curvature, 0.02);
        const Vector_t<double, 3> B = mtfm::curvedField(0.0, 0.0, 0.0, in, table);
        EXPECT_NEAR(B(1), 1.0, 1.0e-6) << "curvature " << curvature;
        EXPECT_DOUBLE_EQ(B(0), 0.0);
        EXPECT_DOUBLE_EQ(B(2), 0.0);
    }
}

// The Frenet x of the curved field continues the entrance x, so a positive
// gradient raises By at positive x for either bend direction.
TEST(MultipoleTFieldModel, CurvedGradientSignFollowsTheEntranceFrame) {
    Table table         = makeTable(5);
    const double length = 1.0;
    for (const double curvature : {+0.5, -0.5}) {
        const mtfm::FieldInputs in = makeInputs({0.0, 1.0}, length, curvature, 0.02);
        const double atPlus  = mtfm::curvedField(+0.01, 0.0, 0.0, in, table)(1);
        const double atMinus = mtfm::curvedField(-0.01, 0.0, 0.0, in, table)(1);
        EXPECT_GT(atPlus, 0.0) << "curvature " << curvature;
        EXPECT_LT(atMinus, 0.0) << "curvature " << curvature;
        EXPECT_NEAR(atPlus, 0.01, 1.0e-3);
    }
}

// A weakly curved body tends to the straight one.
TEST(MultipoleTFieldModel, CurvedTendsToStraightAtSmallCurvature) {
    Table table                  = makeTable(5);
    const double length          = 1.0;
    const mtfm::FieldInputs curved   = makeInputs({0.5, 0.3}, length, 1.0e-6, 0.05);
    const mtfm::FieldInputs straight = makeInputs({0.5, 0.3}, length, 0.0, 0.05);

    for (const double s : {-0.4, 0.0, 0.4}) {
        const Vector_t<double, 3> a = mtfm::curvedField(0.02, 0.01, s, curved, table);
        const Vector_t<double, 3> b = mtfm::straightField(0.02, 0.01, s, straight, table);
        for (unsigned int d = 0; d < 3; ++d) {
            EXPECT_NEAR(a(d), b(d), 1.0e-6) << "component " << d << " at s = " << s;
        }
    }
}

// ---------------------------------------------------------------------------
// Maxwell
// ---------------------------------------------------------------------------

namespace {
    /// Relative divergence and curl by central differences, in units of |B| / dr.
    void divergenceAndCurl(
            const mtfm::FieldInputs& in, const Table& table, const double x, const double y,
            const double s, const double dr, double& divergence, double& curl) {
        const auto value = [&](const double a, const double b, const double c) {
            return (in.curvature == 0.0) ? mtfm::straightField(a, b, c, in, table)
                                         : mtfm::curvedField(a, b, c, in, table);
        };
        const Vector_t<double, 3> xPlus  = value(x + dr, y, s);
        const Vector_t<double, 3> xMinus = value(x - dr, y, s);
        const Vector_t<double, 3> yPlus  = value(x, y + dr, s);
        const Vector_t<double, 3> yMinus = value(x, y - dr, s);
        const Vector_t<double, 3> sPlus  = value(x, y, s + dr);
        const Vector_t<double, 3> sMinus = value(x, y, s - dr);

        divergence = ((xPlus(0) - xMinus(0)) + (yPlus(1) - yMinus(1)) + (sPlus(2) - sMinus(2)))
                     / (2.0 * dr);
        const double curlX = ((yPlus(2) - yMinus(2)) - (sPlus(1) - sMinus(1))) / (2.0 * dr);
        const double curlY = ((sPlus(0) - sMinus(0)) - (xPlus(2) - xMinus(2))) / (2.0 * dr);
        const double curlZ = ((xPlus(1) - xMinus(1)) - (yPlus(0) - yMinus(0))) / (2.0 * dr);
        curl               = std::hypot(curlX, curlY, curlZ);

        const Vector_t<double, 3> centre = value(x, y, s);
        const double magnitude           = std::hypot(centre(0), centre(1), centre(2));
        if (magnitude > 0.0) {
            divergence /= magnitude;
            curl /= magnitude;
        }
    }
}  // namespace

// A field from a scalar potential is divergence and curl free. Checked along the
// body of a straight magnet with every pole populated, fringes included.
TEST(MultipoleTFieldModel, StraightFieldIsSourceFree) {
    Table table                = makeTable(5);
    const double length        = 4.0;
    const mtfm::FieldInputs in = makeInputs({1.0, 1.0, 1.0, 1.0, 1.0}, length, 0.0, 0.3);
    const double dr            = 1.0e-3;

    for (int i = 0; i <= 40; ++i) {
        const double s = -0.75 * length + 1.5 * length * i / 40.0;
        double divergence = 0.0;
        double curl       = 0.0;
        divergenceAndCurl(in, table, 0.05, 0.02, s, dr, divergence, curl);
        EXPECT_NEAR(divergence * dr, 0.0, 1.0e-5) << "at s = " << s;
        EXPECT_NEAR(curl * dr, 0.0, 1.0e-2) << "at s = " << s;
    }
}

// The same for a curved body, where the curl also has to cancel the Frenet scale
// factor. (Vector operators in the curved frame carry metric terms, so this is a
// looser check than the straight one; it still catches a wrong recursion.)
TEST(MultipoleTFieldModel, CurvedFieldIsSourceFree) {
    Table table                = makeTable(5);
    const double length        = 4.0;
    const mtfm::FieldInputs in = makeInputs({1.0, 1.0, 1.0, 1.0, 1.0}, length, 0.25, 0.3);
    const double dr            = 1.0e-3;

    for (int i = 0; i <= 20; ++i) {
        const double s = -0.5 * length + length * i / 20.0;
        double divergence = 0.0;
        double curl       = 0.0;
        divergenceAndCurl(in, table, 0.05, 0.02, s, dr, divergence, curl);
        EXPECT_NEAR(divergence * dr, 0.0, 1.0e-2) << "at s = " << s;
        EXPECT_NEAR(curl * dr, 0.0, 1.0e-2) << "at s = " << s;
    }
}

// ---------------------------------------------------------------------------
// The dispatcher
// ---------------------------------------------------------------------------

// field() takes arc coordinates (arc length from the entrance face) and picks the
// straight or curved recursion from the curvature.
TEST(MultipoleTFieldModel, FieldDispatchesOnCurvatureAndShiftsTheArcLength) {
    Table table         = makeTable(5);
    const double length = 1.0;

    const mtfm::FieldInputs straight = makeInputs({1.0}, length, 0.0, 0.02);
    // Arc length L/2 is the centre of the body.
    EXPECT_NEAR(mtfm::field(arc(0.0, 0.0, 0.5 * length), straight, table)(1), 1.0, 1.0e-6);
    // Arc length 0 is the entrance face.
    EXPECT_NEAR(mtfm::field(arc(0.0, 0.0, 0.0), straight, table)(1), 0.5, 1.0e-6);
    EXPECT_NEAR(mtfm::field(arc(0.0, 0.0, length), straight, table)(1), 0.5, 1.0e-6);

    const mtfm::FieldInputs curved = makeInputs({1.0}, length, 0.5, 0.02);
    EXPECT_NEAR(mtfm::field(arc(0.0, 0.0, 0.5 * length), curved, table)(1), 1.0, 1.0e-6);
}

// toBendArcCoords measures the radial offset away from the centre of curvature
// for either sign, so field() flips it for a negative curvature to get the Frenet
// x that continues the entrance frame.
TEST(MultipoleTFieldModel, FieldFlipsTheRadialOffsetForANegativeCurvature) {
    Table table         = makeTable(5);
    const double length = 1.0;
    const double radial = 0.01;

    const mtfm::FieldInputs positive = makeInputs({0.0, 1.0}, length, +0.5, 0.02);
    const mtfm::FieldInputs negative = makeInputs({0.0, 1.0}, length, -0.5, 0.02);

    // The same radial offset sits at +x when bending one way and at -x the other,
    // so the gradient field comes out with opposite signs.
    const double atPositive = mtfm::field(arc(radial, 0.0, 0.5 * length), positive, table)(1);
    const double atNegative = mtfm::field(arc(radial, 0.0, 0.5 * length), negative, table)(1);
    EXPECT_GT(atPositive, 0.0);
    EXPECT_LT(atNegative, 0.0);
    EXPECT_NEAR(atPositive, -atNegative, 1.0e-6);
}
