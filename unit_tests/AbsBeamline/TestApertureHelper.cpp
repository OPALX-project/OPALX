/**
 * \file TestApertureHelper.cpp
 * \brief Unit tests for the device-callable transverse-aperture predicate.
 *
 * ---------------------------------------------------------------------------
 * Coverage
 * ---------------------------------------------------------------------------
 *
 * 1. Pure geometry (plain TESTs, no IPPL)
 *    - RECTANGULAR: inside/outside, strict boundary
 *    - ELLIPTICAL: inside/outside, strict boundary, sign symmetry
 *    - the default aperture {1e6, 1e6} accepts any realistic point
 *
 * 2. Equivalence with the element-level check (fixture, IPPL)
 *    - probes a DriftRep through the public isInside() for in-window z and
 *      compares against the scalar helper for both aperture types
 */

#include <gtest/gtest.h>

#include "AbsBeamline/ElementBase.h"

#include "Ippl.h"

#include "BeamlineCore/DriftRep.h"

namespace {
    bool inside(double x, double y, ApertureType type, double xLim, double yLim) {
        return ApertureHelper::isInsideAperture(x, y, type, xLim, yLim);
    }
}  // namespace

// ---------------------------------------------------------------------------
// RECTANGULAR
// ---------------------------------------------------------------------------
TEST(ApertureHelperTest, RectangularInsideOutsideAndBoundary) {
    const double xLim = 0.02, yLim = 0.01;
    auto rect = [&](double x, double y) {
        return inside(x, y, ApertureType::RECTANGULAR, xLim, yLim);
    };

    EXPECT_TRUE(rect(0.0, 0.0));
    EXPECT_TRUE(rect(0.019, -0.009));
    EXPECT_TRUE(rect(-0.019, 0.009));

    // one coordinate out
    EXPECT_FALSE(rect(0.021, 0.0));
    EXPECT_FALSE(rect(0.0, -0.011));

    // boundary is strict: exactly at the limit is outside
    EXPECT_FALSE(rect(xLim, 0.0));
    EXPECT_FALSE(rect(0.0, yLim));
    EXPECT_FALSE(rect(-xLim, 0.0));
}

// ---------------------------------------------------------------------------
// ELLIPTICAL
// ---------------------------------------------------------------------------
TEST(ApertureHelperTest, EllipticalInsideOutsideAndBoundary) {
    const double a = 0.03, b = 0.02;
    auto ell = [&](double x, double y) {
        return inside(x, y, ApertureType::ELLIPTICAL, a, b);
    };

    EXPECT_TRUE(ell(0.0, 0.0));
    EXPECT_TRUE(ell(0.029, 0.0));
    EXPECT_TRUE(ell(0.0, -0.019));
    // on-diagonal point inside: (x/a)^2 + (y/b)^2 = 0.5
    EXPECT_TRUE(ell(a / 2.0, b / 2.0));

    // semi-axis points are on the boundary -> outside (strict <)
    EXPECT_FALSE(ell(a, 0.0));
    EXPECT_FALSE(ell(0.0, b));
    EXPECT_FALSE(ell(-a, 0.0));

    // inside the bounding rectangle but outside the ellipse
    EXPECT_FALSE(ell(0.028, 0.014));
}

// ---------------------------------------------------------------------------
// Default aperture
// ---------------------------------------------------------------------------
TEST(ApertureHelperTest, DefaultHugeApertureAcceptsEverything) {
    // ElementBase default: ELLIPTICAL {1e6, 1e6}
    EXPECT_TRUE(inside(10.0, -10.0, ApertureType::ELLIPTICAL, 1e6, 1e6));
    EXPECT_TRUE(inside(-100.0, 100.0, ApertureType::ELLIPTICAL, 1e6, 1e6));
}

// ---------------------------------------------------------------------------
// Equivalence with the element-level aperture check via DriftRep::isInside
// ---------------------------------------------------------------------------
class ApertureHelperElementTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() { ippl::finalize(); }
};

TEST_F(ApertureHelperElementTest, MatchesElementBaseApertureCheck) {
    const double L                                                            = 1.0;
    const std::vector<std::pair<ApertureType, std::vector<double>>> apertures = {
            {ApertureType::RECTANGULAR, {0.02, 0.01}},
            {ApertureType::ELLIPTICAL, {0.03, 0.02}},
    };

    DriftRep drift("test_drift");
    drift.getGeometry().setElementLength(L);

    // probe grid: transverse points around the limits, several z inside [0, L)
    const std::vector<double> coords = {-0.031, -0.02, -0.015, -0.009, 0.0,
                                        0.009,  0.015, 0.02,   0.031};
    const std::vector<double> zs     = {0.0, 0.25, 0.5, 0.999};

    for (const auto& [type, args] : apertures) {
        drift.setAperture(type, args);
        for (double x : coords) {
            for (double y : coords) {
                for (double z : zs) {
                    const Vector_t<double, 3> r(x, y, z);
                    // isInside = z-window && aperture test; z is chosen in-window,
                    // so any disagreement is an aperture-predicate disagreement.
                    EXPECT_EQ(
                            drift.isInside(r),
                            ApertureHelper::isInsideAperture(x, y, type, args[0], args[1]))
                            << "type=" << static_cast<int>(type) << " x=" << x << " y=" << y
                            << " z=" << z;
                }
            }
        }
    }
}
