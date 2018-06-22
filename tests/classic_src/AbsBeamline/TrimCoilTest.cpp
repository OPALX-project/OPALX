#include "gtest/gtest.h"

#include "TrimCoils/TrimCoilBFit.h"
#include "TrimCoils/TrimCoilPhaseFit.h"
#include "TrimCoils/TrimCoilMirrored.h"

#include "opal_test_utilities/SilenceTest.h"

// Test TrimCoilBFit on zeros
TEST(TrimCoil, TrimCoilBFitZeros)
{
    OpalTestUtilities::SilenceTest silencer;

    const double mm2m = 0.001;

    double bmax = 0.0;
    double rmin = 1000; // mm
    double rmax = 2000;

    TrimCoilBFit myTrimCoil(bmax, rmin, rmax, {}, {});

    const double one = 1.0;
    double br = one, bz = one;
    myTrimCoil.applyField((rmin+rmax)*mm2m/2.0, 1, &br, &bz);
    // not changed since bmax 0.0
    EXPECT_NEAR(br, one, 1e-7);
    EXPECT_NEAR(bz, one, 1e-7);

    bmax = 1.0;
    myTrimCoil = TrimCoilBFit(bmax, rmin, rmax, {}, {});

    myTrimCoil.applyField(rmin*mm2m - 1, 1, &br, &bz);
    // not changed since r outside range
    EXPECT_NEAR(br, one, 1e-7);
    EXPECT_NEAR(bz, one, 1e-7);

    myTrimCoil.applyField(rmax*mm2m + 1, 1, &br, &bz);
    // not changed since r outside range
    EXPECT_NEAR(br, one, 1e-7);
    EXPECT_NEAR(bz, one, 1e-7);
}

TEST(TrimCoil, TrimCoilBFit)
{
    OpalTestUtilities::SilenceTest silencer;

    double bmax = 1.0; // 1 T will be converted to 10 kG
    double rmin = 0.0;
    double rmax = 3000.0;

    // polynom 1 + 2*x + 3*x^2
    TrimCoilBFit myTrimCoil(bmax, rmin, rmax, {1.0, 2.0, 3.0}, {});

    double br = 1.0, bz = 1.0;
    myTrimCoil.applyField(2.0, 2.0, &br, &bz);
    // 1 + 10*(1 + 2*2 + 3*2*2) = 171 ; 1 + 10*2*(2 + 6*2) = 281
    EXPECT_NEAR(bz, 171.0, 1e-7);
    EXPECT_NEAR(br, 281.0, 1e-7);

    // rational function (4 + 3*x) / (1 + 2*x)
    myTrimCoil = TrimCoilBFit(bmax, rmin, rmax, {4.0, 3.0}, {1.0, 2.0});

    br = 1.0, bz = 1.0;
    myTrimCoil.applyField(2.0, 2.0, &br, &bz);
    // 1 + 10*(4+3*2) / (1+2*2) = 21.0 ; 1 + 10*2*(-0.2) = -3
    EXPECT_NEAR(bz, 21.0, 1e-7);
    EXPECT_NEAR(br, -3.0, 1e-7);
}

TEST(TrimCoil, TrimCoilPhaseFit)
{
    OpalTestUtilities::SilenceTest silencer;

    double bmax = 1.0; // 1 T will be converted to 10 kG
    double rmin = 0.0;
    double rmax = 3000.0;

    // polynom 1 + 2*x + 3*x^2
    TrimCoilPhaseFit myTrimCoil(bmax, rmin, rmax, {1.0, 2.0, 3.0}, {});

    double br = 1.0, bz = 1.0;
    myTrimCoil.applyField(2.0, 2.0, &br, &bz);
    // 1 - 10*(2 + 6*2) = -139; 1 - 2*10*6 = -119
    EXPECT_NEAR(bz, -139.0, 1e-7);
    EXPECT_NEAR(br, -119.0, 1e-7);

    // rational function (4 + 3*x) / (1 + 2*x)
    myTrimCoil = TrimCoilPhaseFit(bmax, rmin, rmax, {4.0, 3.0}, {1.0, 2.0});

    br = 1.0, bz = 1.0;
    myTrimCoil.applyField(2.0, 2.0, &br, &bz);
    // 1 - 10*(-0.2) = 3; 1 + 2*10*2*(-0.2*2) / 5 = -2.2   
    EXPECT_NEAR(bz, 3.0, 1e-7);
    EXPECT_NEAR(br,-2.2, 1e-7);
}

TEST(TrimCoil, TrimCoilMirrored)
{
    OpalTestUtilities::SilenceTest silencer;

    double bmax = 1.0; // 1 T will be converted to 10 kG
    double rmin = 0.0;
    double rmax = 3000.0;
    double bslope = 1./6.;

    TrimCoilMirrored myTrimCoil(bmax, rmin, rmax, bslope);

    double br = 1.0, bz = 1.0;
    myTrimCoil.applyField(2.0, 2.0, &br, &bz);

    EXPECT_NEAR(bz,-6.1943868603626751, 1e-6);
    EXPECT_NEAR(br, 1.0032755233321968, 1e-6);
}