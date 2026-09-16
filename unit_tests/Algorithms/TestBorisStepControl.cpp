#include "Algorithms/BorisStepControl.h"

#include <limits>

#include "gtest/gtest.h"

TEST(BorisStepControl, CoordinateFloorUsesActualSpeedAndPreservesResolvableDistance) {
    const double c = 299792458.;
    const double beta = 0.07285411307707136;
    for (double scale : {0., 0.001, 1., 8.4, 1e6}) {
        const double fast = boris_step::positionTimeFloor(scale, c);
        const double slow = boris_step::positionTimeFloor(scale, beta * c);
        const double distance = 64 * std::numeric_limits<double>::epsilon() * std::max(1., scale);
        EXPECT_NEAR(slow / fast, 1 / beta, 4e-15);
        EXPECT_NEAR(fast * c, distance, 2 * std::numeric_limits<double>::epsilon() * distance);
        EXPECT_NEAR(slow * beta * c, distance, 2 * std::numeric_limits<double>::epsilon() * distance);
        EXPECT_GT(slow, fast);
    }
}

TEST(BorisStepControl, CoordinateFloorRejectsInvalidScaleAndSpeed) {
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (double speed : {0., -1., infinity, nan})
        EXPECT_THROW(boris_step::positionTimeFloor(1, speed), std::invalid_argument);
    for (double scale : {-1., infinity, nan})
        EXPECT_THROW(boris_step::positionTimeFloor(scale, 1), std::invalid_argument);
    EXPECT_THROW(boris_step::positionTimeFloor(std::numeric_limits<double>::max(),
            std::numeric_limits<double>::min()), std::overflow_error);
}

TEST(BorisStepControl, AccumulatedSolverTranslationKeepsHalfDriftsRepresentable) {
    const double nominalDt = 6.25e-11;
    const double speed = 0.07285411307707136 * 299792458.;
    const double laboratoryScale = 2.;
    const double solverTranslation = 825.; // Accumulated path after many ring turns [m].
    boris_step::Control solverAware(nominalDt, boris_step::positionTimeFloor(
            std::max(laboratoryScale, solverTranslation), speed));
    while (solverAware.canSplit()) solverAware.split();

    // Even the half drift remains well above the solver-coordinate roundoff.
    // Control may bisect once below its floor, so its smallest half drift is
    // bounded by 16 (rather than 64) coordinate roundoff scales.
    const double halfDrift = 0.5 * speed * solverAware.step();
    EXPECT_GE(halfDrift, 16 * std::numeric_limits<double>::epsilon() * solverTranslation);
    EXPECT_GT((solverTranslation + halfDrift) - solverTranslation, 0.);

    // Using only the bounded laboratory ring coordinate reproduces the lost
    // motion that motivated including the PIC frame's path translation.
    boris_step::Control laboratoryOnly(nominalDt,
            boris_step::positionTimeFloor(laboratoryScale, speed));
    while (laboratoryOnly.canSplit()) laboratoryOnly.split();
    const double unresolvedHalfDrift = 0.5 * speed * laboratoryOnly.step();
    EXPECT_DOUBLE_EQ(solverTranslation + unresolvedHalfDrift, solverTranslation);
    EXPECT_GT(solverAware.step(), laboratoryOnly.step());
}

TEST(BorisStepControl, AcceptsOneWholeStepAndRejectsUseAfterCompletion) {
    boris_step::Control control(6.25e-11);
    EXPECT_FALSE(control.done());
    EXPECT_DOUBLE_EQ(control.elapsed(), 0);
    EXPECT_DOUBLE_EQ(control.step(), 6.25e-11);
    control.accept();
    EXPECT_TRUE(control.done());
    EXPECT_FALSE(control.canSplit());
    EXPECT_DOUBLE_EQ(control.elapsed(), 6.25e-11);
    EXPECT_THROW(control.step(), std::logic_error);
    EXPECT_THROW(control.split(), std::logic_error);
    EXPECT_THROW(control.accept(), std::logic_error);
}

TEST(BorisStepControl, RetainsBothHalvesInChronologicalOrder) {
    boris_step::Control control(8e-10);
    control.split();
    control.split();
    control.accept(); // First quarter.
    EXPECT_DOUBLE_EQ(control.elapsed(), 2e-10);
    EXPECT_DOUBLE_EQ(control.step(), 2e-10);
    control.accept(); // Second quarter.
    EXPECT_DOUBLE_EQ(control.elapsed(), 4e-10);
    EXPECT_DOUBLE_EQ(control.step(), 4e-10);
    control.split();
    control.accept(); // Third quarter.
    EXPECT_DOUBLE_EQ(control.elapsed(), 6e-10);
    control.split();
    EXPECT_DOUBLE_EQ(control.step(), 1e-10);
    control.accept();
    control.accept();
    EXPECT_TRUE(control.done());
    EXPECT_DOUBLE_EQ(control.elapsed(), 8e-10);
}

TEST(BorisStepControl, UniformSubdivisionNeverDropsAcceptedTime) {
    const double nominal = 7.123456789e-11;
    const double smallest = std::ldexp(nominal, -16);
    boris_step::Control control(nominal);
    std::size_t accepted = 0;
    while (!control.done()) {
        if (control.step() > smallest) {
            control.split();
        } else {
            EXPECT_DOUBLE_EQ(control.step(), smallest);
            control.accept();
            ++accepted;
            EXPECT_GT(control.elapsed(), 0);
            EXPECT_LE(control.elapsed(), nominal);
        }
    }
    EXPECT_EQ(accepted, 65536u);
    // Multiplication by this power-of-two leaf count is exact. A plain sum
    // would accumulate rounding error (long double also equals double on some
    // supported hosts), whereas the scheduler retains compensated elapsed time.
    EXPECT_DOUBLE_EQ(static_cast<double>(accepted) * smallest, nominal);
    EXPECT_DOUBLE_EQ(control.elapsed(), nominal);
}

TEST(BorisStepControl, FloorIsRelativeToRootRatherThanShrinkingWithChildren) {
    for (double nominal : {1e-100, 6.25e-11, 1e100}) {
        boris_step::Control control(nominal);
        unsigned depth = 0;
        while (control.canSplit()) {
            control.split();
            ++depth;
        }
        EXPECT_EQ(depth, 40u);
        EXPECT_LE(control.step(), 1e-12 * nominal);
        const double stopped = control.step();
        EXPECT_THROW(control.split(), std::runtime_error);
        EXPECT_DOUBLE_EQ(control.step(), stopped);
        unsigned accepted = 0;
        while (!control.done()) {
            control.accept();
            ++accepted;
        }
        EXPECT_EQ(accepted, 41u);
        EXPECT_DOUBLE_EQ(control.elapsed(), nominal);
    }
}

TEST(BorisStepControl, HonorsCallerFloorAndAvoidsUnderflow) {
    boris_step::Control control(8e-10, 2e-10);
    control.split();
    control.split();
    EXPECT_FALSE(control.canSplit());
    EXPECT_THROW(control.split(), std::runtime_error);
    boris_step::Control largerFloor(1e-10, 2e-10);
    EXPECT_FALSE(largerFloor.canSplit());
    largerFloor.accept();
    EXPECT_DOUBLE_EQ(largerFloor.elapsed(), 1e-10);

    const double smallest = std::numeric_limits<double>::denorm_min();
    boris_step::Control subnormal(smallest);
    EXPECT_FALSE(subnormal.canSplit());
    EXPECT_THROW(subnormal.split(), std::runtime_error);
    subnormal.accept();
    EXPECT_DOUBLE_EQ(subnormal.elapsed(), smallest);

    boris_step::Control oddSubnormal(3 * smallest);
    EXPECT_FALSE(oddSubnormal.canSplit());
    EXPECT_THROW(oddSubnormal.split(), std::runtime_error);
    oddSubnormal.accept();
    EXPECT_DOUBLE_EQ(oddSubnormal.elapsed(), 3 * smallest);

    boris_step::Control evenSubnormal(4 * smallest);
    evenSubnormal.split();
    evenSubnormal.split();
    EXPECT_FALSE(evenSubnormal.canSplit());
    EXPECT_DOUBLE_EQ(evenSubnormal.step(), smallest);
    evenSubnormal.accept();
    evenSubnormal.accept();
    EXPECT_DOUBLE_EQ(evenSubnormal.step(), 2 * smallest);
    evenSubnormal.accept();
    EXPECT_DOUBLE_EQ(evenSubnormal.elapsed(), 4 * smallest);
}

TEST(BorisStepControl, BoundsTrialsAndPreservesStateOnBudgetFailure) {
    boris_step::Control one(1, 0, 1);
    EXPECT_THROW(one.split(), std::runtime_error);
    EXPECT_DOUBLE_EQ(one.step(), 1);
    EXPECT_DOUBLE_EQ(one.elapsed(), 0);
    one.accept();
    EXPECT_TRUE(one.done());

    boris_step::Control two(1, 0, 2);
    two.split();
    EXPECT_THROW(two.accept(), std::runtime_error);
    EXPECT_THROW(two.split(), std::runtime_error);
    EXPECT_DOUBLE_EQ(two.step(), 0.5);
    EXPECT_DOUBLE_EQ(two.elapsed(), 0);

    boris_step::Control three(1, 0, 3);
    three.split();
    three.accept();
    three.accept();
    EXPECT_TRUE(three.done());
}

TEST(BorisStepControl, RejectsInvalidConstruction) {
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (double dt : {0., -1., infinity, nan})
        EXPECT_THROW(boris_step::Control control(dt), std::invalid_argument);
    for (double floor : {-1., infinity, nan})
        EXPECT_THROW(boris_step::Control control(1, floor), std::invalid_argument);
    EXPECT_THROW(boris_step::Control control(1, 0, 0), std::invalid_argument);
}
