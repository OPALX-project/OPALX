#include <gtest/gtest.h>

#include "PartBunch/Corrections/ImageChargeCorrection.hpp"
#include "PartBunch/Corrections/ShiftedGreensCorrection.hpp"

TEST(BoundaryCorrectionTest, StepBudgetIsOwnedByCorrection) {
    ImageChargeCorrection<double, 3> correction;
    correction.configure(true, 0.0);

    EXPECT_TRUE(correction.isActiveForStep(0, 0));
    EXPECT_TRUE(correction.isActiveForStep(4, 5));
    EXPECT_FALSE(correction.isActiveForStep(5, 5));

    correction.configure(false, 0.0);
    EXPECT_FALSE(correction.isActiveForStep(0, 0));
}

TEST(BoundaryCorrectionTest, ImageChargePassUsesImageScatterAndNegativeBFieldSign) {
    ImageChargeCorrection<double, 3> correction;
    correction.configure(true, 0.0);

    BinContext<double, 3> context;
    Vector_t<double, 3> origin(0.0);

    const auto pass = correction.makePass(context, origin, 8);
    ASSERT_TRUE(pass.has_value());
    EXPECT_EQ(pass->scatterMode, ImageScatterMode::ImageOnly);
    EXPECT_FALSE(pass->solveRequest.hasShiftedGreens());
    EXPECT_DOUBLE_EQ(pass->bFieldSign, -1.0);
    EXPECT_EQ(pass->flipAxis, -1);
    EXPECT_TRUE(pass->dumpDirichletPlane);
    EXPECT_TRUE(pass->forceSkipFieldDump);
}

TEST(BoundaryCorrectionTest, ShiftedGreensRejectsBackendWithoutCapability) {
    ShiftedGreensCorrection<double, 3> correction;
    SolverCapabilities capabilities;

    EXPECT_THROW(correction.configure(true, 0.0, capabilities), OpalException);
}

TEST(BoundaryCorrectionTest, ShiftedGreensPassBuildsSolveRequestAndZFlip) {
    ShiftedGreensCorrection<double, 3> correction;
    SolverCapabilities capabilities;
    capabilities.supportsShiftedGreens = true;
    correction.configure(true, 0.5, capabilities);

    BinContext<double, 3> context;
    context.solveSpacing = Vector_t<double, 3>(1.0);
    context.solveSpacing[2] = 0.25;

    Vector_t<double, 3> origin(0.0);
    origin[2] = 1.0;

    const auto pass = correction.makePass(context, origin, 8);
    ASSERT_TRUE(pass.has_value());
    ASSERT_TRUE(pass->solveRequest.hasShiftedGreens());

    EXPECT_EQ(pass->scatterMode, ImageScatterMode::PrimaryOnly);
    EXPECT_DOUBLE_EQ((*pass->solveRequest.greensShift)[0], 0.0);
    EXPECT_DOUBLE_EQ((*pass->solveRequest.greensShift)[1], 0.0);
    EXPECT_DOUBLE_EQ((*pass->solveRequest.greensShift)[2], 3.0);
    EXPECT_DOUBLE_EQ(pass->bFieldSign, -1.0);
    EXPECT_EQ(pass->flipAxis, 2);
    EXPECT_FALSE(pass->dumpDirichletPlane);
    EXPECT_FALSE(pass->forceSkipFieldDump);
}
