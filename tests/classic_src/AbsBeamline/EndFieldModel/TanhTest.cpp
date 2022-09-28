#include <cmath>
#include "gtest/gtest.h"
#include "Classic/AbsBeamline/EndFieldModel/Tanh.h"

TEST(TanhTest, FunctionTest) {
    endfieldmodel::Tanh tanh(10, 0.5, 20);
    EXPECT_NEAR(tanh.function(0.0, 0), 1.0, 1e-6);
    EXPECT_NEAR(tanh.function(-10.0, 0), 0.5, 1e-6);
    EXPECT_NEAR(tanh.function(10.0, 0), 0.5, 1e-6);

    double deriv = 0.0;
    double dx = 1e-6;
    for(size_t i = 0; i < 10; ++i) {
        deriv = 0.5*(tanh.function(10.0+dx, i)-tanh.function(10.0-dx, i))/dx;
        double delta = std::fabs(deriv*1e-3);
        EXPECT_NEAR(tanh.function(10.0, i+1), deriv, std::max(delta, 1e-3)) << i;
    }
}

TEST(TanhTest, RescaleTest) {
    endfieldmodel::Tanh tanh(7.323, 0.32, 20);
    endfieldmodel::Tanh* tanh2 = tanh.clone(); 
    ASSERT_NE(tanh2, nullptr);
    tanh2->rescale(0.1);
    double derivativeScale = 1;
    for (size_t i = 0; i < 3; ++i) {
        // if f(x) -> f(10x) then df/dx -> 10 df/dx
        EXPECT_NEAR(tanh.function(9.32, i)*derivativeScale, tanh2->function(0.932, i), 1e-9) << i;
        derivativeScale *= 10;
    }
}
