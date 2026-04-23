/*
 *  Copyright (c) 2025, Jon Thompson
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */


#include <gtest/gtest.h>
#include <vector>
#include "Utilities/LinearSpline.h"

class LinearSplineTest : public testing::Test {
protected:
    void SetUp() override {
        // Test data: y = x^2
        x_data = {0.0, 1.0, 2.0, 3.0, 4.0};
        y_data = {0.0, 1.0, 4.0, 9.0, 16.0};
    }

    std::vector<double> x_data;
    std::vector<double> y_data;
};

TEST_F(LinearSplineTest, BasicInterpolation) {
    const LinearSpline spline(x_data, y_data);
    LinearSpline::Accelerator acc;
    // Test interpolation at data points
    EXPECT_NEAR(spline.eval(0.0, acc), 0.0, 1e-10);
    EXPECT_NEAR(spline.eval(1.0, acc), 1.0, 1e-10);
    EXPECT_NEAR(spline.eval(2.0, acc), 4.0, 1e-10);
    EXPECT_NEAR(spline.eval(3.0, acc), 9.0, 1e-10);
    EXPECT_NEAR(spline.eval(4.0, acc), 16.0, 1e-10);

    // Test interpolation between points
    EXPECT_NEAR(spline.eval(1.5, acc), 2.5, 1e-10);
}

TEST_F(LinearSplineTest, Uninitialised) {
    // Uninitialised
    LinearSpline::Accelerator accel;
    LinearSpline spline_uninit;
    EXPECT_THROW(spline_uninit.eval(0.0, accel), std::runtime_error);
}

TEST_F(LinearSplineTest, Extrapolation) {
    const LinearSpline spline(x_data, y_data);
    LinearSpline::Accelerator acc;

    // Test left extrapolation
    EXPECT_NEAR(spline.eval(-1.0, acc), -1.0, 1e-10);
    EXPECT_NEAR(spline.eval(-0.5, acc), -0.5, 1e-10);

    // Test right extrapolation
    EXPECT_NEAR(spline.eval(5.0, acc), 23.0, 1e-10);
    EXPECT_NEAR(spline.eval(5.5, acc), 26.5, 1e-10);
}

TEST_F(LinearSplineTest, Accelerator) {
    const LinearSpline spline(x_data, y_data);
    LinearSpline::Accelerator acc;

    // First evaluation
    EXPECT_NEAR(spline.eval(1.5, acc), 2.5, 1e-10);

    // Second evaluation in the same interval (should use cached index)
    EXPECT_NEAR(spline.eval(1.5, acc), 2.5, 1e-10);

    // Evaluation in a different interval
    EXPECT_NEAR(spline.eval(2.5, acc), 6.5, 1e-10);
}

TEST_F(LinearSplineTest, InvalidInput) {
    // Too few points
    EXPECT_THROW(LinearSpline spline({0.0}, {0.0}), std::invalid_argument);
    // Differing lengths
    EXPECT_THROW(LinearSpline spline({0.0}, {0.0, 1.0}), std::invalid_argument);
    // Non-increasing x
    EXPECT_THROW(LinearSpline spline({0.0, 1.0, 0.5, 2.0}, {0.0, 1.0, 2.0, 3.0}),
        std::invalid_argument);
}

TEST_F(LinearSplineTest, Integration) {
    LinearSpline spline(x_data, y_data);
    LinearSpline::Accelerator accel;

    // Whole pre-computed intervals
    EXPECT_NEAR( spline.evalIntegral(0, 1, accel), 0.5, 1e-10);
    EXPECT_NEAR( spline.evalIntegral(1, 2, accel), 2.5, 1e-10);
    EXPECT_NEAR( spline.evalIntegral(2, 3, accel), 6.5, 1e-10);
    EXPECT_NEAR( spline.evalIntegral(3, 4, accel), 12.5, 1e-10);
    EXPECT_NEAR( spline.evalIntegral(0, 4, accel), 22.0, 1e-10);

    // Partial first and last intervals
    EXPECT_NEAR( spline.evalIntegral(0.6, 3.2, accel), 11.26, 1e-10);

    // Asking again should use the cache
    EXPECT_NEAR( spline.evalIntegral(0.6, 3.2, accel), 11.26, 1e-10);

    // Extrapolation below the first interval
    EXPECT_NEAR( spline.evalIntegral(-1, 1.5, accel), 0.875, 1e-10);

    // Extrapolation above the last interval
    EXPECT_NEAR( spline.evalIntegral(3, 4.5, accel), 21.375, 1e-10);

    // Swapping integral bounds
    EXPECT_NEAR( spline.evalIntegral(2, 1, accel), 2.5, 1e-10);

    // Uninitialised spline
    LinearSpline spline_uninit;
    EXPECT_THROW(spline_uninit.evalIntegral(2, 1, accel), std::runtime_error);

    // Re-initialise the spline to something different
    spline.init({0.0, 1.0, 2.0, 3.0, 4.0}, {0.0, 2.0, 4.0, 6.0, 8.0});
    EXPECT_NEAR( spline.evalIntegral(0, 4, accel), 16.0, 1e-10);
}

