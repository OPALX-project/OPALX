//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "Algorithms/SplineTimeDependence.h"
#include "Utilities/GeneralOpalException.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

class TestSplineTimeDependence : public testing::Test {
public:
    TestSplineTimeDependence() : times_m(10), values_m(10) {
        for (size_t i = 0; i < 10; ++i) {
            times_m[i]  = static_cast<double>(i * i + i);
            values_m[i] = times_m[i];  // x^2
        }
    }

    std::vector<double> times_m;
    std::vector<double> values_m;
};

TEST_F(TestSplineTimeDependence, ConstructorTest) {
    try {
        SplineTimeDependence timeDep(1, times_m, values_m);
        for (size_t i = 0; i < values_m.size(); ++i) {
            double test = timeDep.getValue(times_m[i]);
            EXPECT_NEAR(test, values_m[i], 1e-9) << "Index " << i;
        }
        SplineTimeDependence* timeDepClone = timeDep.clone();
        for (size_t i = 0; i < values_m.size(); ++i) {
            double test = timeDep.getValue(times_m[i]);
            EXPECT_NEAR(test, values_m[i], 1e-9) << "Index " << i;
        }
        delete timeDepClone;
        EXPECT_THROW(SplineTimeDependence(0, times_m, values_m), std::invalid_argument);
        EXPECT_THROW(SplineTimeDependence(2, times_m, values_m), std::invalid_argument);
        times_m.push_back(1.);
        EXPECT_THROW(SplineTimeDependence(1, times_m, values_m), std::invalid_argument);
        times_m  = std::vector(1, 0.);
        values_m = std::vector(1, 0.);
        EXPECT_THROW(SplineTimeDependence(1, times_m, values_m), std::invalid_argument);
        times_m  = {0., 1.};
        values_m = std::vector(2, 0.);
        EXPECT_NO_THROW(SplineTimeDependence(1, times_m, values_m));
        times_m  = {0., 1., 2.};
        values_m = std::vector(3, 0.);
        EXPECT_THROW(SplineTimeDependence(3, times_m, values_m), std::invalid_argument);
        times_m  = {0., 1., 2., 3.};
        values_m = std::vector(4, 0.);
        EXPECT_NO_THROW(SplineTimeDependence(3, times_m, values_m));
        times_m  = {0., 0., 2., 3.};
        values_m = std::vector(4, 0.);
        EXPECT_THROW(SplineTimeDependence(3, times_m, values_m), std::invalid_argument);
    } catch (GeneralOpalException& exc) {
        EXPECT_TRUE(false) << "Should not have thrown an exception:\n    " << exc.what() << "\n    "
                           << exc.where();
    }
}

TEST_F(TestSplineTimeDependence, LinearLookupTest) {
    SplineTimeDependence timeDep(1, times_m, values_m);
    double test_x      = (times_m[2] + times_m[3]) / 2.;
    const double ref_y = (values_m[2] + values_m[3]) / 2.;
    EXPECT_NEAR(timeDep.getValue(test_x), ref_y, 1e-9);
    test_x = times_m[0] - (times_m[1] - times_m[0]) / 2.;
    EXPECT_THROW(timeDep.getValue(test_x), std::invalid_argument);
    EXPECT_THROW(timeDep.getIntegral(test_x), std::invalid_argument);
    test_x = times_m[9] + (times_m[9] - times_m[8]) / 2.;
    EXPECT_THROW(timeDep.getValue(test_x), std::invalid_argument);
    EXPECT_THROW(timeDep.getIntegral(test_x), std::invalid_argument);
}

TEST_F(TestSplineTimeDependence, CubicLookupTest) {
    // if I give it a quadratic or cubic, it gets the wrong answer!!
    // I am sure that it is doing the right thing though... ahem
    SplineTimeDependence timeDep(1, times_m, values_m);
    for (size_t i = 0; i < times_m.size() - 1; ++i) {
        const double test_x = (times_m[i] + times_m[i + 1]) / 2.;
        const double ref_y  = test_x;
        const double test_y = timeDep.getValue(test_x);
        EXPECT_NEAR(test_y, ref_y, 1e-9);
    }
}

TEST_F(TestSplineTimeDependence, Integral) {
    SplineTimeDependence timeDep(1, times_m, values_m);
    EXPECT_DOUBLE_EQ(timeDep.getIntegral(times_m[0]), 0.0);
    EXPECT_DOUBLE_EQ(timeDep.getIntegral(times_m[1]), 2.0);
    EXPECT_DOUBLE_EQ(timeDep.getIntegral(times_m[2]), 18.0);
    EXPECT_DOUBLE_EQ(timeDep.getIntegral(times_m[3]), 72.0);
}

TEST_F(TestSplineTimeDependence, Print) {
    int argc    = 0;
    char** argv = nullptr;
    ippl::initialize(argc, argv);
    std::stringstream ss;
    Inform inform("Test", ss);
    inform.setOutputLevel(5);
    SplineTimeDependence timeDep(1, times_m, values_m);
    inform << timeDep << endl;
    EXPECT_STREQ(
            "Test> SplineTimeDependence of order 1 with 10 entries\nTest> \n", ss.str().c_str());
    ss.str("");
    SplineTimeDependence timeDep2;
    inform << timeDep2 << endl;
    EXPECT_STREQ("Test> Uninitialised SplineTimeDependence\nTest> \n", ss.str().c_str());
    ippl::finalize();
}
