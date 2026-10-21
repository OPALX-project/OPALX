//
// Unit tests for class Enge
//
// Copyright (c) 2017-2026, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
// All rights reserved.
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

#include "AbsBeamline/EndFieldModel/Enge.h"
#include "AbstractObjects/OpalData.h"
#include "Physics/Physics.h"
#include "gtest/gtest.h"

class TestEnge : public testing::Test {
public:
    TestEnge() {}

    static void SetUpTestSuite() { Kokkos::initialize(); }
    static void TearDownTestSuite() { Kokkos::finalize(); }
};

TEST_F(TestEnge, ConstructorTest) {
    // x0, lambda, max_index
    endfieldmodel::Enge enge0;
    EXPECT_EQ(enge0.getX0(), 0.0);
    EXPECT_EQ(enge0.getLambda(), 0.0);
    EXPECT_EQ(enge0.getCoefficients().size(), 0);

    std::vector<double> a = {5, 6};
    auto enge1            = std::make_shared<endfieldmodel::Enge>(a, 2, 7);
    enge1->initialise();
    EXPECT_EQ(enge1->getX0(), 2.0);
    EXPECT_EQ(enge1->getCentreLength(), 4.0);
    EXPECT_EQ(enge1->getEndLength(), 7.0);
    EXPECT_EQ(enge1->getCoefficients(), a);

    endfieldmodel::Enge* enge2 = enge1->clone();
    enge2->initialise();
    EXPECT_EQ(enge2->getX0(), enge1->getX0());
    EXPECT_EQ(enge2->getCentreLength(), enge1->getCentreLength());
    EXPECT_EQ(enge2->getCoefficients(), enge2->getCoefficients());
    delete enge2;
}

double myEnge(double x, std::vector<double> a, double lambda) {
    double deltaX = x / lambda;
    double xPow   = 1.0;
    double p      = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        p += a[i] * xPow;
        xPow *= deltaX;
    }
    double enge = 1 / (1.0 + exp(p));
    return enge;
}

TEST_F(TestEnge, DerivativeTest) {
    std::vector<double> xVector = {-5.5, -5.0, -3.0, 0.0, 3.0, 5.0, 5.5};
    double dx                   = 1e-4;
    double x0                   = 5.0;
    double lambda               = 0.5;
    std::vector<double> engeA   = {0.0, 1.0, 0.0};
    endfieldmodel::Enge enge(engeA, x0, lambda);
    enge.initialise();
    Kokkos::View<double*> xDevice("x", 1);
    auto xHost = Kokkos::create_mirror_view(xDevice);
    Kokkos::View<double**> yDevice("y", 1, 5);
    for (auto x : xVector) {
        double yTest = enge.function(x, 0);
        double yRef  = -1 + myEnge(x - x0, engeA, lambda) + myEnge(-x - x0, engeA, lambda);
        xHost(0)     = x;
        Kokkos::deep_copy(xDevice, xHost);
        enge.function(xDevice, 4, yDevice);
        auto yHost = Kokkos::create_mirror_view(yDevice);
        Kokkos::deep_copy(yHost, yDevice);
        EXPECT_NEAR(yTest, yRef, 1e-12) << x << " " << yTest << " " << yRef;
        EXPECT_NEAR(yHost(0, 0), yRef, 1e-12) << x << " " << yTest << " " << yRef;
        for (size_t n = 1; n < 5; ++n) {
            double yP     = enge.function(x + dx, n - 1);
            double yM     = enge.function(x - dx, n - 1);
            double dyTest = enge.function(x, n);
            double dyRef  = (yP - yM) / 2 / dx;
            EXPECT_NEAR(dyTest, dyRef, 1e-3) << x << " " << n << " * " << dyTest << " " << dyRef;
            EXPECT_NEAR(yHost(0, n), dyRef, 1e-3)
                    << x << " " << n << " * " << yHost(0, n) << " " << dyRef;
            ;
        }
    }
    EXPECT_NEAR(enge.function(0.0, 0), 1.0, 1e-3);
    EXPECT_NEAR(enge.function(5.0, 0), 0.5, 1e-3);
    EXPECT_NEAR(enge.function(-5.0, 0), 0.5, 1e-3);
}

TEST_F(TestEnge, PrintTest) {
    std::vector<double> a    = {1, 2};
    endfieldmodel::Enge enge = endfieldmodel::Enge(a, 1.0, 0.5);
    enge.initialise();
    enge.print(std::cout);
    std::cout << std::endl;
}

TEST_F(TestEnge, FunctionGpuTest) {
    endfieldmodel::Enge enge(std::vector<double>({0, 1, 0}), 5, 1);
    enge.initialise();
    Kokkos::View<double*> xgpu("xgpu", 3);
    auto xhost = Kokkos::create_mirror_view(xgpu);
    xhost(0)   = 5.0;
    xhost(1)   = -5.0;
    xhost(2)   = 7.0;
    Kokkos::deep_copy(xgpu, xhost);  // target, source
    Kokkos::View<double**> derivativesgpu("d", 3, 5);

    enge.function(xgpu, 4, derivativesgpu);

    auto derivativeshost = Kokkos::create_mirror_view(derivativesgpu);
    Kokkos::deep_copy(derivativeshost, derivativesgpu);  // target, source
    for (auto i : std::vector<int>({0, 1, 2, 3, 4})) {
        EXPECT_NEAR(derivativeshost(0, i), enge.function(5, i), 1e-12);
        EXPECT_NEAR(derivativeshost(1, i), enge.function(-5, i), 1e-12);
        EXPECT_NEAR(derivativeshost(2, i), enge.function(7, i), 1e-12);
        std::cout << i << " " << xhost(0) << " " << xhost(1) << " " << xhost(2)
                  << " gpu-5: " << derivativeshost(0, i) << " cpu-5: " << enge.function(-5, i)
                  << " gpu+5 " << derivativeshost(1, i) << " cpu+5: " << enge.function(5, i)
                  << " gpu+7 " << derivativeshost(2, i) << " cpu+7: " << enge.function(7, i)
                  << std::endl;
    }
}
