//
// Unit tests for class Tanh
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


#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "AbstractObjects/OpalData.h"
#include "Physics/Physics.h"
#include "gtest/gtest.h"

class TestTanh : public testing::Test {
public:
    TestTanh() {
    }

    static void SetUpTestSuite() {
        Kokkos::initialize();
    }
    static void TearDownTestSuite() {
    }
};

TEST_F(TestTanh, ConstructorTest) {
    //x0, lambda, max_index
    auto tanh1 = std::make_shared<endfieldmodel::Tanh>(1, 2, 3);
    EXPECT_EQ(tanh1->getX0(), 1.0);
    EXPECT_EQ(tanh1->getCentreLength(), 2.0);
    EXPECT_EQ(tanh1->getEndLength(), 2.0);

    endfieldmodel::Tanh* tanh2 = tanh1->clone();
    EXPECT_EQ(tanh2->getX0(), 1.0);
    EXPECT_EQ(tanh2->getEndLength(), 2.0);
    delete tanh2;
}

TEST_F(TestTanh, RescaleTest) {
    endfieldmodel::Tanh tanh1(5, 1, 6);
    tanh1.rescale(2);
    EXPECT_EQ(tanh1.getX0(), 10.0);
    EXPECT_EQ(tanh1.getEndLength(), 2.0);
}

TEST_F(TestTanh, FunctionTest) {
    endfieldmodel::Tanh tanh1(5, 1, 6);
    double fCalculated = tanh1.function(5, 0);
    EXPECT_NEAR(fCalculated, 0.5, 1e-6);
    double delta = 1e-9;
    for (auto i: std::vector<int>({1, 2, 3, 4, 5, 6})) {
        double fCalculated = tanh1.function(5, i);
        double fEstimated = (tanh1.function(5+delta, i-1)-tanh1.function(5-delta, i-1))/2/delta;
        EXPECT_NEAR(fCalculated, fEstimated, 1e-6);
    }
}
    // error in GPU land - trying to access HOST memory from DEVICE (or vice
    // versa). Sort it out in the morning.
TEST_F(TestTanh, FunctionGpuTest) {
    endfieldmodel::Tanh tanh1(5, 1, 6);
    Kokkos::View<double*> xgpu("xgpu", 2);
    auto xhost = Kokkos::create_mirror_view(xgpu);
    xhost(0) = -5.0;
    xhost(1) = 5.0;
    Kokkos::deep_copy(xgpu, xhost); // target, source
    Kokkos::View<double**> derivativesgpu("d", 2, 4);

    tanh1.function(xgpu, 4, derivativesgpu);

    auto derivativeshost = Kokkos::create_mirror_view(derivativesgpu);
    Kokkos::deep_copy(derivativeshost, derivativesgpu); // target, source
    for (auto i: std::vector<int>({0, 1, 2, 3})) {
        EXPECT_NEAR(derivativeshost(0, i), tanh1.function(-5, i), 1e-12);
        EXPECT_NEAR(derivativeshost(1, i), tanh1.function(5, i), 1e-12);
    }
}
