//
// Tests for MultipoleTBase
//
// Copyright (c) 2025, Jon Thompson, STFC Rutherford Appleton Laboratory, Didcot, UK
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//
#include "gtest/gtest.h"
#include "AbsBeamline/MultipoleTBase.h"

class MultipoleTBaseTest : public testing::Test, public MultipoleTBase {
public:
    MultipoleTBaseTest() : MultipoleTBase(nullptr) {}

protected:
    static void SetUpTestSuite() {
        int argc = 0;
        char** argv = nullptr;

        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() {
        ippl::finalize();
    }

    void SetUp() override {
        // nothing special
    }

    void TearDown() override {
        // nothing special
    }

    // Overrides of MultipoleTBase
    void initialise() override {}
    BGeometryBase* getGeometry() override { return nullptr; }
    const BGeometryBase* getGeometry() const override { return nullptr; }
    void transformCoords(Vector_t<double>& /*R*/) override {}
    void transformBField(Vector_t<double>& /*B*/, const Vector_t<double>& /*R*/) override {}
    double getScaleFactor(double /*x*/, double /*s*/) override { return 1.0; }
    double getFn(unsigned int /*n*/, double /*x*/, double /*s*/) override { return 0.0; }
    Vector_t<double> localCartesianToOpalCartesian(const Vector_t<double>& r) override { return r; }

    // Tanh derivative coefficient test helper
    std::vector<double> tanhCoefficients(const unsigned int derivative) const {
        std::vector<double> coefficients;
        const auto numCoefficients = tanhCoefficients_m.extent(1);
        coefficients.resize(numCoefficients);
        for(unsigned int i = 0; i < numCoefficients ; ++i) {
            coefficients[i] = tanhCoefficients_m(derivative, i);
        }
        return coefficients;
    }
};

TEST_F(MultipoleTBaseTest, TransverseDerivatives) {
    constexpr std::array poles = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    std::array<double, MaxDerivatives> derivatives;
    calcTransverseDerivatives(poles, 6, 1.0, derivatives);
    EXPECT_DOUBLE_EQ(derivatives[0], 55.0);
    EXPECT_DOUBLE_EQ(derivatives[1], 170.0);
    EXPECT_DOUBLE_EQ(derivatives[2], 414.0);
    EXPECT_DOUBLE_EQ(derivatives[3], 696.0);
    EXPECT_DOUBLE_EQ(derivatives[4], 600.0);
    EXPECT_DOUBLE_EQ(derivatives[5], 0.0);
}

TEST_F(MultipoleTBaseTest, FringeDerivativeConstants) {
    generateTanhCoefficients(5);
    EXPECT_EQ(tanhCoefficients(0), (std::vector<double>{0, 1, 0, 0, 0, 0, 0}));
    EXPECT_EQ(tanhCoefficients(1), (std::vector<double>{1, 0, -1, 0, 0, 0, 0}));
    EXPECT_EQ(tanhCoefficients(2), (std::vector<double>{0, -2, 0, 2, 0, 0, 0}));
    EXPECT_EQ(tanhCoefficients(3), (std::vector<double>{-2, 0, 8, 0, -6, 0, 0}));
    EXPECT_EQ(tanhCoefficients(4), (std::vector<double>{0, 16, 0, -40, 0, 24, 0}));
    EXPECT_EQ(tanhCoefficients(5), (std::vector<double>{16, 0, -136, 0, 240, 0, -120}));
}

TEST_F(MultipoleTBaseTest, FringeDerivatives) {
    std::array<double, MaxDerivatives> derivativesNeg;
    std::array<double, MaxDerivatives> derivativesPos;
    generateTanhCoefficients(6);
    calcFringeDerivatives(2.0, 1.0, 1.0, -2.0, tanhCoefficients_m, derivativesNeg);
    calcFringeDerivatives(2.0, 1.0, 1.0, 2.0, tanhCoefficients_m, derivativesPos);
    EXPECT_NEAR(derivativesNeg[0], 0.49966464986953352, 1e-10);
    EXPECT_NEAR(derivativesPos[0], 0.49966464986953352, 1e-10);
    EXPECT_NEAR(derivativesNeg[1], 0.49932952465848707, 1e-10);
    EXPECT_NEAR(derivativesPos[1], -0.49932952465848707, 1e-10);
    EXPECT_NEAR(derivativesNeg[2], -0.0013400513070529474, 1e-10);
    EXPECT_NEAR(derivativesPos[2], -0.0013400513070529474, 1e-10);
    EXPECT_NEAR(derivativesNeg[3], -1.0026765069198489, 1e-10);
    EXPECT_NEAR(derivativesPos[3], 1.0026765069198489, 1e-10);
    EXPECT_NEAR(derivativesNeg[4], -0.0053386419156264964, 1e-10);
    EXPECT_NEAR(derivativesPos[4], -0.0053386419156264964, 1e-10);
    EXPECT_NEAR(derivativesNeg[5], 7.9893801387861263, 1e-10);
    EXPECT_NEAR(derivativesPos[5], -7.9893801387861263, 1e-10);
    EXPECT_NEAR(derivativesNeg[6], -0.021010422121266359, 1e-10);
    EXPECT_NEAR(derivativesPos[6], -0.021010422121266359, 1e-10);
}
