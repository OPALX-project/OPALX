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
#include "gtest/gtest.h"
#include "AbsBeamline/MultipoleT.h"
#include <vector>

#include "AbsBeamline/Component.h"

class TestMultipoleTStraight : public testing::Test, public MultipoleT {
public:
    TestMultipoleTStraight() : MultipoleT("Magnet") {
    }

protected:
    // Overrides of testing::Test
    static void SetUpTestSuite() {
        int argc    = 0;
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

    // Test helper functions
    static Vector_t<double, 3> rotateBy(const Vector_t<double, 3>& center,
            const Vector_t<double, 3>& point, const double theta) {
        const auto x_prime = point[0] - center[0];
        const auto z_prime = point[2] - center[2];
        const auto x_rotated = x_prime * cos(theta) - z_prime * sin(theta) + center[0];
        const auto z_rotated = x_prime * sin(theta) + z_prime * cos(theta) + center[2];
        return {x_rotated, point[1], z_rotated};
    }

    static std::vector<std::vector<double>> partialsDerivB(const Vector_t<double, 3>& R,
        const Vector_t<double, 3>& /*B*/, const double stepSize, Component* dummyField,
        const double theta = 0.0)
    {
        // builds a matrix of all partial derivatives of B -> dx_i B_j
        std::vector allPartials(3, std::vector<double>(3));
        for(int i = 0; i < 3; i++) {
            Vector_t<double, 3> P, E;
            double t = 0;
            // B at the previous and next grid points R_prev, R_next
            Vector_t<double, 3> R_pprev = R, R_prev = R, R_next = R, R_nnext = R;
            R_pprev(i) -= 2 * stepSize;
            R_nnext(i) += 2 * stepSize;
            R_prev(i) -= stepSize;
            R_next(i) += stepSize;
            // Rotate these points by the angle
            R_pprev = rotateBy(R, R_pprev, theta);
            R_nnext = rotateBy(R, R_nnext, theta);
            R_prev = rotateBy(R, R_prev, theta);
            R_next = rotateBy(R, R_next, theta);
            // Get the magnetic fields and derivatives
            Vector_t<double, 3> B_prev, B_next, B_pprev, B_nnext;
            dummyField->apply(R_prev, P, t, E, B_prev);
            dummyField->apply(R_next, P, t, E, B_next);
            dummyField->apply(R_pprev, P, t, E, B_pprev);
            dummyField->apply(R_nnext, P, t, E, B_nnext);
            for(int j = 0; j < 3; j++) {
                allPartials[i][j] =
                        (B_pprev[j] - 8 * B_prev[j] + 8 * B_next[j] - B_nnext[j]) / (12 * stepSize);
            }
        }
        return allPartials;
    }

    static double calcDivB(const Vector_t<double, 3>& R, const Vector_t<double, 3>& B,
            const double stepSize, Component* dummyField, const double theta = 0.0) {
        double div = 0;
        std::vector partials(3, std::vector<double>(3));
        partials = partialsDerivB(R, B, stepSize, dummyField, theta);
        for(int i = 0; i < 3; i++) {
            div += partials[i][i];
        }
        return div;
    }

    static std::vector<double> calcCurlB(const Vector_t<double, 3>& R, const Vector_t<double, 3>& B,
            const double stepSize, Component* dummyField, const double theta = 0.0) {
        std::vector<double> curl(3);
        std::vector partials(3, std::vector<double>(3));
        partials = partialsDerivB(R, B, stepSize, dummyField, theta);
        curl[0] = partials[1][2] - partials[2][1];
        curl[1] = partials[2][0] - partials[0][2];
        curl[2] = partials[0][1] - partials[1][0];
        return curl;
    }

    void grabDataLine(std::vector<double>& line, Vector_t<double, 3> pos) {
        constexpr double stepSize = 3.0/100.0;
        const Vector_t<double, 3> P(3);
        Vector_t<double, 3> E(3);
        for(size_t i = 0; i < line.size(); ++i) {
            double t = 0.0;
            Vector_t<double, 3> R{pos[0], pos[1] - 1.5 + static_cast<double>(i) * stepSize, pos[2]};
            Vector_t<double, 3> B{};
            apply(R, P, t, E, B);
            line[i] = std::hypot(B[0], B[1], B[2]);
            std::cout << line[i] << " ";
        }
        std::cout << std::endl;
    }

    void grabDataLineParallel(std::vector<double>& line, Vector_t<double, 3> pos) {
        constexpr double stepSize = 3.0/100.0;
        // Create the views
        Kokkos::View<Vector_t<double, 3>*> R;
        Kokkos::View<Vector_t<double, 3>*> E;
        Kokkos::View<Vector_t<double, 3>*> B;
        Kokkos::resize(R, line.size());
        Kokkos::resize(E, line.size());
        Kokkos::resize(B, line.size());
        auto hostR = Kokkos::create_mirror_view(R);
        auto hostB = Kokkos::create_mirror_view(B);
        // Set the particle positions
        for(size_t i = 0; i < line.size(); ++i) {
            hostR(i) = {pos[0], pos[1] - 1.5 + static_cast<double>(i) * stepSize, pos[2]};
        }
        Kokkos::deep_copy(R, hostR);
        // Get the fields
        apply(R, E, B, 0.0);
        // Return the fields
        Kokkos::deep_copy(hostB, B);
        for(size_t i = 0; i < line.size(); ++i) {
            line[i] = std::hypot(hostB(i)[0], hostB(i)[1], hostB(i)[2]);
            std::cout << line[i] << " ";
        }
        std::cout << std::endl;
    }
};

TEST_F(TestMultipoleTStraight, StraightShape) {
    std::vector<double> line;
    line.resize(101);
    // Set up the magnet
    constexpr double length = 4.4;
    setBendAngle(0.0, false);
    setElementLength(length);
    setAperture(3.5, 3.5);
    setFringeField(2.2, 0.3, 0.3);
    setRotation(0.0);
    setEntranceAngle(0.0);
    setMaxOrder(5, 10);
    const auto pos = localCartesianToOpalCartesian({0,0,0});
    // Check dipole has constant field magnitude
    setTransProfile({1.0});
    grabDataLine(line, pos);
    for(size_t i = 0; i < line.size(); ++i) {
        EXPECT_NEAR(line[i], 1.0, 1e-2);
    }
    // Check quadrupole has linear field magnitude
    setTransProfile({0.0, 1.0});
    grabDataLine(line, pos);
    double expected = 0;
    double delta = line[51] - line[50];
    EXPECT_NEAR(line[50], expected, 1e-2);
    for(size_t i = 0; i < 50; ++i) {
        expected += delta;
        EXPECT_NEAR(line[50 - i - 1], expected, 1e-2);
        EXPECT_NEAR(line[50 + i + 1], expected, 1e-2);
    }
    // Check the sextupole has quadratic field magnitude
    setTransProfile({0.0, 0.0, 1.0});
    grabDataLine(line, pos);
    delta = std::sqrt(line[51] - line[50]);
    EXPECT_NEAR(line[50], 0, 1e-2);
    for(size_t i = 0; i < 50; ++i) {
        expected = std::pow(static_cast<double>(i + 1) * delta, 2);
        EXPECT_NEAR(line[50 - i - 1], expected, 1e-2) << i;
        EXPECT_NEAR(line[50 + i + 1], expected, 1e-2) << i;
    }
    // Check the octupole has cubic field magnitude
    setTransProfile({0.0, 0.0, 0.0, 1.0});
    grabDataLine(line, pos);
    delta = std::cbrt(line[51] - line[50]);
    EXPECT_NEAR(line[50], 0, 1e-2);
    for(size_t i = 0; i < 50; ++i) {
        expected = std::pow(static_cast<double>(i + 1) * delta, 3);
        EXPECT_NEAR(line[50 - i - 1], expected, 1e-2) << i;
        EXPECT_NEAR(line[50 + i + 1], expected, 1e-2) << i;
    }
}

TEST_F(TestMultipoleTStraight, StraightShapeParallel) {
    std::vector<double> line;
    line.resize(101);
    // Set up the magnet
    constexpr double length = 4.4;
    setBendAngle(0.0, false);
    setElementLength(length);
    setAperture(3.5, 3.5);
    setFringeField(2.2, 0.3, 0.3);
    setRotation(0.0);
    setEntranceAngle(0.0);
    setMaxOrder(5, 10);
    const auto pos = localCartesianToOpalCartesian({0,0,0});
    // Check dipole has constant field magnitude
    setTransProfile({1.0});
    grabDataLineParallel(line, pos);
    for(size_t i = 0; i < line.size(); ++i) {
        EXPECT_NEAR(line[i], 1.0, 1e-2);
    }
    // Check quadrupole has linear field magnitude
    setTransProfile({0.0, 1.0});
    grabDataLineParallel(line, pos);
    double expected = 0;
    double delta = line[51] - line[50];
    EXPECT_NEAR(line[50], expected, 1e-2);
    for(size_t i = 0; i < 50; ++i) {
        expected += delta;
        EXPECT_NEAR(line[50 - i - 1], expected, 1e-2);
        EXPECT_NEAR(line[50 + i + 1], expected, 1e-2);
    }
    // Check the sextupole has quadratic field magnitude
    setTransProfile({0.0, 0.0, 1.0});
    grabDataLineParallel(line, pos);
    delta = std::sqrt(line[51] - line[50]);
    EXPECT_NEAR(line[50], 0, 1e-2);
    for(size_t i = 0; i < 50; ++i) {
        expected = std::pow(static_cast<double>(i + 1) * delta, 2);
        EXPECT_NEAR(line[50 - i - 1], expected, 1e-2) << i;
        EXPECT_NEAR(line[50 + i + 1], expected, 1e-2) << i;
    }
    // Check the octupole has cubic field magnitude
    setTransProfile({0.0, 0.0, 0.0, 1.0});
    grabDataLineParallel(line, pos);
    delta = std::cbrt(line[51] - line[50]);
    EXPECT_NEAR(line[50], 0, 1e-2);
    for(size_t i = 0; i < 50; ++i) {
        expected = std::pow(static_cast<double>(i + 1) * delta, 3);
        EXPECT_NEAR(line[50 - i - 1], expected, 1e-2) << i;
        EXPECT_NEAR(line[50 + i + 1], expected, 1e-2) << i;
    }
}

#if 0
TEST_F(TestMultipoleTStraight, Straight) {
    constexpr double length = 4.4;
    setBendAngle(0.0, false);
    setElementLength(length);
    setAperture(3.5, 3.5);
    setFringeField(2.2, 0.3, 0.3);
    setRotation(0.0);
    setEntranceAngle(0.0);
    setTransProfile({1.0, 1.0});
    setMaxOrder(5, 20);
    Vector_t<double, 3> R(0.0, 0.0, 0.0), E(3);
    const Vector_t<double, 3> P(3);
    for(size_t i = 0; i < 6; ++i) {
        double x = -0.3 + static_cast<double>(i) * 0.1;
        for(size_t j = 0; j < 6; ++j) {
            double z = -0.3;
            double t = 0.0;
            constexpr double stepSize = 1e-3;
            double y = -3.0 + static_cast<double>(j) * 1.0;
            R[0] = x;
            R[1] = z;
            R[2] = y - length / 2.0;
            Vector_t<double, 3> B(0., 0., 0.);
            apply(R, P, t, E, B);
            std::cout << "R: " << R[0] << " " << R[1] << " " << R[2] << " = "
                      << B[0] << " " << B[1] << " " << B[2] << std::endl;
            double div = calcDivB(R, B, stepSize, this);
            std::vector<double> curl = calcCurlB(R, B, stepSize, this);
            double curlMag = 0.0;
            curlMag += MultipoleTBase::powerInteger(curl[0], 2);
            curlMag += MultipoleTBase::powerInteger(curl[1], 2);
            curlMag += MultipoleTBase::powerInteger(curl[2], 2);
            curlMag = sqrt(curlMag);
            EXPECT_NEAR(div, 0, 1e-1)
                << "R: " << x << " " << z << " " << y << std::endl
                << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                << "Del: " << div << " " << curl[0] << " " << curl[1]
                << " " << curl[2] << std::endl;
            EXPECT_NEAR(curlMag, 0, 1e-1)
                << "R: " << x << " " << z << " " << y << std::endl
                << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                << "Del: " << div << " " << curl[0] << " " << curl[1] << " "
                << curl[2] << std::endl;
        }
    }
}
#endif
