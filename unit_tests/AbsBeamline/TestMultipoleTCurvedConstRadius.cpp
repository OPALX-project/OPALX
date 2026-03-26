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
#include <vector>
#include "AbsBeamline/Component.h"
#include "AbsBeamline/MultipoleT.h"
#include "gtest/gtest.h"

class TestMultipoleTCurvedConstRadius : public testing::Test, public MultipoleT {
public:
    TestMultipoleTCurvedConstRadius() : MultipoleT("Magnet") {}

protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }
    static void TearDownTestSuite() { ippl::finalize(); }

    // Test helper functions
    static Vector_t<double, 3> curvilinearToGlobal(
        const Vector_t<double, 3>& local, const Vector_t<double, 3>& elementEntry,
        const double elementLength, const double bendAngle) {
        const auto radius  = elementLength / bendAngle;
        const double s     = local[2] + elementLength / 2.0;
        const auto cosSbyR = std::cos(s / radius);
        const auto sinSbyR = std::sin(s / radius);
        const double x     = radius - (local[0] + radius) * cosSbyR + elementEntry[0];
        const double y     = local[1] + elementEntry[1];
        const double z     = (local[0] + radius) * sinSbyR + elementEntry[2];
        return {x, y, z};
    }

    void grabTransverseDataLine(
        std::vector<double>& line, const double s, const double width,
        const Vector_t<double, 3>& elementEntry, const double elementLength,
        const double bendAngle) const {
        // Create the views
        std::vector<Vector_t<double, 3>> local;
        Kokkos::View<Vector_t<double, 3>*> R;
        Kokkos::View<Vector_t<double, 3>*> E;
        Kokkos::View<Vector_t<double, 3>*> B;
        local.resize(line.size());
        Kokkos::resize(R, line.size());
        Kokkos::resize(E, line.size());
        Kokkos::resize(B, line.size());
        const auto hostR = Kokkos::create_mirror_view(R);
        const auto hostB = Kokkos::create_mirror_view(B);
        // Set the particle positions
        const double stepSize = width / static_cast<double>(line.size() - 1);
        for (size_t i = 0; i < line.size(); ++i) {
            local[i] = {static_cast<double>(i) * stepSize - width / 2, 0, s};
            hostR(i) = curvilinearToGlobal(local[i], elementEntry, elementLength, bendAngle);
        }
        Kokkos::deep_copy(R, hostR);
        // Get the fields
        apply(R, E, B, 0.0, line.size());
        // Return the fields
        Kokkos::deep_copy(hostB, B);
        for (size_t i = 0; i < line.size(); ++i) {
            line[i] = std::hypot(hostB(i)[0], hostB(i)[1], hostB(i)[2]);
            //std::cout << i << ": Local=" << local[i] << ", Global=" << R[i]
            //          << ", mag(B)=" << line[i] << std::endl;
        }
        //std::cout << std::endl;
    }

    void grabLongitudinalDivCurlLine(
        std::vector<double>& fieldLine, std::vector<double>& divLine,
        std::vector<Vector_t<double, 3>>& curlLine, const double x, const double length,
        const Vector_t<double, 3>& elementEntry, const double elementLength, const double bendAngle,
        const double dr) {
        const double stepSize = length / static_cast<double>(divLine.size() - 1);
        const double startS   = 0 - (elementLength - length) / 2;
        for (size_t i = 0; i < divLine.size(); ++i) {
            // Get the sourrounding 6 B fields
            Vector_t<double, 3> local{x, 0, static_cast<double>(i) * stepSize + startS};
            Vector_t<double, 3> R =
                curvilinearToGlobal(local, elementEntry, elementLength, bendAngle);
            Vector_t<double, 3> B;
            Vector_t<double, 3> Bxp;
            Vector_t<double, 3> Bxm;
            Vector_t<double, 3> Byp;
            Vector_t<double, 3> Bym;
            Vector_t<double, 3> Bzp;
            Vector_t<double, 3> Bzm;
            Vector_t<double, 3> E;
            apply(R, {}, 0.0, E, B);
            apply(R + Vector_t<double, 3>{dr, 0, 0}, {}, 0.0, E, Bxp);
            apply(R - Vector_t<double, 3>{dr, 0, 0}, {}, 0.0, E, Bxm);
            apply(R + Vector_t<double, 3>{0, dr, 0}, {}, 0.0, E, Byp);
            apply(R - Vector_t<double, 3>{0, dr, 0}, {}, 0.0, E, Bym);
            apply(R + Vector_t<double, 3>{0, 0, dr}, {}, 0.0, E, Bzp);
            apply(R - Vector_t<double, 3>{0, 0, dr}, {}, 0.0, E, Bzm);
            // Calculate the divergence and curl using central differences
            divLine[i]     = (Bxp[0] - Bxm[0] + Byp[1] - Bym[1] + Bzp[2] - Bzm[2]) / 2 / dr;
            curlLine[i][0] = (Byp[2] - Bym[2] - Bzp[1] + Bzm[1]) / 2 / dr;
            curlLine[i][1] = (Bzp[0] - Bzm[0] - Bxp[2] + Bxm[2]) / 2 / dr;
            curlLine[i][2] = (Bxp[1] - Bxm[1] - Byp[0] + Bym[0]) / 2 / dr;
            fieldLine[i]   = std::hypot(B[0], B[1], B[2]);
            //std::cout << i << ": Local=" << local << ", Global=" << R << ", div(B)=" << divLine[i]
            //          << ", curl(B)=" << curlLine[i] << ", B=" << fieldLine[i] << std::endl;
        }
    }
};

TEST_F(TestMultipoleTCurvedConstRadius, Dipole) {
    std::vector<double> line(101);
    // Set up the magnet
    constexpr double length      = 4.4;
    constexpr double bendAngle   = M_PI / 8.0;
    constexpr double dipoleField = 1.0;
    setBendAngle(bendAngle, false);
    setElementLength(length);
    setAperture(3.5, 3.5);
    setFringeField(length / 2, 0.3, 0.3);
    setRotation(0.0);
    setEntranceAngle(0.0);
    setMaxOrder(5, 10);
    // Check dipole has the correct constant transverse field magnitude at the center
    setTransProfile({dipoleField});
    grabTransverseDataLine(line, 0, 3.0, {0, 0, 0}, length, bendAngle);
    for (const double val : line) {
        EXPECT_NEAR(val, dipoleField, 1e-2);
    }
    // Check dipole has the correct constant transverse field magnitude at the edge
    setTransProfile({dipoleField});
    grabTransverseDataLine(line, -length / 2, 3.0, {0, 0, 0}, length, bendAngle);
    for (const double val : line) {
        EXPECT_NEAR(val, dipoleField / 2, 1e-2);
    }
}

TEST_F(TestMultipoleTCurvedConstRadius, Quadrupole) {
    constexpr unsigned int samplesPerSide = 50;
    std::vector<double> line(2 * samplesPerSide + 1);
    // Set up the magnet
    constexpr double length          = 4.4;
    constexpr double bendAngle       = M_PI / 8.0;
    constexpr double quadrupoleField = 1.0;
    setBendAngle(bendAngle, false);
    setElementLength(length);
    setAperture(3.5, 3.5);
    setFringeField(length / 2, 0.3, 0.3);
    setRotation(0.0);
    setEntranceAngle(0.0);
    setMaxOrder(5, 10);
    // Check quadrupole has linear field magnitude at the center
    setTransProfile({0, quadrupoleField});
    grabTransverseDataLine(line, 0, 3.0, {0, 0, 0}, length, bendAngle);
    double expected = 0;
    double delta    = line[samplesPerSide + 1] - line[samplesPerSide];
    EXPECT_NEAR(delta, 0.03, 1e-5);
    EXPECT_NEAR(line[samplesPerSide], expected, 1e-5);
    for (size_t i = 0; i < samplesPerSide; ++i) {
        expected += delta;
        EXPECT_NEAR(line[samplesPerSide - i - 1], expected, 1e-2);
        EXPECT_NEAR(line[samplesPerSide + i + 1], expected, 1e-2);
    }
    // Check quadrupole has linear field magnitude at the edge with half the slope
    setTransProfile({0, quadrupoleField});
    grabTransverseDataLine(line, -length / 2, 3.0, {0, 0, 0}, length, bendAngle);
    expected = 0;
    delta    = line[samplesPerSide + 1] - line[samplesPerSide];
    EXPECT_NEAR(delta, 0.015, 1e-5);
    EXPECT_NEAR(line[samplesPerSide], expected, 1e-5);
    for (size_t i = 0; i < samplesPerSide; ++i) {
        expected += delta;
        EXPECT_NEAR(line[samplesPerSide - i - 1], expected, 1e-2);
        EXPECT_NEAR(line[samplesPerSide + i + 1], expected, 1e-2);
    }
}

TEST_F(TestMultipoleTCurvedConstRadius, Sextupole) {
    constexpr unsigned int samplesPerSide = 50;
    std::vector<double> line(2 * samplesPerSide + 1);
    // Set up the magnet
    constexpr double length         = 4.4;
    constexpr double bendAngle      = M_PI / 8.0;
    constexpr double sextupoleField = 1.0;
    setBendAngle(bendAngle, false);
    setElementLength(length);
    setAperture(3.5, 3.5);
    setFringeField(length / 2, 0.3, 0.3);
    setRotation(0.0);
    setEntranceAngle(0.0);
    setMaxOrder(5, 10);
    // Check sextupole has quadratic field magnitude at the center
    setTransProfile({0, 0, sextupoleField});
    grabTransverseDataLine(line, 0, 3.0, {0, 0, 0}, length, bendAngle);
    double expected = 0;
    double delta    = std::sqrt(line[samplesPerSide + 1] - line[samplesPerSide]);
    EXPECT_NEAR(delta, 0.03, 1e-5);
    EXPECT_NEAR(line[samplesPerSide], expected, 1e-5);
    for (size_t i = 0; i < samplesPerSide; ++i) {
        expected = std::pow(static_cast<double>(i + 1) * delta, 2);
        EXPECT_NEAR(line[samplesPerSide - i - 1], expected, 1e-2);
        EXPECT_NEAR(line[samplesPerSide + i + 1], expected, 1e-2);
    }
    // Check sextupole has quadratic field magnitude at the edge with half the slope
    setTransProfile({0, 0, sextupoleField});
    grabTransverseDataLine(line, -length / 2, 3.0, {0, 0, 0}, length, bendAngle);
    expected = 0;
    delta    = std::sqrt(line[samplesPerSide + 1] - line[samplesPerSide]);
    EXPECT_NEAR(delta, 0.02121, 1e-5);
    EXPECT_NEAR(line[samplesPerSide], expected, 1e-5);
    for (size_t i = 0; i < samplesPerSide; ++i) {
        expected = std::pow(static_cast<double>(i + 1) * delta, 2);
        EXPECT_NEAR(line[samplesPerSide - i - 1], expected, 1e-2);
        EXPECT_NEAR(line[samplesPerSide + i + 1], expected, 1e-2);
    }
}

TEST_F(TestMultipoleTCurvedConstRadius, Octupole) {
    constexpr unsigned int samplesPerSide = 50;
    std::vector<double> line(2 * samplesPerSide + 1);
    // Set up the magnet
    constexpr double length         = 4.4;
    constexpr double bendAngle      = M_PI / 8.0;
    constexpr double sextupoleField = 1.0;
    setBendAngle(bendAngle, false);
    setElementLength(length);
    setAperture(3.5, 3.5);
    setFringeField(length / 2, 0.3, 0.3);
    setRotation(0.0);
    setEntranceAngle(0.0);
    setMaxOrder(5, 10);
    // Check octupole has cubic field magnitude at the center
    setTransProfile({0, 0, 0, sextupoleField});
    grabTransverseDataLine(line, 0, 3.0, {0, 0, 0}, length, bendAngle);
    double expected = 0;
    double delta    = std::cbrt(line[samplesPerSide + 1] - line[samplesPerSide]);
    EXPECT_NEAR(delta, 0.03, 1e-5);
    EXPECT_NEAR(line[samplesPerSide], expected, 1e-5);
    for (size_t i = 0; i < samplesPerSide; ++i) {
        expected = std::pow(static_cast<double>(i + 1) * delta, 3);
        EXPECT_NEAR(line[samplesPerSide - i - 1], expected, 1e-2);
        EXPECT_NEAR(line[samplesPerSide + i + 1], expected, 1e-2);
    }
    // Check octupole has cubic field magnitude at the edge with half the slope
    setTransProfile({0, 0, 0, sextupoleField});
    grabTransverseDataLine(line, -length / 2, 3.0, {0, 0, 0}, length, bendAngle);
    expected = 0;
    delta    = std::cbrt(line[samplesPerSide + 1] - line[samplesPerSide]);
    EXPECT_NEAR(delta, 0.023811, 1e-5);
    EXPECT_NEAR(line[samplesPerSide], expected, 1e-5);
    for (size_t i = 0; i < samplesPerSide; ++i) {
        expected = std::pow(static_cast<double>(i + 1) * delta, 3);
        EXPECT_NEAR(line[samplesPerSide - i - 1], expected, 1e-2);
        EXPECT_NEAR(line[samplesPerSide + i + 1], expected, 1e-2);
    }
}

TEST_F(TestMultipoleTCurvedConstRadius, Decapole) {
    constexpr unsigned int samplesPerSide = 50;
    std::vector<double> line(2 * samplesPerSide + 1);
    // Set up the magnet
    constexpr double length        = 4.4;
    constexpr double bendAngle     = M_PI / 8.0;
    constexpr double decapoleField = 1.0;
    setBendAngle(bendAngle, false);
    setElementLength(length);
    setAperture(3.5, 3.5);
    setFringeField(length / 2, 0.3, 0.3);
    setRotation(0.0);
    setEntranceAngle(0.0);
    setMaxOrder(5, 10);
    // Check decapole has quartic field magnitude at the center
    setTransProfile({0, 0, 0, 0, decapoleField});
    grabTransverseDataLine(line, 0, 3.0, {0, 0, 0}, length, bendAngle);
    double expected = 0;
    double delta    = std::pow(line[samplesPerSide + 1] - line[samplesPerSide], 1.0 / 4.0);
    EXPECT_NEAR(delta, 0.03, 1e-5);
    EXPECT_NEAR(line[samplesPerSide], expected, 1e-5);
    for (size_t i = 0; i < samplesPerSide; ++i) {
        expected = std::pow(static_cast<double>(i + 1) * delta, 4);
        EXPECT_NEAR(line[samplesPerSide - i - 1], expected, 1e-2);
        EXPECT_NEAR(line[samplesPerSide + i + 1], expected, 1e-2);
    }
    // Check decapole has quartic field magnitude at the edge with half the slope
    setTransProfile({0, 0, 0, 0, decapoleField});
    grabTransverseDataLine(line, -length / 2, 3.0, {0, 0, 0}, length, bendAngle);
    expected = 0;
    delta    = std::pow(line[samplesPerSide + 1] - line[samplesPerSide], 1.0 / 4.0);
    EXPECT_NEAR(delta, 0.025227, 1e-5);
    EXPECT_NEAR(line[samplesPerSide], expected, 1e-5);
    for (size_t i = 0; i < samplesPerSide; ++i) {
        expected = std::pow(static_cast<double>(i + 1) * delta, 4);
        EXPECT_NEAR(line[samplesPerSide - i - 1], expected, 1e-2);
        EXPECT_NEAR(line[samplesPerSide + i + 1], expected, 1e-2);
    }
}

TEST_F(TestMultipoleTCurvedConstRadius, DivCurl) {
    constexpr unsigned int samplesPerSide = 20;
    std::vector<double> divLine(2 * samplesPerSide + 1);
    std::vector<double> fieldLine(2 * samplesPerSide + 1);
    std::vector<Vector_t<double, 3>> curlLine(2 * samplesPerSide + 1);
    constexpr double dr = 0.001;
    // Set up the magnet
    constexpr double length    = 4.4;
    constexpr double bendAngle = M_PI / 8.0;
    setBendAngle(bendAngle, false);
    setElementLength(length);
    setAperture(3.5, 3.5);
    setFringeField(length / 2, 0.3, 0.3);
    setRotation(0.0);
    setEntranceAngle(0.0);
    setMaxOrder(5, 10);
    setTransProfile({1, 1, 1, 1, 1});
    // Get the data
    grabLongitudinalDivCurlLine(
        fieldLine, divLine, curlLine, 0, length * 1.5, {}, length, bendAngle, dr);
    // Analyse the results
    for (size_t i = 0; i < divLine.size(); ++i) {
        const double divError = std::abs(divLine[i]) * dr / fieldLine[i];
        const double curlError =
            std::hypot(curlLine[i][0], curlLine[i][1], curlLine[i][2]) * dr / fieldLine[i];
        EXPECT_NEAR(divError, 0, 1e-5);
        EXPECT_NEAR(curlError, 0, 1e-2);
    }
}

TEST_F(TestMultipoleTCurvedConstRadius, Geometry) {
    TestMultipoleTCurvedConstRadius* magnet = this;
    auto* geom = &magnet->getGeometry();
    EXPECT_NE(geom, nullptr);
    auto* constGeom = &const_cast<const TestMultipoleTCurvedConstRadius*>(magnet)->getGeometry();
    EXPECT_NE(constGeom, nullptr);
    EXPECT_EQ(geom, constGeom);
}

TEST_F(TestMultipoleTCurvedConstRadius, BoundingBox) {
    std::vector<double> line(3);
    // Set up the magnet
    constexpr double length      = 4;
    constexpr double bendAngle   = M_PI / 8.0;
    constexpr double dipoleField = 1.0;
    setBendAngle(bendAngle, false);
    setElementLength(length);
    setAperture(3.5, 3.5);
    setFringeField(length / 2, 3, 3);
    setRotation(0.0);
    setEntranceAngle(0.0);
    setMaxOrder(5, 10);
    setBoundingBoxLength(6);
    setTransProfile({dipoleField});
    // Check field vanishes outside the bounding box
    grabTransverseDataLine(line, -3.1, 3.0, {0, 0, 0}, length, bendAngle);
    for (const double val : line) {
        EXPECT_EQ(val, 0.0);
    }
    // Check field is present inside the bounding box
    grabTransverseDataLine(line, -3.0, 3.0, {0, 0, 0}, length, bendAngle);
    for (const double val : line) {
        EXPECT_NE(val, 0.0);
    }
    // Check field is present inside the bounding box
    grabTransverseDataLine(line, 0.0, 3.0, {0, 0, 0}, length, bendAngle);
    for (const double val : line) {
        EXPECT_NE(val, 0.0);
    }
    // Check field is present inside the bounding box
    grabTransverseDataLine(line, 3.0, 3.0, {0, 0, 0}, length, bendAngle);
    for (const double val : line) {
        EXPECT_NE(val, 0.0);
    }
    // Check field vanishes outside the bounding box
    grabTransverseDataLine(line, 3.1, 3.0, {0, 0, 0}, length, bendAngle);
    for (const double val : line) {
        EXPECT_EQ(val, 0.0);
    }
}

TEST_F(TestMultipoleTCurvedConstRadius, Aperture) {
    std::vector<double> line(9);
    // Set up the magnet
    constexpr double length      = 4;
    constexpr double bendAngle   = M_PI / 8.0;
    constexpr double dipoleField = 1.0;
    setBendAngle(bendAngle, false);
    setElementLength(length);
    setAperture(3.5, 3.5);
    setFringeField(length / 2, 3, 3);
    setRotation(0.0);
    setEntranceAngle(0.0);
    setMaxOrder(5, 10);
    setTransProfile({dipoleField});
    // Check field vanishes outside the aperture
    grabTransverseDataLine(line, 0, 4.0, {0, 0, 0}, length, bendAngle);
    EXPECT_EQ(line[0], 0.0);  // -4.0
    EXPECT_NE(line[1], 0.0);  // -3.0
    EXPECT_NE(line[2], 0.0);  // -2.0
    EXPECT_NE(line[3], 0.0);  // -1.0
    EXPECT_NE(line[4], 0.0);  // 0.0
    EXPECT_NE(line[5], 0.0);  // 1.0
    EXPECT_NE(line[6], 0.0);  // 2.0
    EXPECT_NE(line[7], 0.0);  // 3.0
    EXPECT_EQ(line[8], 0.0);  // 3.0
}

