/**
 * \file TestAstra1DDynamic.cpp
 * \brief Unit tests for Astra1DDynamic fieldmap class
 *
 * Tests cover:
 *
 * Parsing & metadata:
 * - File parsing
 * - Frequency parsing and conversion (MHz -> angular frequency)
 * - Field dimension retrieval (z-bounds)
 *
 * Field evaluation:
 * - Point field evaluation via getFieldstrength()
 * - Uniform on-axis field correctness
 * - Static computeField() correctness for simple coefficients
 * - Reconstructed on-axis field via getOnaxisEz()
 * - Device-safe traveling-wave application via applyTravelingWave()
 *
 * Geometry & bounds:
 * - Behaviour outside z-range
 * - isInside() helper (boundary semantics)
 *
 * API behaviour:
 * - getFieldDerivative() for simple constant field
 * - getFrequency() correctness
 * - getInfo() does not crash
 * - swap() does not crash
 *
 * Robustness:
 * - Missing / invalid file handling
 * - Fieldmap dictionary caching
 * - readMap() / freeMap() lifecycle correctness
 *
 * Semantics:
 * - Field accumulation: E and B are accumulated, not overwritten
 *
 *  * Notes:
 * - The current Astra1DDynamic preprocessing / Fourier reconstruction yields
 *   an amplitude factor of (nz - 1) / nz for uniform input fields. The tests
 *   reflect the current implementation behaviour.
 */

#include <gtest/gtest.h>

#include "Fields/Astra1DDynamic.h"
#include "Fields/Fieldmap.h"
#include "Ippl.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Utilities/GeneralOpalException.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace {

    // ---------------------------------------------------------------------------
    // Helper: write 1D Astra dynamic fieldmap
    // Format:
    // 1DDynamic accuracy [TRUE|FALSE]
    // freq_MHz
    // z Ez
    // z Ez
    // ...
    //
    // z is written in meters.
    // Ez is written in MV/m.
    // ---------------------------------------------------------------------------
    std::string writeAstra1DFieldmap(
            const std::string& path, const std::vector<double>& z_m,
            const std::vector<double>& ez_MVpm, int accuracy = 8, double freq_MHz = 100.0,
            bool normalize = true) {
        EXPECT_EQ(z_m.size(), ez_MVpm.size());

        std::ofstream f(path);
        f << "AstraDynamic " << accuracy << " " << (normalize ? "TRUE" : "FALSE") << "\n";
        f << freq_MHz << "\n";

        for (size_t i = 0; i < z_m.size(); ++i) {
            f << z_m[i] << " " << ez_MVpm[i] << "\n";
        }

        return path;
    }

    // ---------------------------------------------------------------------------
    // Helper: constant 1D map
    // ---------------------------------------------------------------------------
    std::string writeConstantAstra1DFieldmap(
            const std::string& path, double zbegin_m, double zend_m, int nz, double ez_MVpm,
            int accuracy = 8, double freq_MHz = 100.0, bool normalize = true) {
        std::vector<double> z(nz);
        std::vector<double> ez(nz, ez_MVpm);

        const double dz = (zend_m - zbegin_m) / (nz - 1);
        for (int i = 0; i < nz; ++i) {
            z[i] = zbegin_m + i * dz;
        }

        return writeAstra1DFieldmap(path, z, ez, accuracy, freq_MHz, normalize);
    }

}  // namespace

// ===========================================================================
// Fixture
// ===========================================================================
class Astra1DDynamicTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() {
        Fieldmap::clearDictionary();
        ippl::finalize();
    }

    void SetUp() override {
        tmpDir_ = std::filesystem::temp_directory_path() / "opalx_astra1d_test";
        std::filesystem::create_directories(tmpDir_);
    }

    void TearDown() override {
        Fieldmap::clearDictionary();
        std::filesystem::remove_all(tmpDir_);
    }

    std::string tmpFile(const std::string& name) const { return (tmpDir_ / name).string(); }

    std::filesystem::path tmpDir_;
};

// ===========================================================================
// Test: Parse fieldmap and verify field dimensions
// ===========================================================================
TEST_F(Astra1DDynamicTest, ParseAndDimensions) {
    const double zb   = 0.0;
    const double ze   = 0.10;
    const double freq = 100.0;
    const int nz      = 5;

    std::string fname =
            writeConstantAstra1DFieldmap(tmpFile("parse.map"), zb, ze, nz, 1.0, 8, freq, true);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);
    EXPECT_EQ(fm->getType(), TAstraDynamic);

    Fieldmap::readMap(fname);

    double zBegin = 0.0, zEnd = 0.0;
    fm->getFieldDimensions(zBegin, zEnd);

    EXPECT_NEAR(zBegin, zb, 1e-12);
    EXPECT_NEAR(zEnd, ze, 1e-12);
}

// ===========================================================================
// Test: Frequency parsing and conversion
// ===========================================================================
TEST_F(Astra1DDynamicTest, FrequencyParsing) {
    const double freq = 100.0;

    std::string fname =
            writeConstantAstra1DFieldmap(tmpFile("freq.map"), 0.0, 0.10, 5, 1.0, 8, freq, true);

    auto* fm = dynamic_cast<Astra1DDynamic*>(Fieldmap::getFieldmap(fname));
    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    const double expected = freq * Physics::two_pi * Units::MHz2Hz;
    EXPECT_NEAR(fm->getFrequency(), expected, 1e-6);
}

// ===========================================================================
// Test: isInside() helper
// ===========================================================================
TEST_F(Astra1DDynamicTest, IsInside) {
    std::string fname =
            writeConstantAstra1DFieldmap(tmpFile("inside.map"), 0.0, 0.10, 5, 1.0, 8, 100.0, true);

    auto* fm = dynamic_cast<Astra1DDynamic*>(Fieldmap::getFieldmap(fname));
    ASSERT_NE(fm, nullptr);

    EXPECT_TRUE(fm->isInside({0.0, 0.0, 0.00}));
    EXPECT_TRUE(fm->isInside({0.0, 0.0, 0.05}));
    EXPECT_FALSE(fm->isInside({0.0, 0.0, 0.10}));
    EXPECT_FALSE(fm->isInside({0.0, 0.0, -0.01}));
}

// ===========================================================================
// Test: Uniform field – on-axis field strength
// ===========================================================================
TEST_F(Astra1DDynamicTest, UniformFieldStrengthOnAxis) {
    // With normalization = true and uniform Ez = 3 MV/m,
    // the reconstructed on-axis field follows the current Astra1DDynamic
    // Fourier/mirroring convention, which for nz samples gives a factor
    // (nz - 1) / nz on a uniform input field.

    const int nz = 9;

    std::string fname = writeConstantAstra1DFieldmap(
            tmpFile("uniform.map"), 0.0, 0.10, nz, 3.0, 8, 100.0, true);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool outside = fm->getFieldstrength(R, E, B);

    EXPECT_FALSE(outside);

    EXPECT_NEAR(E[0], 0.0, 1e-10);
    EXPECT_NEAR(E[1], 0.0, 1e-10);

    // const double expected_scale = static_cast<double>(nz - 1) / nz;
    EXPECT_NEAR(E[2], 1.0e6, 1e-4);

    // on-axis, transverse magnetic field projects to zero
    EXPECT_NEAR(B[0], 0.0, 1e-10);
    EXPECT_NEAR(B[1], 0.0, 1e-10);
    EXPECT_NEAR(B[2], 0.0, 1e-10);
}

// ===========================================================================
// Test: No normalization
// ===========================================================================
TEST_F(Astra1DDynamicTest, NoNormalization) {
    const int nz = 9;

    std::string fname = writeConstantAstra1DFieldmap(
            tmpFile("nonorm.map"), 0.0, 0.10, nz, 3.0, 8, 100.0, false);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool outside = fm->getFieldstrength(R, E, B);

    EXPECT_FALSE(outside);

    // const double expected_scale = static_cast<double>(nz - 1) / nz;
    EXPECT_NEAR(E[2], 3.0e6, 1e-3);
}

// ===========================================================================
// Test: Outside z range returns outside flag and no modification
// ===========================================================================
TEST_F(Astra1DDynamicTest, OutsideZRange) {
    std::string fname =
            writeConstantAstra1DFieldmap(tmpFile("outside.map"), 0.0, 0.10, 9, 1.0, 8, 100.0, true);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    {
        Vector_t<double, 3> R = {0.0, 0.0, -0.01};
        Vector_t<double, 3> E = {1.0, 2.0, 3.0};
        Vector_t<double, 3> B = {4.0, 5.0, 6.0};

        bool outside = fm->getFieldstrength(R, E, B);

        EXPECT_TRUE(outside);
        EXPECT_NEAR(E[0], 1.0, 1e-15);
        EXPECT_NEAR(E[1], 2.0, 1e-15);
        EXPECT_NEAR(E[2], 3.0, 1e-15);
        EXPECT_NEAR(B[0], 4.0, 1e-15);
        EXPECT_NEAR(B[1], 5.0, 1e-15);
        EXPECT_NEAR(B[2], 6.0, 1e-15);
    }

    {
        Vector_t<double, 3> R = {0.0, 0.0, 0.10};
        Vector_t<double, 3> E = {1.0, 2.0, 3.0};
        Vector_t<double, 3> B = {4.0, 5.0, 6.0};

        bool outside = fm->getFieldstrength(R, E, B);

        EXPECT_TRUE(outside);
        EXPECT_NEAR(E[0], 1.0, 1e-15);
        EXPECT_NEAR(E[1], 2.0, 1e-15);
        EXPECT_NEAR(E[2], 3.0, 1e-15);
        EXPECT_NEAR(B[0], 4.0, 1e-15);
        EXPECT_NEAR(B[1], 5.0, 1e-15);
        EXPECT_NEAR(B[2], 6.0, 1e-15);
    }
}

// ===========================================================================
// Test: Accumulation semantics – E and B are accumulated, not overwritten
// ===========================================================================
TEST_F(Astra1DDynamicTest, FieldAccumulation) {
    const int nz = 9;

    std::string fname =
            writeConstantAstra1DFieldmap(tmpFile("accum.map"), 0.0, 0.10, nz, 1.0, 8, 100.0, true);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {1.0, 2.0, 3.0};
    Vector_t<double, 3> B = {4.0, 5.0, 6.0};

    fm->getFieldstrength(R, E, B);

    EXPECT_NEAR(E[0], 1.0, 1e-10);
    EXPECT_NEAR(E[1], 2.0, 1e-10);

    // const double expected_scale = static_cast<double>(nz - 1) / nz;
    EXPECT_NEAR(E[2], 3.0 + 1.0e6, 1e-4);

    // on-axis, no transverse magnetic contribution
    EXPECT_NEAR(B[0], 4.0, 1e-10);
    EXPECT_NEAR(B[1], 5.0, 1e-10);
    EXPECT_NEAR(B[2], 6.0, 1e-10);
}

// ===========================================================================
// Test: getFieldDerivative for constant field should be zero on-axis
// ===========================================================================
TEST_F(Astra1DDynamicTest, FieldDerivativeConstantField) {
    std::string fname =
            writeConstantAstra1DFieldmap(tmpFile("deriv.map"), 0.0, 0.10, 9, 2.0, 8, 100.0, true);

    auto* fm = dynamic_cast<Astra1DDynamic*>(Fieldmap::getFieldmap(fname));
    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool outside = fm->getFieldDerivative(R, E, B, DZ);

    EXPECT_FALSE(outside);
    EXPECT_NEAR(E[2], 0.0, 1e-6);
}

// ===========================================================================
// Test: getOnaxisEz returns the normalized raw on-axis profile
// for the legacy 2-token Astra1D header format.
// ===========================================================================
TEST_F(Astra1DDynamicTest, GetOnaxisEzReturnsNormalizedRawProfile) {
    const int nz            = 9;
    const std::string fname = tmpFile("onaxis_legacy.map");

    {
        std::ofstream out(fname);
        ASSERT_TRUE(out.good());

        // Legacy 2-token header accepted by both:
        // - Fieldmap::getFieldmap()
        // - current Astra1DDynamic::getOnaxisEz()
        out << "AstraDynamic 8\n";
        out << "100.0\n";

        for (int i = 0; i < nz; ++i) {
            const double z = 0.10 * double(i) / double(nz - 1);
            out << std::setprecision(16) << z << " " << 3.0 << "\n";
        }
    }

    auto* fm = dynamic_cast<Astra1DDynamic*>(Fieldmap::getFieldmap(fname));
    ASSERT_NE(fm, nullptr);

    ASSERT_NO_THROW(Fieldmap::readMap(fname));

    std::vector<std::pair<double, double>> F;
    ASSERT_NO_THROW(fm->getOnaxisEz(F));

    ASSERT_EQ(F.size(), nz);

    EXPECT_NEAR(F.front().first, 0.0, 1e-12);
    EXPECT_NEAR(F.back().first, 0.10, 1e-12);

    for (int i = 0; i < nz; ++i) {
        EXPECT_NEAR(F[i].second, 1.0, 1e-12);
    }
}

// ===========================================================================
// Test: computeField static function directly with host view
// ===========================================================================
TEST_F(Astra1DDynamicTest, ComputeFieldStaticConstantCoefficient) {
    // accuracy = 2  => coefficients [a0, a1, b1]
    // choose only a0 = 1.0, so ez = 1.0 everywhere
    Kokkos::View<double*, Kokkos::HostSpace> coefs("coefs", 3);
    coefs(0) = 1.0;
    coefs(1) = 0.0;
    coefs(2) = 0.0;

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    Astra1DDynamic::computeField(
            R, E, B, coefs,
            0.0,   // zbegin
            0.20,  // length
            1.0,   // xlrep
            2      // accuracy
    );

    EXPECT_NEAR(E[0], 0.0, 1e-12);
    EXPECT_NEAR(E[1], 0.0, 1e-12);
    EXPECT_NEAR(E[2], 1.0, 1e-12);

    EXPECT_NEAR(B[0], 0.0, 1e-12);
    EXPECT_NEAR(B[1], 0.0, 1e-12);
}

// ===========================================================================
// Test: swap() does not crash
// ===========================================================================
TEST_F(Astra1DDynamicTest, SwapNoCrash) {
    std::string fname =
            writeConstantAstra1DFieldmap(tmpFile("swap.map"), 0.0, 0.10, 5, 1.0, 8, 100.0, true);

    auto* fm = dynamic_cast<Astra1DDynamic*>(Fieldmap::getFieldmap(fname));
    ASSERT_NE(fm, nullptr);

    EXPECT_NO_THROW(fm->swap());
}

// ===========================================================================
// Test: getInfo does not crash
// ===========================================================================
TEST_F(Astra1DDynamicTest, GetInfoNoCrash) {
    std::string fname =
            writeConstantAstra1DFieldmap(tmpFile("info.map"), 0.0, 0.10, 5, 1.0, 8, 100.0, true);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);

    Inform msg("test");
    EXPECT_NO_THROW(fm->getInfo(&msg));
}

// ===========================================================================
// Test: Missing file
// ===========================================================================
TEST_F(Astra1DDynamicTest, MissingFile) {
    std::string fname = tmpFile("nonexistent.map");
    EXPECT_THROW(Fieldmap::getFieldmap(fname), GeneralOpalException);
}

// ===========================================================================
// Test: Dictionary caching
// ===========================================================================
TEST_F(Astra1DDynamicTest, DictionaryCaching) {
    std::string fname =
            writeConstantAstra1DFieldmap(tmpFile("cache.map"), 0.0, 0.10, 5, 1.0, 8, 100.0, true);

    Fieldmap* fm1 = Fieldmap::getFieldmap(fname);
    Fieldmap* fm2 = Fieldmap::getFieldmap(fname);

    EXPECT_EQ(fm1, fm2);
}

// ===========================================================================
// Test: readMap / freeMap cycle
// ===========================================================================
TEST_F(Astra1DDynamicTest, ReadFreeCycle) {
    std::string fname =
            writeConstantAstra1DFieldmap(tmpFile("cycle.map"), 0.0, 0.10, 9, 1.0, 8, 100.0, true);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);

    fm->readMap();

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    fm->getFieldstrength(R, E, B);
    EXPECT_NE(E[2], 0.0);

    fm->freeMap();
    fm->readMap();

    E = {0.0, 0.0, 0.0};
    B = {0.0, 0.0, 0.0};

    fm->getFieldstrength(R, E, B);
    EXPECT_NE(E[2], 0.0);
}

// ===========================================================================
// Test: applyTravelingWave entry region
// ===========================================================================
TEST_F(Astra1DDynamicTest, ComputeTravelingWaveFieldEntryRegion) {
    Kokkos::View<double*, Kokkos::HostSpace> coefs("coefs", 3);
    coefs(0) = 1.0;
    coefs(1) = 0.0;
    coefs(2) = 0.0;

    Vector_t<double, 3> R = {0.0, 0.0, -0.01};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    Astra1DDynamic::computeTravelingWaveField(
            R, E, B, coefs,
            0.0,       // zbegin
            0.10,      // zend
            0.20,      // length
            1.0,       // xlrep
            2,         // accuracy
            2.0, 0.0,  // entry
            0.0, 0.0,  // core1
            0.0, 0.0,  // core2
            0.0, 0.0,  // exit
            0.05,      // startCoreField
            0.20,      // startExitField
            0.0,       // mappedStartExitField
            0.10,      // periodLength
            0.05,      // cellLength
            0.30       // elementLength
    );

    EXPECT_NEAR(E[2], 2.0, 1e-12);
}

// ===========================================================================
// Test: computeTravelingWaveField core region accumulates two contributions
// ===========================================================================
TEST_F(Astra1DDynamicTest, ComputeTravelingWaveFieldCoreRegion) {
    Kokkos::View<double*, Kokkos::HostSpace> coefs("coefs", 3);
    coefs(0) = 1.0;
    coefs(1) = 0.0;
    coefs(2) = 0.0;

    Vector_t<double, 3> R = {0.0, 0.0, 0.03};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    // Consistent geometry for a fieldmap spanning z in [0.0, 0.10)
    const double periodLength    = 0.05;
    const double startCoreField  = 0.025;
    const double cellLength      = 0.05;
    const double startExitField  = 0.125;
    const double mappedStartExit = 0.05;
    const double elementLength   = 0.175;

    Astra1DDynamic::computeTravelingWaveField(
            R, E, B, coefs,
            0.0,       // zbegin
            0.10,      // zend
            0.20,      // length
            1.0,       // xlrep
            2,         // accuracy
            0.0, 0.0,  // entry
            1.5, 0.0,  // core1
            2.5, 0.0,  // core2
            0.0, 0.0,  // exit
            startCoreField, startExitField, mappedStartExit, periodLength, cellLength,
            elementLength);

    // Uniform constant coefficient => each core contribution gives Ez = 1
    EXPECT_NEAR(E[2], 4.0, 1e-12);
}

// ===========================================================================
// Test: computeTravelingWaveField exit region
// ===========================================================================
TEST_F(Astra1DDynamicTest, ComputeTravelingWaveFieldExitRegion) {
    Kokkos::View<double*, Kokkos::HostSpace> coefs("coefs", 3);
    coefs(0) = 1.0;
    coefs(1) = 0.0;
    coefs(2) = 0.0;

    Vector_t<double, 3> R = {0.0, 0.0, 0.11};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    // Same consistent TW geometry as above
    const double periodLength    = 0.05;
    const double startCoreField  = 0.025;
    const double cellLength      = 0.05;
    const double startExitField  = 0.125;
    const double mappedStartExit = 0.05;
    const double elementLength   = 0.175;

    Astra1DDynamic::computeTravelingWaveField(
            R, E, B, coefs,
            0.0,       // zbegin
            0.10,      // zend
            0.20,      // length
            1.0,       // xlrep
            2,         // accuracy
            0.0, 0.0,  // entry
            0.0, 0.0,  // core1
            0.0, 0.0,  // core2
            3.0, 0.0,  // exit
            startCoreField, startExitField, mappedStartExit, periodLength, cellLength,
            elementLength);

    EXPECT_NEAR(E[2], 3.0, 1e-12);
}

// ===========================================================================
// Test: applyTravelingWave does nothing outside TW longitudinal range
// ===========================================================================
TEST_F(Astra1DDynamicTest, ComputeTravelingWaveFieldOutsideRange) {
    Kokkos::View<double*, Kokkos::HostSpace> coefs("coefs", 3);
    coefs(0) = 1.0;
    coefs(1) = 0.0;
    coefs(2) = 0.0;

    Vector_t<double, 3> R = {0.0, 0.0, -0.06};
    Vector_t<double, 3> E = {1.0, 2.0, 3.0};
    Vector_t<double, 3> B = {4.0, 5.0, 6.0};

    Astra1DDynamic::computeTravelingWaveField(
            R, E, B, coefs, 0.0, 0.10, 0.20, 1.0, 2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.05,
            0.20, 0.0, 0.10, 0.05, 0.30);

    EXPECT_NEAR(E[0], 1.0, 1e-12);
    EXPECT_NEAR(E[1], 2.0, 1e-12);
    EXPECT_NEAR(E[2], 3.0, 1e-12);
    EXPECT_NEAR(B[0], 4.0, 1e-12);
    EXPECT_NEAR(B[1], 5.0, 1e-12);
    EXPECT_NEAR(B[2], 6.0, 1e-12);
}
