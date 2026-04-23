/**
 * \file TestFM2DDynamic.cpp
 * \brief Unit tests for FM2DDynamic fieldmap class
 *
 * Tests cover:
 *
 * Parsing & metadata:
 * - File parsing (XZ and ZX orientations)
 * - Frequency parsing and conversion (MHz → angular frequency)
 * - Field dimension retrieval (z-bounds)
 *
 * Field evaluation:
 * - Point field evaluation via getFieldstrength()
 * - Uniform field correctness (on-axis behaviour)
 * - Bilinear interpolation for spatially varying fields
 * - Exact values at grid points (GridPointAccuracy)
 *
 * Geometry & bounds:
 * - Behaviour outside z-range (no modification, outside flag)
 * - Behaviour outside radial range
 * - isInside() helper (including boundary semantics)
 *
 * Physical behaviour:
 * - Normalization and scaling:
 *   * Electric field: MV/m → V/m
 *   * Magnetic field: Bt → μ0 scaling
 * - Projection of cylindrical components (Er, Bt) into Cartesian (x,y)
 * - Consistency between XZ and ZX orientations
 *
 * Low-level kernel:
 * - Static computeField() correctness (host views, interpolation, projection)
 *
 * API behaviour:
 * - swap() toggle (orientation handling)
 * - getFieldDerivative() throws (not implemented)
 * - getFieldDimensions(6-arg) throws
 * - getFrequency() correctness
 * - getInfo() does not crash
 *
 * Robustness:
 * - Missing / invalid file handling
 * - Fieldmap dictionary caching (singleton behaviour)
 * - readMap() / freeMap() lifecycle correctness
 *
 * Semantics:
 * - Field accumulation: E and B are accumulated, not overwritten
 */

#include <gtest/gtest.h>

#include "Fields/Fieldmap.h"
#include "Fields/FM2DDynamic.h"
#include "Physics/Units.h"
#include "Utilities/GeneralOpalException.h"
#include "Ippl.h"

#include <cmath>
#include <fstream>
#include <filesystem>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Helper: write dynamic XZ fieldmap
// Format:
// 2DDynamic XZ
// zbegin zend nz
// freq (MHz)
// rbegin rend nr
// data: Ez Er dummy Bt
// ---------------------------------------------------------------------------
std::string writeXZFieldmap(const std::string& path,
                            double zbegin_cm, double zend_cm, int nz,
                            double rbegin_cm, double rend_cm, int nr,
                            double freq_MHz = 100.0,
                            double uniformEz = 1.0,
                            double uniformEr = 0.0,
                            double uniformBt = 0.0)
{
    std::ofstream f(path);
    f << "2DDynamic XZ\n";
    f << zbegin_cm << " " << zend_cm << " " << nz << "\n";
    f << freq_MHz << "\n";
    f << rbegin_cm << " " << rend_cm << " " << nr << "\n";

    for (int jr = 0; jr <= nr; ++jr) {
        for (int iz = 0; iz <= nz; ++iz) {
            f << uniformEz << " " << uniformEr << " 0.0 " << uniformBt << "\n";
        }
    }
    return path;
}

// ---------------------------------------------------------------------------
// Helper: write dynamic fieldmap in ZX orientation
// Format:
// 2DDynamic ZX
// rbegin rend nr
// freq
// zbegin zend nz
// data: Er Ez Bt dummy
// ---------------------------------------------------------------------------
std::string writeZXFieldmap(const std::string& path,
                            double zbegin_cm, double zend_cm, int nz,
                            double rbegin_cm, double rend_cm, int nr,
                            double freq_MHz = 100.0,
                            double uniformEz = 1.0,
                            double uniformEr = 0.0,
                            double uniformBt = 0.0)
{
    std::ofstream f(path);
    f << "2DDynamic ZX\n";
    f << rbegin_cm << " " << rend_cm << " " << nr << "\n";
    f << freq_MHz << "\n";
    f << zbegin_cm << " " << zend_cm << " " << nz << "\n";

    for (int i = 0; i <= nz; ++i) {
        for (int j = 0; j <= nr; ++j) {
            f << uniformEr << " " << uniformEz << " " << uniformBt << " 0.0\n";
        }
    }
    return path;
}

// ---------------------------------------------------------------------------
// Helper: Write a fieldmap with spatially varying field for interpolation tests.
// Uses XZ orientation. Field: Ez = z [MV/m], Er = r [MV/m], Bt = r [T], 
// where r and z are in meters.
// ---------------------------------------------------------------------------
std::string writeVaryingXZFieldmap(const std::string& path, 
                                   double zbegin_cm, double zend_cm, int nz,
                                   double rbegin_cm, double rend_cm, int nr,
                                   double freq_MHz)
{
    const double hz_cm = (zend_cm - zbegin_cm) / nz;
    const double hr_cm = (rend_cm - rbegin_cm) / nr;

    std::ofstream f(path);
    f << "2DDynamic XZ\n";
    f << zbegin_cm << " " << zend_cm << " " << nz << "\n";
    f << freq_MHz << "\n";
    f << rbegin_cm << " " << rend_cm << " " << nr << "\n";

    for (int jr = 0; jr <= nr; ++jr) {
        double r_m = (rbegin_cm + jr * hr_cm) * Units::cm2m;
        for (int iz = 0; iz <= nz; ++iz) {
            double z_m = (zbegin_cm + iz * hz_cm) * Units::cm2m;
            f << z_m << " " << r_m << " 0.0 " << r_m << "\n";
        }
    }
    return path;
}

} // namespace

// ===========================================================================
// Fixture
// ===========================================================================
class FM2DDynamicTest : public ::testing::Test {
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
        tmpDir_ = std::filesystem::temp_directory_path() / "opalx_fm2d_test";
        std::filesystem::create_directories(tmpDir_);
    }

    void TearDown() override {
        Fieldmap::clearDictionary();
        std::filesystem::remove_all(tmpDir_);
    }

    std::string tmpFile(const std::string& name) const {
        return (tmpDir_ / name).string();
    }

    std::filesystem::path tmpDir_;
};

// ===========================================================================
// Test: Parse XZ orientation and verify field dimensions
// ===========================================================================
TEST_F(FM2DDynamicTest, ParseXZOrientation) {
    const double zb = 0.0, ze = 10.0;  // cm
    const double rb = 0.0, re = 5.0;   // cm
    const double freq = 100.0;         // MHz
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("xz.map"), 
                                        zb, ze, nz, 
                                        rb, re, nr, 
                                        freq);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);
    EXPECT_EQ(fm->getType(), T2DDynamic);

    Fieldmap::readMap(fname);

    double zBegin, zEnd;
    fm->getFieldDimensions(zBegin, zEnd);

    EXPECT_NEAR(zBegin, zb * Units::cm2m, 1e-12);
    EXPECT_NEAR(zEnd,   ze * Units::cm2m, 1e-12);
}

// ===========================================================================
// Test: Parse ZX orientation and verify field dimensions
// ===========================================================================
TEST_F(FM2DDynamicTest, ParseZXOrientation) {
    const double zb = -5.0, ze = 15.0;  // cm
    const double rb = 0.0,  re = 3.0;   // cm
    const double freq = 100.0;          // MHz
    const int nz = 3, nr = 2;

    std::string fname = writeZXFieldmap(tmpFile("zx.map"), 
                                        zb, ze, nz, 
                                        rb, re, nr, 
                                        freq);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    double zBegin, zEnd;
    fm->getFieldDimensions(zBegin, zEnd);

    EXPECT_NEAR(zBegin, zb * Units::cm2m, 1e-12);
    EXPECT_NEAR(zEnd,   ze * Units::cm2m, 1e-12);
}

// ===========================================================================
// Test: Frequency parsing and conversion
// ===========================================================================
TEST_F(FM2DDynamicTest, FrequencyParsing) {
    const double zb = 0.0, ze = 10.0;  // cm
    const double rb = 0.0, re = 5.0;   // cm
    const double freq = 100.0;         // MHz
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(
        tmpFile("freq.map"), zb, ze, nz, rb, re, nr, freq);

    auto* fm = dynamic_cast<FM2DDynamic*>(Fieldmap::getFieldmap(fname));
    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    // Constructor converts: MHz → Hz → angular frequency
    double expected = freq * Physics::two_pi * Units::MHz2Hz; // Should be angular frequency in rad/s

    EXPECT_NEAR(fm->getFrequency(), expected, 1e-6);
}

// ===========================================================================
// Test: Uniform field – getFieldstrength returns correct value at grid centre
// ===========================================================================
TEST_F(FM2DDynamicTest, UniformFieldStrength) {
    const double zb = 0.0, ze = 20.0;
    const double rb = 0.0, re = 5.0;
    const double freq = 100.0;
    const int nz = 4, nr = 2;

    const double Ez_val = 1.0;
    const double Er_val = 0.0;
    const double Bt_val = 0.0;

    std::string fname = writeXZFieldmap(tmpFile("uni.map"), 
                                        zb, ze, nz, 
                                        rb, re, nr, 
                                        freq, 
                                        Ez_val, Er_val, Bt_val);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    // Query at centre of field (on-axis)
    Vector_t<double, 3> R = {0.0, 0.0, 0.05};  // z=0.05m
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool outside = fm->getFieldstrength(R, E, B);

    EXPECT_FALSE(outside);

    // On-axis: Er = 0 → no transverse E
    EXPECT_NEAR(E[0], 0.0, 1e-10);
    EXPECT_NEAR(E[1], 0.0, 1e-10);

    // Longitudinal electric field
    EXPECT_NEAR(E[2], Ez_val * 1e6, 1e-10);  // MV/m → V/m scaling

    // No magnetic field
    EXPECT_NEAR(B[0], 0.0, 1e-10);
    EXPECT_NEAR(B[1], 0.0, 1e-10);
    EXPECT_NEAR(B[2], 0.0, 1e-10);
}

// ===========================================================================
// Test: Field outside z range returns without modifying E/B
// ===========================================================================
TEST_F(FM2DDynamicTest, OutsideZRange) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const double freq = 100.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("oz.map"), 
                                        zb, ze, nz, 
                                        rb, re, nr, 
                                        freq);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    // Before z start
    {
        Vector_t<double, 3> R = {0.0, 0.0, -0.01};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        bool outside = fm->getFieldstrength(R, E, B);

        EXPECT_TRUE(outside);

        // E and B should remain unchanged
        EXPECT_NEAR(E[0], 0.0, 1e-15);
        EXPECT_NEAR(E[1], 0.0, 1e-15);
        EXPECT_NEAR(E[2], 0.0, 1e-15);

        EXPECT_NEAR(B[0], 0.0, 1e-15);
        EXPECT_NEAR(B[1], 0.0, 1e-15);
        EXPECT_NEAR(B[2], 0.0, 1e-15);
    }

    // After z end
    {
        Vector_t<double, 3> R = {0.0, 0.0, 0.20};  // 20cm > 10cm
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        bool outside = fm->getFieldstrength(R, E, B);

        EXPECT_TRUE(outside);

        EXPECT_NEAR(E[0], 0.0, 1e-15);
        EXPECT_NEAR(E[1], 0.0, 1e-15);
        EXPECT_NEAR(E[2], 0.0, 1e-15);

        EXPECT_NEAR(B[0], 0.0, 1e-15);
        EXPECT_NEAR(B[1], 0.0, 1e-15);
        EXPECT_NEAR(B[2], 0.0, 1e-15);
    }
}

// ===========================================================================
// Test: Field outside radial range returns true (outside flag)
// ===========================================================================
TEST_F(FM2DDynamicTest, OutsideRadialRange) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 2.0;
    const double freq = 100.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("or.map"), 
                                        zb, ze, nz, 
                                        rb, re, nr, 
                                        freq);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    // r = 0.05 m = 5 cm > re = 2 cm
    Vector_t<double, 3> R = {0.05, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool outside = fm->getFieldstrength(R, E, B);

    EXPECT_TRUE(outside);

    EXPECT_NEAR(E[0], 0.0, 1e-15);
    EXPECT_NEAR(E[1], 0.0, 1e-15);
    EXPECT_NEAR(E[2], 0.0, 1e-15);

    EXPECT_NEAR(B[0], 0.0, 1e-15);
    EXPECT_NEAR(B[1], 0.0, 1e-15);
    EXPECT_NEAR(B[2], 0.0, 1e-15);
}


// ===========================================================================
// Test: isInside() helper
// ===========================================================================
TEST_F(FM2DDynamicTest, IsInside) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const double freq = 100.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("ins.map"),
                                        zb, ze, nz,
                                        rb, re, nr, 
                                        freq);

    auto* fm = dynamic_cast<FM2DDynamic*>(Fieldmap::getFieldmap(fname));
    ASSERT_NE(fm, nullptr);

    // Inside
    EXPECT_TRUE(fm->isInside({0.0, 0.0, 0.05}));

    // On zbegin boundary
    EXPECT_TRUE(fm->isInside({0.0, 0.0, 0.0}));

    // Just before zend (zend = 0.10 m)
    EXPECT_TRUE(fm->isInside({0.0, 0.0, 0.099}));

    // At zend – should be outside (< zend, not <=)
    EXPECT_FALSE(fm->isInside({0.0, 0.0, 0.10}));

    // Before zbegin
    EXPECT_FALSE(fm->isInside({0.0, 0.0, -0.01}));

    // Outside radial
    EXPECT_FALSE(fm->isInside({0.10, 0.0, 0.05}));
}

// ===========================================================================
// Test: swap() toggles the swap state
// ===========================================================================
TEST_F(FM2DDynamicTest, SwapToggle) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const double freq = 100.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("sw.map"),
                                       zb, ze, nz,
                                       rb, re, nr, 
                                       freq);

    auto* fm = dynamic_cast<FM2DDynamic*>(Fieldmap::getFieldmap(fname));
    ASSERT_NE(fm, nullptr);

    // swap_m starts false for XZ
    fm->swap();   // now true
    fm->swap();   // back to false

    // No crash – just verify the toggle doesn't explode
}

// ===========================================================================
// Test: Normalization flag
// ===========================================================================
TEST_F(FM2DDynamicTest, Normalization) {
    // normalize_m defaults to true.
    // With uniform Ez = 3.0 (MV/m), Ezmax = 3.0,
    // so after normalization Ez -> 1.0 * 1e6 V/m

    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const double freq = 100.0;
    const int nz = 4, nr = 2;

    const double Ez_val = 3.0; // MV/m

    std::string fname = writeXZFieldmap(tmpFile("norm.map"),
                                        zb, ze, nz,
                                        rb, re, nr,
                                        freq,
                                        /*Ez=*/Ez_val,
                                        /*Er=*/0.0,
                                        /*Bt=*/0.0);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    fm->getFieldstrength(R, E, B);

    // After normalization:
    // Ez_normalized = Ez / Ezmax = 1.0
    // Then scaled: 1.0 * 1e6 V/m
    EXPECT_NEAR(E[2], 1e6, 1e-6);
}

// ===========================================================================
// Test: Bilinear interpolation with spatially varying field
// ===========================================================================
TEST_F(FM2DDynamicTest, BilinearInterpolation) {
    // Field before normalization:
    // Ez(r,z) = z [MV/m], Er(r,z) = r [MV/m], Bt(r,z) = r
    //
    // After normalization + scaling:
    // Ez = (z / Ezmax) * 1e6
    // Er = (r / Ezmax) * 1e6
    // Bt = (r / Ezmax) * mu0

    const double zb = 0.0, ze = 10.0;   // cm
    const double rb = 0.0, re = 5.0;    // cm
    const double freq = 100.0;
    const int nz = 4, nr = 4;

    const double zend_m = ze * Units::cm2m;  // 0.10 m
    const double Ezmax  = zend_m;            // same logic as magnetostatic

    std::string fname = writeVaryingXZFieldmap(tmpFile("var.map"), 
                                               zb, ze, nz, 
                                               rb, re, nr, 
                                               freq);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    Fieldmap::readMap(fname);

    const double scaleE = 1e6 / Ezmax;

    // ----------------------------------------------------------------------
    // On-axis: z = 5cm = 0.05m
    // ----------------------------------------------------------------------
    {
        Vector_t<double, 3> R = {0.0, 0.0, 0.05};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        fm->getFieldstrength(R, E, B);

        EXPECT_NEAR(E[2], 0.05 * scaleE, 1e-6);
        EXPECT_NEAR(E[0], 0.0, 1e-10);
        EXPECT_NEAR(E[1], 0.0, 1e-10);
    }

    // ----------------------------------------------------------------------
    // Off-axis along x
    // ----------------------------------------------------------------------
    {
        double r_m = 0.025;  // 2.5 cm
        Vector_t<double, 3> R = {r_m, 0.0, 0.05};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        fm->getFieldstrength(R, E, B);

        EXPECT_NEAR(E[2], 0.05 * scaleE, 1e-6);

        // Er projected onto x
        EXPECT_NEAR(E[0], r_m * scaleE, 1e-6);
        EXPECT_NEAR(E[1], 0.0, 1e-10);

        // Bt produces transverse B
        EXPECT_NEAR(B[1], r_m * Physics::mu_0 / Ezmax, 1e-10);
    }

    // ----------------------------------------------------------------------
    // Off-axis along y
    // ----------------------------------------------------------------------
    {
        double r_m = 0.025;
        Vector_t<double, 3> R = {0.0, r_m, 0.05};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        fm->getFieldstrength(R, E, B);

        EXPECT_NEAR(E[2], 0.05 * scaleE, 1e-6);

        EXPECT_NEAR(E[0], 0.0, 1e-10);
        EXPECT_NEAR(E[1], r_m * scaleE, 1e-6);

        EXPECT_NEAR(B[0], -r_m * Physics::mu_0 / Ezmax, 1e-10);
    }

    // ----------------------------------------------------------------------
    // Non-grid z value
    // ----------------------------------------------------------------------
    {
        Vector_t<double, 3> R = {0.0, 0.0, 0.03};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        fm->getFieldstrength(R, E, B);

        EXPECT_NEAR(E[2], 0.03 * scaleE, 1e-6);
    }
}

// ===========================================================================
// Test: computeField static function directly with host views
// ===========================================================================
TEST_F(FM2DDynamicTest, ComputeFieldStatic) {
    // Build a small grid manually: 3 z-points, 2 r-points
    // Ez = 1.0 everywhere, Er = 0.5 everywhere, Bt = 0.25 everywhere

    const int npz = 3, npr = 2;
    const double hz = 0.05, hr = 0.025;
    const double zbegin = 0.0;

    Kokkos::View<double*, Kokkos::HostSpace> Ez("Ez", npz * npr);
    Kokkos::View<double*, Kokkos::HostSpace> Er("Er", npz * npr);
    Kokkos::View<double*, Kokkos::HostSpace> Bt("Bt", npz * npr);

    for (int j = 0; j < npr; ++j)
        for (int i = 0; i < npz; ++i) {
            Ez(j * npz + i) = 1.0;
            Er(j * npz + i) = 0.5;
            Bt(j * npz + i) = 0.25;
        }

    // ----------------------------------------------------------------------
    // On-axis: R = (0, 0, 0.025)
    // ----------------------------------------------------------------------
    {
        Vector_t<double, 3> R = {0.0, 0.0, 0.025};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        FM2DDynamic::computeField(R, E, B, Ez, Er, Bt,
                                  hr, hz, zbegin, npr, npz);

        EXPECT_NEAR(E[2], 1.0, 1e-10);

        // On-axis RR=0 → no transverse components
        EXPECT_NEAR(E[0], 0.0, 1e-10);
        EXPECT_NEAR(E[1], 0.0, 1e-10);
        EXPECT_NEAR(B[0], 0.0, 1e-10);
        EXPECT_NEAR(B[1], 0.0, 1e-10);
    }

    // ----------------------------------------------------------------------
    // Off-axis: R = (0.01, 0, 0.025)
    // ----------------------------------------------------------------------
    {
        Vector_t<double, 3> R = {0.01, 0.0, 0.025};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        FM2DDynamic::computeField(R, E, B, Ez, Er, Bt,
                                  hr, hz, zbegin, npr, npz);

        EXPECT_NEAR(E[2], 1.0, 1e-10);

        // Er projected → x
        EXPECT_NEAR(E[0], 0.5, 1e-10);
        EXPECT_NEAR(E[1], 0.0, 1e-10);

        // Bt projected → y
        EXPECT_NEAR(B[1], 0.25, 1e-10);
        EXPECT_NEAR(B[0], 0.0, 1e-10);
    }

    // ----------------------------------------------------------------------
    // Outside z range (before): no modification
    // ----------------------------------------------------------------------
    {
        Vector_t<double, 3> R = {0.0, 0.0, -0.01};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        FM2DDynamic::computeField(R, E, B, Ez, Er, Bt,
                                  hr, hz, zbegin, npr, npz);

        EXPECT_NEAR(E[2], 0.0, 1e-15);
        EXPECT_NEAR(B[0], 0.0, 1e-15);
    }

    // ----------------------------------------------------------------------
    // Outside r range: no modification
    // ----------------------------------------------------------------------
    {
        Vector_t<double, 3> R = {0.05, 0.0, 0.025};  // outside radial grid
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        FM2DDynamic::computeField(R, E, B, Ez, Er, Bt,
                                  hr, hz, zbegin, npr, npz);

        EXPECT_NEAR(E[0], 0.0, 1e-15);
        EXPECT_NEAR(E[1], 0.0, 1e-15);
        EXPECT_NEAR(E[2], 0.0, 1e-15);
        EXPECT_NEAR(B[0], 0.0, 1e-15);
        EXPECT_NEAR(B[1], 0.0, 1e-15);
    }
}

// ===========================================================================
// Test: getFieldDerivative throws
// ===========================================================================
TEST_F(FM2DDynamicTest, GetFieldDerivativeThrows) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const double freq = 100.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("deriv.map"),
                                        zb, ze, nz,
                                        rb, re, nr,
                                        freq);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    EXPECT_THROW(fm->getFieldDerivative(R, E, B, DX),
                 GeneralOpalException);
}

// ===========================================================================
// Test: getFieldDimensions 6-arg throws
// ===========================================================================
TEST_F(FM2DDynamicTest, GetFieldDimensions6ArgThrows) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const double freq = 100.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("dim6.map"),
                                        zb, ze, nz,
                                        rb, re, nr, 
                                        freq);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);

    double a, b, c, d, e, f;

    EXPECT_THROW(fm->getFieldDimensions(a, b, c, d, e, f),
                 GeneralOpalException);
}

// ===========================================================================
// Test: getInfo does not crash
// ===========================================================================
TEST_F(FM2DDynamicTest, GetInfoNoCrash) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const double freq = 100.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("info.map"), 
                                        zb, ze, nz,
                                        rb, re, nr,
                                        freq);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);

    Inform msg("test");
    EXPECT_NO_THROW(fm->getInfo(&msg));
}

// ===========================================================================
// Test: Missing file
// ===========================================================================
TEST_F(FM2DDynamicTest, MissingFile) {
    std::string fname = tmpFile("nonexistent.map");
    EXPECT_THROW(Fieldmap::getFieldmap(fname),
                 GeneralOpalException);
}

// ===========================================================================
// Test: Fieldmap dictionary caching – second getFieldmap returns same pointer
// ===========================================================================
TEST_F(FM2DDynamicTest, DictionaryCaching) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const double freq = 100.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("cache.map"),
                                        zb, ze, nz,
                                        rb, re, nr,
                                        freq);

    Fieldmap* fm1 = Fieldmap::getFieldmap(fname);
    Fieldmap* fm2 = Fieldmap::getFieldmap(fname);

    EXPECT_EQ(fm1, fm2);
}

// ===========================================================================
// Test: readMap / freeMap cycle – reading after free should work
// ===========================================================================
TEST_F(FM2DDynamicTest, ReadFreeCycle) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const double freq = 100.0;
    const int nz = 4, nr = 2;


    std::string fname = writeXZFieldmap(tmpFile("cycle.map"),
                                        zb, ze, nz,
                                        rb, re, nr,
                                        freq);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    fm->getFieldstrength(R, E, B);

    // sanity check: field exists (not zero / not crashing)
    EXPECT_NE(E[2], 0.0);

    fm->freeMap();

    fm->readMap();

    E = {0.0, 0.0, 0.0};
    B = {0.0, 0.0, 0.0};

    fm->getFieldstrength(R, E, B);

    EXPECT_NE(E[2], 0.0);
}

// ===========================================================================
// Test: On-axis at grid points matches exact data (dynamic fieldmap)
// ===========================================================================
TEST_F(FM2DDynamicTest, GridPointAccuracy) {
    // Field before normalization: Ez(r,z) = z, Er = r
    // With normalization: Ezmax = zend_m, so normalised Ez = z / zend_m

    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const double freq = 100.0;
    const int nz = 10, nr = 5;

    const double zb_m = zb * Units::cm2m;
    const double ze_m = ze * Units::cm2m;
    //const double zend_m = ze_m;

    std::string fname = writeVaryingXZFieldmap(tmpFile("grid.map"),
                                               zb, ze, nz,
                                               rb, re, nr,
                                               freq);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    const double hz_m = (ze_m - zb_m) / nz; 

    double Ezmax = 1.0;
    double scaleE = 1e6 / Ezmax;

    for (int iz = 0; iz < nz; ++iz) {
        double z = zb_m + iz * hz_m;

        Vector_t<double, 3> R = {0.0, 0.0, z};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        fm->getFieldstrength(R, E, B);

        // Dynamic: electric field carries the longitudinal component
        EXPECT_NEAR(E[2], (z / ze_m) * scaleE, 1e-9)
            << "Mismatch at iz=" << iz << " z=" << z;

        // No transverse components on-axis
        EXPECT_NEAR(E[0], 0.0, 1e-12);
        EXPECT_NEAR(E[1], 0.0, 1e-12);
        EXPECT_NEAR(B[0], 0.0, 1e-12);
        EXPECT_NEAR(B[1], 0.0, 1e-12);
    }
}

// ===========================================================================
// Test: ZX orientation gives same results as XZ for identical physical field
// ===========================================================================
TEST_F(FM2DDynamicTest, XZvsZXConsistency) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const double freq = 100.0;
    const int nz = 4, nr = 2;
    const double Ez_val = 2.0;

    std::string fnameXZ = writeXZFieldmap(tmpFile("xz2.map"), 
                                          zb, ze, nz, 
                                          rb, re, nr,
                                          freq,
                                          Ez_val);
    std::string fnameZX = writeZXFieldmap(tmpFile("zx2.map"), 
                                          zb, ze, nz, 
                                          rb, re, nr,
                                          freq,
                                          Ez_val);

    Fieldmap* fmXZ = Fieldmap::getFieldmap(fnameXZ);
    Fieldmap* fmZX = Fieldmap::getFieldmap(fnameZX);
    Fieldmap::readMap(fnameXZ);
    Fieldmap::readMap(fnameZX);

    ASSERT_NE(fmXZ, nullptr);
    ASSERT_NE(fmZX, nullptr);

    Vector_t<double, 3> R    = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E_xz = {0.0, 0.0, 0.0};
    Vector_t<double, 3> E_zx = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B_xz = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B_zx = {0.0, 0.0, 0.0};

    fmXZ->getFieldstrength(R, E_xz, B_xz);
    fmZX->getFieldstrength(R, E_zx, B_zx);

    EXPECT_NEAR(E_xz[0], E_zx[0], 1e-10);
    EXPECT_NEAR(E_xz[1], E_zx[1], 1e-10);
    EXPECT_NEAR(E_xz[2], E_zx[2], 1e-10);

    EXPECT_NEAR(B_xz[0], B_zx[0], 1e-10);
    EXPECT_NEAR(B_xz[1], B_zx[1], 1e-10);
    EXPECT_NEAR(B_xz[2], B_zx[2], 1e-10);
}

// ===========================================================================
// Test: Er and Bt project correctly into Cartesian components (Dynamic)
// ===========================================================================
TEST_F(FM2DDynamicTest, FieldProjection) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const double freq = 100.0;
    const int nz = 4, nr = 2;

    // Define uniform fields
    const double Ez_val = 1.0;
    const double Er_val = 1.0;
    const double Bt_val = 1.0;

    std::string fname = writeXZFieldmap(tmpFile("proj_dyn.map"), 
                                        zb, ze, nz, 
                                        rb, re, nr, 
                                        freq,
                                        Ez_val, Er_val, Bt_val);

    Fieldmap* fm = dynamic_cast<FM2DDynamic*>(Fieldmap::getFieldmap(fname));
    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    // 45 degrees in transverse plane
    double r = 0.02;
    double x = r / std::sqrt(2.0);
    double y = r / std::sqrt(2.0);

    Vector_t<double, 3> R = {x, y, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    fm->getFieldstrength(R, E, B);

    double Ezmax = Ez_val; // uniform field → max = value
    double scaleE = 1e6 / Ezmax;
    double scaleB = Physics::mu_0 / Ezmax;

    double expected_Er = (Er_val / std::sqrt(2.0)) * scaleE;
    double expected_Ez = Ez_val * scaleE;
    double expected_Bt = (Bt_val / std::sqrt(2.0)) * scaleB;

    // --- Electric field (from Er) ---
    EXPECT_NEAR(E[0], expected_Er, 1e-10);
    EXPECT_NEAR(E[1], expected_Er, 1e-10);

    // --- Longitudinal electric ---
    EXPECT_NEAR(E[2], expected_Ez, 1e-10);

    // --- Magnetic field (from Bt) ---
    // Note the rotation!
    EXPECT_NEAR(B[0], -expected_Bt, 1e-10);  // -Bt * y/R
    EXPECT_NEAR(B[1],  expected_Bt, 1e-10);  // +Bt * x/R

    // No longitudinal magnetic component
    EXPECT_NEAR(B[2], 0.0, 1e-10);
}

// ===========================================================================
// Test: Accumulation semantics – E and B are accumulated, not overwritten
// ===========================================================================
TEST_F(FM2DDynamicTest, FieldAccumulation) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const double freq = 100.0;
    const int nz = 4, nr = 2;

    // Uniform fields
    const double Ez_val = 1.0;
    const double Er_val = 0.0;
    const double Bt_val = 0.0;

    std::string fname = writeXZFieldmap(tmpFile("accum_dyn.map"),
                                        zb, ze, nz,
                                        rb, re, nr,
                                        freq,
                                        Ez_val, Er_val, Bt_val);

    Fieldmap* fm = dynamic_cast<FM2DDynamic*>(Fieldmap::getFieldmap(fname));
    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 5.0}; // Pre-existing Efields
    Vector_t<double, 3> B = {0.0, 0.0, 3.0}; // Pre-existing Bfields

    fm->getFieldstrength(R, E, B);

    // --- Expected scaling ---
    double Ezmax = Ez_val;
    double scaleE = 1e6 / Ezmax;

    // Only Ez contributes here
    double expected_added_Ez = Ez_val * scaleE;

    // --- Check accumulation ---
    EXPECT_NEAR(E[0], 0.0, 1e-10);
    EXPECT_NEAR(E[1], 0.0, 1e-10);
    EXPECT_NEAR(E[2], 5.0 + expected_added_Ez, 1e-10);

    // No magnetic contribution (Bt = 0)
    EXPECT_NEAR(B[0], 0.0, 1e-10);
    EXPECT_NEAR(B[1], 0.0, 1e-10);
    EXPECT_NEAR(B[2], 3.0, 1e-10);
}
