/**
 * \file TestFM2DMagnetoStatic.cpp
 * \brief Unit tests for FM2DMagnetoStatic fieldmap class
 *
 * Tests cover:
 * - File parsing (XZ and ZX orientations)
 * - Field dimensions retrieval
 * - Point field evaluation via getFieldstrength (bilinear interpolation)
 * - Boundary / out-of-range behaviour
 * - isInside() helper
 * - swap() toggle
 * - Normalization flag
 * - Static computeField() kernel
 * - getFieldDerivative / getFieldDimensions(6-arg) / getFrequency / setFrequency throw
 * - getInfo prints without crashing
 * - Invalid / missing file handling
 */

#include <gtest/gtest.h>

#include "Fields/Fieldmap.h"
#include "Fields/FM2DMagnetoStatic.h"
#include "Physics/Units.h"
#include "Utilities/GeneralOpalException.h"
#include "Ippl.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>

namespace {

// ---------------------------------------------------------------------------
// Helper: write a minimal 2D magnetostatic fieldmap file in XZ orientation.
//
// Grid: nr+1 radial points from rbegin..rend (cm)
//        nz+1 longitudinal points from zbegin..zend (cm)
//
// NOTE: normalize_m defaults to true in the Fieldmap base class when the
// header has only two strings ("2DMagnetoStatic XZ").  The three-string
// form ("2DMagnetoStatic XZ TRUE/FALSE") is broken in the current parser
// (all three output references alias the same variable, clobbering the
// orientation).  Therefore this helper always writes a two-string header
// and normalization is always active.
//
// Default: uniform Bz = 1.0 T, Br = 0.0 T  (unless overridden)
// ---------------------------------------------------------------------------
std::string writeXZFieldmap(const std::string& path,
                            double zbegin_cm, double zend_cm, int nz,
                            double rbegin_cm, double rend_cm, int nr,
                            double uniformBz = 1.0,
                            double uniformBr = 0.0) {
    std::ofstream f(path);
    f << "2DMagnetoStatic XZ\n";
    f << zbegin_cm << " " << zend_cm << " " << nz << "\n";
    f << rbegin_cm << " " << rend_cm << " " << nr << "\n";
    // XZ: outer loop r (nr+1), inner loop z (nz+1), columns: Bz Br
    for (int jr = 0; jr <= nr; ++jr) {
        for (int iz = 0; iz <= nz; ++iz) {
            f << uniformBz << " " << uniformBr << "\n";
        }
    }
    f.close();
    return path;
}

// ---------------------------------------------------------------------------
// Helper: write a minimal 2D magnetostatic fieldmap file in ZX orientation.
// ---------------------------------------------------------------------------
std::string writeZXFieldmap(const std::string& path,
                            double zbegin_cm, double zend_cm, int nz,
                            double rbegin_cm, double rend_cm, int nr,
                            double uniformBz = 1.0,
                            double uniformBr = 0.0) {
    std::ofstream f(path);
    f << "2DMagnetoStatic ZX\n";
    f << rbegin_cm << " " << rend_cm << " " << nr << "\n";
    f << zbegin_cm << " " << zend_cm << " " << nz << "\n";
    // ZX: outer loop z (nz+1), inner loop r (nr+1), columns: Br Bz
    for (int iz = 0; iz <= nz; ++iz) {
        for (int jr = 0; jr <= nr; ++jr) {
            f << uniformBr << " " << uniformBz << "\n";
        }
    }
    f.close();
    return path;
}

// ---------------------------------------------------------------------------
// Helper: write a fieldmap with spatially varying field for interpolation tests.
// Uses XZ orientation.  Field: Bz(r,z) = z [T], Br(r,z) = r [T]
// where r and z are in metres.
// ---------------------------------------------------------------------------
std::string writeVaryingXZFieldmap(const std::string& path,
                                   double zbegin_cm, double zend_cm, int nz,
                                   double rbegin_cm, double rend_cm, int nr) {
    const double hz_cm = (zend_cm - zbegin_cm) / nz;
    const double hr_cm = (rend_cm - rbegin_cm) / nr;

    std::ofstream f(path);
    f << "2DMagnetoStatic XZ\n";
    f << zbegin_cm << " " << zend_cm << " " << nz << "\n";
    f << rbegin_cm << " " << rend_cm << " " << nr << "\n";
    for (int jr = 0; jr <= nr; ++jr) {
        double r_m = (rbegin_cm + jr * hr_cm) * Units::cm2m;
        for (int iz = 0; iz <= nz; ++iz) {
            double z_m = (zbegin_cm + iz * hz_cm) * Units::cm2m;
            // Bz = z, Br = r
            f << z_m << " " << r_m << "\n";
        }
    }
    f.close();
    return path;
}

}  // anonymous namespace

// ===========================================================================
// Test fixture
// ===========================================================================
class FM2DMagnetoStaticTest : public ::testing::Test {
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

    /// Helper to build a full path inside the temp directory
    std::string tmpFile(const std::string& name) const {
        return (tmpDir_ / name).string();
    }

    std::filesystem::path tmpDir_;
};

// ===========================================================================
// Test: Parse XZ orientation and verify field dimensions
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, ParseXZOrientation) {
    const double zb = 0.0, ze = 10.0;  // cm
    const double rb = 0.0, re = 5.0;   // cm
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("xz.map"), zb, ze, nz, rb, re, nr);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);
    EXPECT_EQ(fm->getType(), T2DMagnetoStatic);

    Fieldmap::readMap(fname);

    double zBegin, zEnd;
    fm->getFieldDimensions(zBegin, zEnd);
    EXPECT_NEAR(zBegin, zb * Units::cm2m, 1e-12);
    EXPECT_NEAR(zEnd,   ze * Units::cm2m, 1e-12);
}

// ===========================================================================
// Test: Parse ZX orientation and verify field dimensions
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, ParseZXOrientation) {
    const double zb = -5.0, ze = 15.0;  // cm
    const double rb = 0.0,  re = 3.0;   // cm
    const int nz = 3, nr = 2;

    std::string fname = writeZXFieldmap(tmpFile("zx.map"), zb, ze, nz, rb, re, nr);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    double zBegin, zEnd;
    fm->getFieldDimensions(zBegin, zEnd);
    EXPECT_NEAR(zBegin, zb * Units::cm2m, 1e-12);
    EXPECT_NEAR(zEnd,   ze * Units::cm2m, 1e-12);
}

// ===========================================================================
// Test: Uniform field – getFieldstrength returns correct value at grid centre
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, UniformFieldStrength) {
    const double zb = 0.0, ze = 20.0;
    const double rb = 0.0, re = 5.0;
    const int nz = 4, nr = 2;
    // Use Bz_val = 1.0 so that normalization (Bz /= Bzmax) is a no-op.
    const double Bz_val = 1.0, Br_val = 0.0;

    std::string fname = writeXZFieldmap(tmpFile("uni.map"), zb, ze, nz, rb, re, nr,
                                        Bz_val, Br_val);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    Fieldmap::readMap(fname);

    // Query at centre of field (on-axis)
    Vector_t<double, 3> R = {0.0, 0.0, 0.05};  // z=5cm=0.05m, on axis
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool outside = fm->getFieldstrength(R, E, B);
    EXPECT_FALSE(outside);
    EXPECT_NEAR(B[2], Bz_val, 1e-10);
    EXPECT_NEAR(B[0], 0.0, 1e-10);
    EXPECT_NEAR(B[1], 0.0, 1e-10);
}

// ===========================================================================
// Test: Field outside z range returns without modifying B (early return 
//       from computeField when indexz out of bounds)
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, OutsideZRange) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("oz.map"), zb, ze, nz, rb, re, nr);
    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    Fieldmap::readMap(fname);

    // Before z start
    {
        Vector_t<double, 3> R = {0.0, 0.0, -0.01};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};
        fm->getFieldstrength(R, E, B);
        // B should remain 0 – computeField returns without writing
        EXPECT_NEAR(B[0], 0.0, 1e-15);
        EXPECT_NEAR(B[1], 0.0, 1e-15);
        EXPECT_NEAR(B[2], 0.0, 1e-15);
    }

    // After z end
    {
        Vector_t<double, 3> R = {0.0, 0.0, 0.20};  // 20cm > 10cm
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};
        fm->getFieldstrength(R, E, B);
        EXPECT_NEAR(B[0], 0.0, 1e-15);
        EXPECT_NEAR(B[1], 0.0, 1e-15);
        EXPECT_NEAR(B[2], 0.0, 1e-15);
    }
}

// ===========================================================================
// Test: Field outside radial range returns true (outside flag)
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, OutsideRadialRange) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 2.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("or.map"), zb, ze, nz, rb, re, nr);
    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    Fieldmap::readMap(fname);

    // r = 0.05 m = 5 cm > re = 2 cm
    Vector_t<double, 3> R = {0.05, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool outside = fm->getFieldstrength(R, E, B);
    EXPECT_TRUE(outside);
}

// ===========================================================================
// Test: isInside() helper
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, IsInside) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("ins.map"), zb, ze, nz, rb, re, nr);
    auto* fm = dynamic_cast<FM2DMagnetoStatic*>(Fieldmap::getFieldmap(fname));
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
TEST_F(FM2DMagnetoStaticTest, SwapToggle) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("sw.map"), zb, ze, nz, rb, re, nr);
    auto* fm = dynamic_cast<FM2DMagnetoStatic*>(Fieldmap::getFieldmap(fname));
    ASSERT_NE(fm, nullptr);
    // swap_m starts false for XZ
    fm->swap();   // now true
    fm->swap();   // back to false
    // No crash – just verify the toggle doesn't explode
}

// ===========================================================================
// Test: Normalization flag
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, Normalization) {
    // normalize_m defaults to true.  With a uniform field of Bz_val = 3.0,
    // Bzmax = 3.0, so after normalization every Bz value becomes 1.0.
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const int nz = 4, nr = 2;
    const double Bz_val = 3.0;

    std::string fname = writeXZFieldmap(tmpFile("norm.map"), zb, ze, nz, rb, re, nr,
                                        Bz_val);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};
    fm->getFieldstrength(R, E, B);

    // All Bz values are the same (uniform), so Bzmax = Bz_val.
    // Normalized: every Bz /= Bzmax => 1.0
    EXPECT_NEAR(B[2], 1.0, 1e-10);
}

// ===========================================================================
// Test: Bilinear interpolation with spatially varying field
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, BilinearInterpolation) {
    // Field before normalization: Bz(r,z) = z [T],  Br(r,z) = r [T]
    // Normalization (always active) divides all values by Bzmax.
    // On-axis z values range from zbegin_m=0 to zend_m, so Bzmax = zend_m.
    const double zb = 0.0, ze = 10.0;   // cm
    const double rb = 0.0, re = 5.0;    // cm
    const int nz = 4, nr = 4;
    const double zend_m = ze * Units::cm2m;  // 0.10 m
    const double Bzmax = zend_m;             // normalization factor

    std::string fname = writeVaryingXZFieldmap(tmpFile("var.map"), zb, ze, nz, rb, re, nr);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    Fieldmap::readMap(fname);

    // Test at a mid-grid point: z = 5cm = 0.05m, on-axis r = 0
    {
        Vector_t<double, 3> R = {0.0, 0.0, 0.05};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};
        fm->getFieldstrength(R, E, B);
        EXPECT_NEAR(B[2], 0.05 / Bzmax, 1e-10);
        EXPECT_NEAR(B[0], 0.0, 1e-10);
        EXPECT_NEAR(B[1], 0.0, 1e-10);
    }

    // Test off-axis: r = 2.5cm = 0.025m along x, z = 5cm = 0.05m
    {
        double r_m = 0.025;
        Vector_t<double, 3> R = {r_m, 0.0, 0.05};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};
        fm->getFieldstrength(R, E, B);
        EXPECT_NEAR(B[2], 0.05 / Bzmax, 1e-10);
        // Br(r=0.025) = r/Bzmax, projected onto x: Br*x/R = r_m/Bzmax
        EXPECT_NEAR(B[0], r_m / Bzmax, 1e-10);
        EXPECT_NEAR(B[1], 0.0, 1e-10);
    }

    // Test off-axis along y
    {
        double r_m = 0.025;
        Vector_t<double, 3> R = {0.0, r_m, 0.05};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};
        fm->getFieldstrength(R, E, B);
        EXPECT_NEAR(B[2], 0.05 / Bzmax, 1e-10);
        EXPECT_NEAR(B[0], 0.0, 1e-10);
        EXPECT_NEAR(B[1], r_m / Bzmax, 1e-10);
    }

    // Test at a non-grid z value: z = 3cm = 0.03m
    {
        Vector_t<double, 3> R = {0.0, 0.0, 0.03};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};
        fm->getFieldstrength(R, E, B);
        EXPECT_NEAR(B[2], 0.03 / Bzmax, 1e-10);
    }
}

// ===========================================================================
// Test: computeField static function directly with host views
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, ComputeFieldStatic) {
    // Build a small grid manually: 3 z-points, 2 r-points
    // Bz = 1.0 everywhere, Br = 0.5 everywhere
    const int npz = 3, npr = 2;
    const double hz = 0.05, hr = 0.025;
    const double zbegin = 0.0;

    Kokkos::View<double*, Kokkos::HostSpace> Bz("Bz", npz * npr);
    Kokkos::View<double*, Kokkos::HostSpace> Br("Br", npz * npr);

    for (int j = 0; j < npr; ++j)
        for (int i = 0; i < npz; ++i) {
            Bz(j * npz + i) = 1.0;
            Br(j * npz + i) = 0.5;
        }

    // On-axis evaluation: R = (0, 0, 0.025)
    {
        Vector_t<double, 3> R = {0.0, 0.0, 0.025};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};
        FM2DMagnetoStatic::computeField(R, B, Bz, Br, hr, hz, zbegin, npr, npz);
        EXPECT_NEAR(B[2], 1.0, 1e-10);
        // On-axis RR=0, so Br contribution skipped
        EXPECT_NEAR(B[0], 0.0, 1e-10);
    }

    // Off-axis: R = (0.01, 0, 0.025)
    {
        Vector_t<double, 3> R = {0.01, 0.0, 0.025};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};
        FM2DMagnetoStatic::computeField(R, B, Bz, Br, hr, hz, zbegin, npr, npz);
        EXPECT_NEAR(B[2], 1.0, 1e-10);
        EXPECT_NEAR(B[0], 0.5, 1e-10);   // Br * x/R = 0.5 * 1.0
        EXPECT_NEAR(B[1], 0.0, 1e-10);
    }

    // Outside z range (before): should return false with B unchanged
    {
        Vector_t<double, 3> R = {0.0, 0.0, -0.01};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};
        FM2DMagnetoStatic::computeField(R, B, Bz, Br, hr, hz, zbegin, npr, npz);
        EXPECT_NEAR(B[2], 0.0, 1e-15);
    }

    // Outside r range: should return true
    {
        Vector_t<double, 3> R = {0.05, 0.0, 0.025};  // r = 0.05 > hr*(npr-1)=0.025
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};
        FM2DMagnetoStatic::computeField(R, B, Bz, Br, hr, hz, zbegin, npr, npz);
        EXPECT_NEAR(B[0], 0.0, 1e-15);
        EXPECT_NEAR(B[1], 0.0, 1e-15);
        EXPECT_NEAR(B[2], 0.0, 1e-15);
    }
}

// ===========================================================================
// Test: getFieldDerivative throws
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, GetFieldDerivativeThrows) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("deriv.map"), zb, ze, nz, rb, re, nr);
    Fieldmap* fm = Fieldmap::getFieldmap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    EXPECT_THROW(fm->getFieldDerivative(R, E, B, DX), GeneralOpalException);
}

// ===========================================================================
// Test: getFieldDimensions 6-arg returns the full cylindrical bounding box
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, GetFieldDimensions6ArgReturnsBoundingBox) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("dim6.map"), zb, ze, nz, rb, re, nr);
    Fieldmap* fm = Fieldmap::getFieldmap(fname);

    double xIni = 0.0, xFinal = 0.0, yIni = 0.0, yFinal = 0.0, zIni = 0.0, zFinal = 0.0;
    fm->getFieldDimensions(xIni, xFinal, yIni, yFinal, zIni, zFinal);

    EXPECT_NEAR(xIni, -re * Units::cm2m, 1e-12);
    EXPECT_NEAR(xFinal, re * Units::cm2m, 1e-12);
    EXPECT_NEAR(yIni, -re * Units::cm2m, 1e-12);
    EXPECT_NEAR(yFinal, re * Units::cm2m, 1e-12);
    EXPECT_NEAR(zIni, zb * Units::cm2m, 1e-12);
    EXPECT_NEAR(zFinal, ze * Units::cm2m, 1e-12);
}

// ===========================================================================
// Test: getFrequency / setFrequency throw
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, FrequencyThrows) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("freq.map"), zb, ze, nz, rb, re, nr);
    Fieldmap* fm = Fieldmap::getFieldmap(fname);

    EXPECT_THROW(fm->getFrequency(), GeneralOpalException);
    EXPECT_THROW(fm->setFrequency(100.0), GeneralOpalException);
}

// ===========================================================================
// Test: getInfo does not crash
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, GetInfoNoCrash) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("info.map"), zb, ze, nz, rb, re, nr);
    Fieldmap* fm = Fieldmap::getFieldmap(fname);

    Inform msg("test");
    EXPECT_NO_THROW(fm->getInfo(&msg));
}

// ===========================================================================
// Test: Missing file
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, MissingFile) {
    std::string fname = tmpFile("nonexistent.map");
    EXPECT_THROW(Fieldmap::getFieldmap(fname), GeneralOpalException);
}

// ===========================================================================
// Test: Fieldmap dictionary caching – second getFieldmap returns same pointer
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, DictionaryCaching) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("cache.map"), zb, ze, nz, rb, re, nr);
    Fieldmap* fm1 = Fieldmap::getFieldmap(fname);
    Fieldmap* fm2 = Fieldmap::getFieldmap(fname);
    EXPECT_EQ(fm1, fm2);
}

// ===========================================================================
// Test: readMap / freeMap cycle – reading after free should work
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, ReadFreeCycle) {
    // Use Bz_val = 1.0 so normalization is a no-op (uniform: Bzmax = 1.0)
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const int nz = 4, nr = 2;
    const double Bz_val = 1.0;

    std::string fname = writeXZFieldmap(tmpFile("cycle.map"), zb, ze, nz, rb, re, nr,
                                        Bz_val);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    Fieldmap::readMap(fname);

    // Verify field works
    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};
    fm->getFieldstrength(R, E, B);
    EXPECT_NEAR(B[2], Bz_val, 1e-10);

    // Free and re-read
    fm->freeMap();

    fm->readMap();
    B = {0.0, 0.0, 0.0};
    fm->getFieldstrength(R, E, B);
    EXPECT_NEAR(B[2], Bz_val, 1e-10);
}

// ===========================================================================
// Test: On-axis at grid points matches exact data
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, GridPointAccuracy) {
    // Field before normalization: Bz(r,z) = z, Br = r
    // With normalization: Bzmax = zend_m, so normalised Bz = z / zend_m
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const int nz = 10, nr = 5;
    const double zend_m = ze * Units::cm2m;
    const double Bzmax  = zend_m;

    std::string fname = writeVaryingXZFieldmap(tmpFile("grid.map"), zb, ze, nz, rb, re, nr);
    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    Fieldmap::readMap(fname);

    const double hz_m = (ze - zb) * Units::cm2m / nz;

    // Check at each z grid point on axis
    for (int iz = 0; iz < nz; ++iz) {
        double z = zb * Units::cm2m + iz * hz_m;
        Vector_t<double, 3> R = {0.0, 0.0, z};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        fm->getFieldstrength(R, E, B);
        EXPECT_NEAR(B[2], z / Bzmax, 1e-10) << "Mismatch at iz=" << iz << " z=" << z;
    }
}

// ===========================================================================
// Test: ZX orientation gives same results as XZ for identical physical field
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, XZvsZXConsistency) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const int nz = 4, nr = 2;
    const double Bz_val = 2.0;

    std::string fnameXZ = writeXZFieldmap(tmpFile("xz2.map"), zb, ze, nz, rb, re, nr,
                                          Bz_val);
    std::string fnameZX = writeZXFieldmap(tmpFile("zx2.map"), zb, ze, nz, rb, re, nr,
                                          Bz_val);

    Fieldmap* fmXZ = Fieldmap::getFieldmap(fnameXZ);
    Fieldmap* fmZX = Fieldmap::getFieldmap(fnameZX);
    Fieldmap::readMap(fnameXZ);
    Fieldmap::readMap(fnameZX);

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B_xz = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B_zx = {0.0, 0.0, 0.0};

    fmXZ->getFieldstrength(R, E, B_xz);
    fmZX->getFieldstrength(R, E, B_zx);

    EXPECT_NEAR(B_xz[0], B_zx[0], 1e-10);
    EXPECT_NEAR(B_xz[1], B_zx[1], 1e-10);
    EXPECT_NEAR(B_xz[2], B_zx[2], 1e-10);
}

// ===========================================================================
// Test: Field with both Br and Bz components projects correctly in x-y
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, RadialFieldProjection) {
    // Uniform Bz = 1.0 (needed to avoid division by zero in normalization),
    // Br = 1.0.  After normalization: Bzmax = 1.0, so values unchanged.
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const int nz = 4, nr = 2;

    std::string fname = writeXZFieldmap(tmpFile("proj.map"), zb, ze, nz, rb, re, nr,
                                        /*Bz=*/1.0, /*Br=*/1.0);
    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    Fieldmap::readMap(fname);

    // At 45 degrees in x-y plane
    double r = 0.02;  // 2 cm
    double x = r / std::sqrt(2.0);
    double y = r / std::sqrt(2.0);
    Vector_t<double, 3> R = {x, y, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    fm->getFieldstrength(R, E, B);
    // Br = 1.0, projected: Bx = Br * x/R, By = Br * y/R
    double expected_Bxy = 1.0 / std::sqrt(2.0);
    EXPECT_NEAR(B[0], expected_Bxy, 1e-10);
    EXPECT_NEAR(B[1], expected_Bxy, 1e-10);
    EXPECT_NEAR(B[2], 1.0, 1e-10);
}

// ===========================================================================
// Test: Accumulation semantics – B is accumulated, not overwritten
// ===========================================================================
TEST_F(FM2DMagnetoStaticTest, FieldAccumulation) {
    const double zb = 0.0, ze = 10.0;
    const double rb = 0.0, re = 5.0;
    const int nz = 4, nr = 2;
    const double Bz_val = 1.0;

    std::string fname = writeXZFieldmap(tmpFile("accum.map"), zb, ze, nz, rb, re, nr,
                                        Bz_val);
    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 5.0};  // pre-existing B

    fm->getFieldstrength(R, E, B);
    // Bz should be 5.0 (initial) + 1.0 (field) = 6.0
    EXPECT_NEAR(B[2], 6.0, 1e-10);
}
