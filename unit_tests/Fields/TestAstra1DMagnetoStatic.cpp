/**
 * ile TestAstra1DMagnetoStatic.cpp
 * rief Unit tests for Astra1DMagnetoStatic fieldmap class.
 *
 * Tests cover:
 *
 * Parsing & metadata:
 * - File parsing
 * - Field dimension retrieval
 *
 * Field evaluation:
 * - Point field evaluation via getFieldstrength()
 * - Uniform on-axis Bz correctness
 * - Static computeField() correctness for simple coefficients
 * - Off-axis field expansion for a simple Fourier mode
 *
 * Geometry & bounds:
 * - Behaviour outside z-range
 * - isInside() helper boundary semantics
 *
 * API behaviour:
 * - getFrequency() / setFrequency() throw
 * - getFieldDerivative() current no-op behaviour
 * - getFieldDimensions(6-arg) throws
 * - getInfo() does not crash
 * - swap() does not crash
 *
 * Robustness:
 * - Missing file handling
 * - Fieldmap dictionary caching
 * - readMap() / freeMap() lifecycle correctness
 * - Rejection of too few samples
 * - Rejection of zero-amplitude maps
 *
 * Semantics:
 * - Field accumulation: B is accumulated, not overwritten
 */

#include <gtest/gtest.h>

#include "Fields/Astra1DMagnetoStatic.h"
#include "Fields/Fieldmap.h"
#include "Ippl.h"
#include "Physics/Physics.h"
#include "Utilities/GeneralOpalException.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace {

    constexpr const char* kAstra1DMagnetoStaticType = "AstraMagnetoStatic";

    std::string writeAstra1DMagnetoStaticFieldmap(
            const std::string& path, const std::vector<double>& z_m, const std::vector<double>& bz,
            int accuracy = 8, bool normalize = true) {
        EXPECT_EQ(z_m.size(), bz.size());

        std::ofstream f(path);

        f << kAstra1DMagnetoStaticType << " " << accuracy << " " << (normalize ? "TRUE" : "FALSE")
          << "\n";

        for (std::size_t i = 0; i < z_m.size(); ++i) {
            f << std::setprecision(16) << z_m[i] << " " << bz[i] << "\n";
        }

        return path;
    }

    std::string writeConstantAstra1DMagnetoStaticFieldmap(
            const std::string& path, double zbegin_m, double zend_m, int nz, double bz,
            int accuracy = 8, bool normalize = true) {
        std::vector<double> z(nz);
        std::vector<double> b(nz, bz);

        const double dz = (zend_m - zbegin_m) / double(nz - 1);

        for (int i = 0; i < nz; ++i) {
            z[i] = zbegin_m + i * dz;
        }

        return writeAstra1DMagnetoStaticFieldmap(path, z, b, accuracy, normalize);
    }

    std::string writeRampAstra1DMagnetoStaticFieldmap(
            const std::string& path, double zbegin_m, double zend_m, int nz, int accuracy = 8,
            bool normalize = false) {
        std::vector<double> z(nz);
        std::vector<double> b(nz);

        const double dz = (zend_m - zbegin_m) / double(nz - 1);

        for (int i = 0; i < nz; ++i) {
            z[i] = zbegin_m + i * dz;
            b[i] = 1.0 + z[i] / zend_m;
        }

        return writeAstra1DMagnetoStaticFieldmap(path, z, b, accuracy, normalize);
    }

}  // namespace

class Astra1DMagnetoStaticTest : public ::testing::Test {
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
        tmpDir_ = std::filesystem::temp_directory_path() / "opalx_astra1dmagnetostatic_test";

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
TEST_F(Astra1DMagnetoStaticTest, ParseAndDimensions) {
    const double zb = 0.0;
    const double ze = 0.10;
    const int nz    = 5;

    std::string fname =
            writeConstantAstra1DMagnetoStaticFieldmap(tmpFile("parse.map"), zb, ze, nz, 1.0);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);

    ASSERT_NE(fm, nullptr);
    EXPECT_EQ(fm->getType(), TAstraMagnetoStatic);

    Fieldmap::readMap(fname);

    double zBegin = 0.0;
    double zEnd   = 0.0;

    fm->getFieldDimensions(zBegin, zEnd);

    EXPECT_NEAR(zBegin, zb, 1e-12);
    EXPECT_NEAR(zEnd, ze, 1e-12);
}

// ===========================================================================
// Test: isInside() helper
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, IsInside) {
    std::string fname =
            writeConstantAstra1DMagnetoStaticFieldmap(tmpFile("inside.map"), 0.0, 0.10, 5, 1.0);

    auto* fm = dynamic_cast<Astra1DMagnetoStatic*>(Fieldmap::getFieldmap(fname));

    ASSERT_NE(fm, nullptr);

    EXPECT_TRUE(fm->isInside({0.0, 0.0, 0.00}));
    EXPECT_TRUE(fm->isInside({0.0, 0.0, 0.05}));
    EXPECT_FALSE(fm->isInside({0.0, 0.0, 0.10}));
    EXPECT_FALSE(fm->isInside({0.0, 0.0, -0.01}));
}

// ===========================================================================
// Test: Uniform field with normalization, on-axis Bz
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, UniformFieldStrengthOnAxis) {
    const int nz = 9;

    std::string fname = writeConstantAstra1DMagnetoStaticFieldmap(
            tmpFile("uniform.map"), 0.0, 0.10, nz, 3.0, 8, true);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);

    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {1.0, 2.0, 3.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool outside = fm->getFieldstrength(R, E, B);

    EXPECT_FALSE(outside);

    EXPECT_NEAR(E[0], 1.0, 1e-15);
    EXPECT_NEAR(E[1], 2.0, 1e-15);
    EXPECT_NEAR(E[2], 3.0, 1e-15);

    EXPECT_NEAR(B[0], 0.0, 1e-10);
    EXPECT_NEAR(B[1], 0.0, 1e-10);
    EXPECT_NEAR(B[2], 1.0, 1e-10);
}

// ===========================================================================
// Test: Uniform field without normalization
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, NoNormalization) {
    const int nz = 9;

    std::string fname = writeConstantAstra1DMagnetoStaticFieldmap(
            tmpFile("nonorm.map"), 0.0, 0.10, nz, 3.0, 8, false);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);

    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool outside = fm->getFieldstrength(R, E, B);

    EXPECT_FALSE(outside);

    EXPECT_NEAR(B[0], 0.0, 1e-10);
    EXPECT_NEAR(B[1], 0.0, 1e-10);
    EXPECT_NEAR(B[2], 3.0, 1e-10);
}

// ===========================================================================
// Test: Outside z range returns outside flag and does not modify fields
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, OutsideZRange) {
    std::string fname =
            writeConstantAstra1DMagnetoStaticFieldmap(tmpFile("outside.map"), 0.0, 0.10, 9, 1.0);

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
// Test: Accumulation semantics — B is accumulated, not overwritten
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, FieldAccumulation) {
    std::string fname =
            writeConstantAstra1DMagnetoStaticFieldmap(tmpFile("accum.map"), 0.0, 0.10, 9, 1.0);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);

    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {4.0, 5.0, 6.0};

    fm->getFieldstrength(R, E, B);

    EXPECT_NEAR(B[0], 4.0, 1e-10);
    EXPECT_NEAR(B[1], 5.0, 1e-10);
    EXPECT_NEAR(B[2], 7.0, 1e-10);
}

// ===========================================================================
// Test: computeField static function directly with a constant coefficient
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, ComputeFieldStaticConstantCoefficient) {
    Kokkos::View<double*, Kokkos::HostSpace> coefs("coefs", 3);

    coefs(0) = 1.0;
    coefs(1) = 0.0;
    coefs(2) = 0.0;

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    Astra1DMagnetoStatic::computeField(R, E, B, coefs, 0.0, 0.20, 2);

    EXPECT_NEAR(E[0], 0.0, 1e-12);
    EXPECT_NEAR(E[1], 0.0, 1e-12);
    EXPECT_NEAR(E[2], 0.0, 1e-12);

    EXPECT_NEAR(B[0], 0.0, 1e-12);
    EXPECT_NEAR(B[1], 0.0, 1e-12);
    EXPECT_NEAR(B[2], 1.0, 1e-12);
}

// ===========================================================================
// Test: computeField off-axis with a constant coefficient
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, ComputeFieldStaticOffAxisConstantCoefficient) {
    Kokkos::View<double*, Kokkos::HostSpace> coefs("coefs", 3);

    coefs(0) = 2.0;
    coefs(1) = 0.0;
    coefs(2) = 0.0;

    Vector_t<double, 3> R = {0.01, 0.02, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    Astra1DMagnetoStatic::computeField(R, E, B, coefs, 0.0, 0.20, 2);

    EXPECT_NEAR(B[0], 0.0, 1e-12);
    EXPECT_NEAR(B[1], 0.0, 1e-12);
    EXPECT_NEAR(B[2], 2.0, 1e-12);
}

// ===========================================================================
// Test: computeField off-axis with one Fourier mode
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, ComputeFieldStaticSingleModeOffAxis) {
    Kokkos::View<double*, Kokkos::HostSpace> coefs("coefs", 3);

    coefs(0) = 0.0;
    coefs(1) = 1.0;
    coefs(2) = 0.0;

    const double zbegin = 0.0;
    const double length = 0.20;
    const int accuracy  = 2;
    const double base   = Physics::two_pi / length;

    Vector_t<double, 3> R = {0.01, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    Astra1DMagnetoStatic::computeField(R, E, B, coefs, zbegin, length, accuracy);

    const double rr2   = R(0) * R(0) + R(1) * R(1);
    const double bzp   = base;
    const double bzpp  = 0.0;
    const double bzppp = -base * base * base;

    const double expectedBfieldR = -bzp / 2.0 + bzppp * rr2 / 16.0;

    EXPECT_NEAR(B[0], expectedBfieldR * R(0), 1e-12);
    EXPECT_NEAR(B[1], 0.0, 1e-12);
    EXPECT_NEAR(B[2], 0.0 - bzpp * rr2 / 4.0, 1e-12);
}

// ===========================================================================
// Test: Non-uniform fieldmap produces different Bz values at different z
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, NonUniformMapProducesNonUniformField) {
    std::string fname =
            writeRampAstra1DMagnetoStaticFieldmap(tmpFile("ramp.map"), 0.0, 0.10, 17, 12, false);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);

    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    Vector_t<double, 3> E1 = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B1 = {0.0, 0.0, 0.0};
    Vector_t<double, 3> E2 = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B2 = {0.0, 0.0, 0.0};

    fm->getFieldstrength({0.0, 0.0, 0.025}, E1, B1);
    fm->getFieldstrength({0.0, 0.0, 0.075}, E2, B2);

    EXPECT_NE(B1[2], B2[2]);
}

// ===========================================================================
// Test: getFieldDerivative current no-op behaviour
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, FieldDerivativeNoOp) {
    std::string fname =
            writeConstantAstra1DMagnetoStaticFieldmap(tmpFile("deriv.map"), 0.0, 0.10, 9, 2.0);

    auto* fm = dynamic_cast<Astra1DMagnetoStatic*>(Fieldmap::getFieldmap(fname));

    ASSERT_NE(fm, nullptr);

    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {1.0, 2.0, 3.0};
    Vector_t<double, 3> B = {4.0, 5.0, 6.0};

    bool outside = fm->getFieldDerivative(R, E, B, DZ);

    EXPECT_FALSE(outside);

    EXPECT_NEAR(E[0], 1.0, 1e-15);
    EXPECT_NEAR(E[1], 2.0, 1e-15);
    EXPECT_NEAR(E[2], 3.0, 1e-15);

    EXPECT_NEAR(B[0], 4.0, 1e-15);
    EXPECT_NEAR(B[1], 5.0, 1e-15);
    EXPECT_NEAR(B[2], 6.0, 1e-15);
}

// ===========================================================================
// Test: getFrequency / setFrequency throw for magnetostatic fieldmap
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, FrequencyThrows) {
    std::string fname =
            writeConstantAstra1DMagnetoStaticFieldmap(tmpFile("freq.map"), 0.0, 0.10, 5, 1.0);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);

    ASSERT_NE(fm, nullptr);

    EXPECT_THROW(fm->getFrequency(), GeneralOpalException);
    EXPECT_THROW(fm->setFrequency(100.0), GeneralOpalException);
}

// ===========================================================================
// Test: getFieldDimensions 6-arg overload throws
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, GetFieldDimensions6ArgThrows) {
    std::string fname =
            writeConstantAstra1DMagnetoStaticFieldmap(tmpFile("dim6.map"), 0.0, 0.10, 5, 1.0);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);

    ASSERT_NE(fm, nullptr);

    double xIni   = 0.0;
    double xFinal = 0.0;
    double yIni   = 0.0;
    double yFinal = 0.0;
    double zIni   = 0.0;
    double zFinal = 0.0;

    EXPECT_THROW(
            fm->getFieldDimensions(xIni, xFinal, yIni, yFinal, zIni, zFinal), GeneralOpalException);
}

// ===========================================================================
// Test: swap() does not crash
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, SwapNoCrash) {
    std::string fname =
            writeConstantAstra1DMagnetoStaticFieldmap(tmpFile("swap.map"), 0.0, 0.10, 5, 1.0);

    auto* fm = dynamic_cast<Astra1DMagnetoStatic*>(Fieldmap::getFieldmap(fname));

    ASSERT_NE(fm, nullptr);

    EXPECT_NO_THROW(fm->swap());
}

// ===========================================================================
// Test: getInfo() does not crash
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, GetInfoNoCrash) {
    std::string fname =
            writeConstantAstra1DMagnetoStaticFieldmap(tmpFile("info.map"), 0.0, 0.10, 5, 1.0);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);

    ASSERT_NE(fm, nullptr);

    Inform msg("test");

    EXPECT_NO_THROW(fm->getInfo(&msg));
}

// ===========================================================================
// Test: Missing file
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, MissingFile) {
    std::string fname = tmpFile("nonexistent.map");

    EXPECT_THROW(Fieldmap::getFieldmap(fname), GeneralOpalException);
}

// ===========================================================================
// Test: Fieldmap dictionary caching
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, DictionaryCaching) {
    std::string fname =
            writeConstantAstra1DMagnetoStaticFieldmap(tmpFile("cache.map"), 0.0, 0.10, 5, 1.0);

    Fieldmap* fm1 = Fieldmap::getFieldmap(fname);
    Fieldmap* fm2 = Fieldmap::getFieldmap(fname);

    EXPECT_EQ(fm1, fm2);
}

// ===========================================================================
// Test: readMap / freeMap cycle
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, ReadFreeCycle) {
    std::string fname =
            writeConstantAstra1DMagnetoStaticFieldmap(tmpFile("cycle.map"), 0.0, 0.10, 9, 1.0);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);

    ASSERT_NE(fm, nullptr);

    fm->readMap();

    Vector_t<double, 3> R = {0.0, 0.0, 0.05};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    fm->getFieldstrength(R, E, B);

    EXPECT_NEAR(B[2], 1.0, 1e-10);

    fm->freeMap();
    fm->readMap();

    B = {0.0, 0.0, 0.0};

    fm->getFieldstrength(R, E, B);

    EXPECT_NEAR(B[2], 1.0, 1e-10);
}

// ===========================================================================
// Test: Invalid fieldmap with fewer than two valid samples
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, RejectsTooFewSamples) {
    const std::string fname = tmpFile("too_few.map");

    {
        std::ofstream f(fname);
        f << kAstra1DMagnetoStaticType << " 8 TRUE\n";
        f << "0.0 1.0\n";
    }

    EXPECT_THROW(Fieldmap::getFieldmap(fname), GeneralOpalException);
}

// ===========================================================================
// Test: Zero-amplitude map is rejected during readMap()
// ===========================================================================
TEST_F(Astra1DMagnetoStaticTest, RejectsZeroAmplitudeMap) {
    std::string fname =
            writeConstantAstra1DMagnetoStaticFieldmap(tmpFile("zero.map"), 0.0, 0.10, 5, 0.0);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);

    ASSERT_NE(fm, nullptr);

    EXPECT_THROW(Fieldmap::readMap(fname), GeneralOpalException);
}
