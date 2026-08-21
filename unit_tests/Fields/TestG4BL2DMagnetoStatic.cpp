/**
 * \file TestG4BL2DMagnetoStatic.cpp
 * \brief Unit tests for the G4beamline `cylinder` fieldmap reader
 *
 * Tests cover:
 * - Header parsing (param / cylinder, key order, optional param line)
 * - Field dimensions, mm -> m conversion, isInside()
 * - Bilinear interpolation and Br projection into x-y
 * - Absolute Tesla: the map is NOT normalized, unlike the other readers
 * - normB from the param line
 * - ZREVERSE: mirror in z, Bz negated, Br unchanged
 * - Rows longer than READ_BUFFER_LENGTH (256 chars)
 * - Rejection of grid / data / extend* / E-field blocks and malformed rows
 * - Missing Br block is zero-filled
 * - Dictionary caching, cache conflict on mismatched ZREVERSE
 * - getFieldDerivative / getFrequency / setFrequency throw
 */

#include <gtest/gtest.h>

#include "Fields/Fieldmap.h"
#include "Fields/G4BL2DMagnetoStatic.h"
#include "Ippl.h"
#include "Physics/Units.h"
#include "Utilities/GeneralOpalException.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>

namespace {

    /// Field sample as a function of the grid indices, in Tesla
    using FieldFn = std::function<double(int /*iz*/, int /*ir*/)>;

    // -----------------------------------------------------------------------
    // Helper: write a G4beamline cylinder map.
    //
    // Grid: nZ points from Z0 in steps of dZ, nR points from 0 in steps of dR,
    // all in mm.  Each block is nZ rows of nR tab-separated values, row index
    // = z, column index = r.
    //
    // paramLine may be empty to omit the optional `param` line entirely.
    // -----------------------------------------------------------------------
    std::string writeCylinderMap(
            const std::string& path, double Z0, int nZ, double dZ, int nR, double dR, FieldFn bz,
            FieldFn br, const std::string& paramLine = "param normB=1. current=1.") {
        std::ofstream f(path);
        if (!paramLine.empty()) {
            f << paramLine << "\n";
        }
        f << "cylinder Z0=" << Z0 << " nZ=" << nZ << " dZ=" << dZ << " nR=" << nR << " dR=" << dR
          << "\n";

        f << "Bz\n";
        for (int iz = 0; iz < nZ; ++iz) {
            for (int ir = 0; ir < nR; ++ir) {
                f << bz(iz, ir) << (ir + 1 == nR ? '\n' : '\t');
            }
        }
        if (br) {
            f << "Br\n";
            for (int iz = 0; iz < nZ; ++iz) {
                for (int ir = 0; ir < nR; ++ir) {
                    f << br(iz, ir) << (ir + 1 == nR ? '\n' : '\t');
                }
            }
        }
        f.close();
        return path;
    }

    /// Uniform value, for use as a FieldFn
    FieldFn constantField(double value) {
        return [value](int, int) {
            return value;
        };
    }

}  // namespace

// ===========================================================================
// Test fixture
// ===========================================================================
class G4BL2DMagnetoStaticTest : public ::testing::Test {
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
        tmpDir_ = std::filesystem::temp_directory_path() / "opalx_g4bl2d_test";
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
// Test: header parsing and field dimensions, mm -> m
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, ParseHeaderAndDimensions) {
    // z from -100 to +100 mm, r from 0 to 40 mm
    std::string fname = writeCylinderMap(
            tmpFile("basic.g4blmap"), -100.0, 5, 50.0, 5, 10.0, constantField(1.0),
            constantField(0.0));

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(fm, nullptr);
    EXPECT_EQ(fm->getType(), TG4BL2DMagnetoStatic);
    EXPECT_FALSE(fm->isZReversed());

    Fieldmap::readMap(fname);

    double zBegin = 0.0, zEnd = 0.0;
    fm->getFieldDimensions(zBegin, zEnd);
    EXPECT_NEAR(zBegin, -0.100, 1e-12);
    EXPECT_NEAR(zEnd, 0.100, 1e-12);

    double xIni = 0.0, xFinal = 0.0, yIni = 0.0, yFinal = 0.0, zIni = 0.0, zFinal = 0.0;
    fm->getFieldDimensions(xIni, xFinal, yIni, yFinal, zIni, zFinal);
    EXPECT_NEAR(xIni, -0.040, 1e-12);
    EXPECT_NEAR(xFinal, 0.040, 1e-12);
    EXPECT_NEAR(yIni, -0.040, 1e-12);
    EXPECT_NEAR(yFinal, 0.040, 1e-12);
    EXPECT_NEAR(zIni, -0.100, 1e-12);
    EXPECT_NEAR(zFinal, 0.100, 1e-12);
}

// ===========================================================================
// Test: the `param` line is optional and key order is free
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, ParamLineOptionalAndKeyOrderFree) {
    {
        std::string fname = writeCylinderMap(
                tmpFile("noparam.g4blmap"), 0.0, 4, 25.0, 3, 10.0, constantField(2.0),
                constantField(0.0), /*paramLine=*/"");
        Fieldmap* fm = Fieldmap::getFieldmap(fname);
        ASSERT_NE(fm, nullptr);
        EXPECT_EQ(fm->getType(), TG4BL2DMagnetoStatic);
    }

    // Same grid, keys shuffled, plus the X0/Y0/tolerance keys the reader ignores
    {
        const std::string path = tmpFile("shuffled.g4blmap");
        std::ofstream f(path);
        f << "cylinder nR=3 dR=10.0 X0=0.0 Y0=0.0 dZ=25.0 nZ=4 Z0=0.0 tolerance=0.01\n";
        f << "Bz\n";
        for (int iz = 0; iz < 4; ++iz) {
            f << "2.0\t2.0\t2.0\n";
        }
        f.close();

        Fieldmap* fm = Fieldmap::getFieldmap(path);
        ASSERT_NE(fm, nullptr);
        Fieldmap::readMap(path);

        double zBegin = 0.0, zEnd = 0.0;
        fm->getFieldDimensions(zBegin, zEnd);
        EXPECT_NEAR(zBegin, 0.0, 1e-12);
        EXPECT_NEAR(zEnd, 0.075, 1e-12);
    }
}

// ===========================================================================
// Test: values are absolute Tesla -- NOT normalized to a 1 T peak
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, NotNormalized) {
    std::string fname = writeCylinderMap(
            tmpFile("abs.g4blmap"), 0.0, 5, 50.0, 4, 10.0, constantField(3.0), constantField(0.0));

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.100};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    EXPECT_FALSE(fm->getFieldstrength(R, E, B));
    // FM2DMagnetoStatic would return 1.0 here after normalizing.
    EXPECT_NEAR(B[2], 3.0, 1e-12);
}

// ===========================================================================
// Test: the param line scales the field by normB / current
//
// G4beamline evaluates B * normB * current_deck / current_file.  The file's
// half of that is applied on load, so KS is exactly the deck-side `current=`.
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, ParamLineScalesField) {
    struct Case {
        const char* name;
        const char* paramLine;
        double factor;
    };
    const Case cases[] = {
            {"normb.g4blmap", "param normB=2. current=1.", 2.0},
            {"current.g4blmap", "param current=4.", 0.25},
            {"both.g4blmap", "param normB=3. current=2.", 1.5},
            {"neither.g4blmap", "param maxline=1024", 1.0},
    };

    for (const Case& c : cases) {
        std::string fname = writeCylinderMap(
                tmpFile(c.name), 0.0, 5, 50.0, 4, 10.0, constantField(1.5), constantField(0.25),
                c.paramLine);

        Fieldmap* fm = Fieldmap::getFieldmap(fname);
        Fieldmap::readMap(fname);

        Vector_t<double, 3> R = {0.010, 0.0, 0.100};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        EXPECT_FALSE(fm->getFieldstrength(R, E, B)) << c.name;
        EXPECT_NEAR(B[2], 1.5 * c.factor, 1e-12) << c.name;
        EXPECT_NEAR(B[0], 0.25 * c.factor, 1e-12) << c.name;
    }
}

// ===========================================================================
// Test: current=0 is rejected rather than dividing by zero
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, ZeroCurrentThrows) {
    std::string fname = writeCylinderMap(
            tmpFile("zerocurrent.g4blmap"), 0.0, 4, 25.0, 3, 10.0, constantField(1.0),
            constantField(0.0), "param normB=1. current=0.");

    EXPECT_THROW(Fieldmap::getFieldmap(fname), GeneralOpalException);
}

// ===========================================================================
// Test: bilinear interpolation and Br projection into x-y
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, InterpolationAndRadialProjection) {
    // Bz = z [m], Br = r [m], on a 0..100 mm x 0..40 mm grid
    const double dZ = 25.0, dR = 10.0;
    auto bz = [dZ](int iz, int) {
        return iz * dZ * Units::mm2m;
    };
    auto br = [dR](int, int ir) {
        return ir * dR * Units::mm2m;
    };

    std::string fname = writeCylinderMap(tmpFile("interp.g4blmap"), 0.0, 5, dZ, 5, dR, bz, br);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    Fieldmap::readMap(fname);

    // On a grid point, on axis
    {
        Vector_t<double, 3> R = {0.0, 0.0, 0.050};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};
        EXPECT_FALSE(fm->getFieldstrength(R, E, B));
        EXPECT_NEAR(B[2], 0.050, 1e-12);
        EXPECT_NEAR(B[0], 0.0, 1e-12);
        EXPECT_NEAR(B[1], 0.0, 1e-12);
    }

    // Between grid points in z, still exact because the field is linear in z
    {
        Vector_t<double, 3> R = {0.0, 0.0, 0.030};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};
        EXPECT_FALSE(fm->getFieldstrength(R, E, B));
        EXPECT_NEAR(B[2], 0.030, 1e-12);
    }

    // Off axis at 45 degrees: Br splits evenly between x and y
    {
        const double r        = 0.020;
        const double xy       = r / std::sqrt(2.0);
        Vector_t<double, 3> R = {xy, xy, 0.050};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};
        EXPECT_FALSE(fm->getFieldstrength(R, E, B));
        EXPECT_NEAR(B[0], r / std::sqrt(2.0), 1e-12);
        EXPECT_NEAR(B[1], r / std::sqrt(2.0), 1e-12);
        EXPECT_NEAR(B[2], 0.050, 1e-12);
    }
}

// ===========================================================================
// Test: isInside() and out-of-bounds behaviour
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, Bounds) {
    std::string fname = writeCylinderMap(
            tmpFile("bounds.g4blmap"), -100.0, 5, 50.0, 5, 10.0, constantField(1.0),
            constantField(0.0));

    auto* fm = dynamic_cast<G4BL2DMagnetoStatic*>(Fieldmap::getFieldmap(fname));
    ASSERT_NE(fm, nullptr);
    Fieldmap::readMap(fname);

    EXPECT_TRUE(fm->isInside({0.0, 0.0, 0.0}));
    EXPECT_TRUE(fm->isInside({0.0, 0.0, -0.100}));   // on zbegin
    EXPECT_FALSE(fm->isInside({0.0, 0.0, 0.100}));   // zend is exclusive
    EXPECT_FALSE(fm->isInside({0.0, 0.0, -0.101}));  // before zbegin
    EXPECT_FALSE(fm->isInside({0.050, 0.0, 0.0}));   // r > rend = 40 mm

    // Outside leaves B untouched and reports true
    Vector_t<double, 3> R = {0.0, 0.0, 0.200};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};
    EXPECT_TRUE(fm->getFieldstrength(R, E, B));
    EXPECT_NEAR(B[2], 0.0, 1e-15);
}

// ===========================================================================
// Test: B is accumulated, not overwritten
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, FieldAccumulation) {
    std::string fname = writeCylinderMap(
            tmpFile("accum.g4blmap"), 0.0, 5, 50.0, 4, 10.0, constantField(1.0),
            constantField(0.0));

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.100};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 5.0};

    fm->getFieldstrength(R, E, B);
    EXPECT_NEAR(B[2], 6.0, 1e-12);
}

// ===========================================================================
// Test: ZREVERSE mirrors in z and negates Bz, leaving Br alone
//
// This is what G4beamline's `place ... rotation=Y180` does to an axisymmetric
// map.  A deliberately z-asymmetric field is needed, otherwise the mirror is
// invisible.
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, ZReverse) {
    // z from 0 to 100 mm.  Bz = 1 + iz, Br = 10 + iz: both rise with z.
    const int nZ = 5, nR = 4;
    auto bz = [](int iz, int) {
        return 1.0 + iz;
    };
    auto br = [](int iz, int) {
        return 10.0 + iz;
    };

    std::string forwardName =
            writeCylinderMap(tmpFile("fwd.g4blmap"), 0.0, nZ, 25.0, nR, 10.0, bz, br);
    std::string reverseName =
            writeCylinderMap(tmpFile("rev.g4blmap"), 0.0, nZ, 25.0, nR, 10.0, bz, br);

    Fieldmap* forward = Fieldmap::getFieldmap(forwardName);
    Fieldmap* reverse = Fieldmap::getFieldmap(reverseName, false, true);
    ASSERT_NE(reverse, nullptr);
    EXPECT_FALSE(forward->isZReversed());
    EXPECT_TRUE(reverse->isZReversed());
    Fieldmap::readMap(forwardName);
    Fieldmap::readMap(reverseName);

    // Bounds are mirrored: 0..100 mm becomes -100..0 mm
    double zBegin = 0.0, zEnd = 0.0;
    reverse->getFieldDimensions(zBegin, zEnd);
    EXPECT_NEAR(zBegin, -0.100, 1e-12);
    EXPECT_NEAR(zEnd, 0.0, 1e-12);

    // Bz(-z) = -Bz(z), Br(-z) = +Br(z).  Sampled at cell centres so that both
    // z and -z stay strictly inside their maps -- zend is exclusive on each.
    for (int iz = 0; iz < nZ - 1; ++iz) {
        const double z = (iz + 0.5) * 25.0 * Units::mm2m;

        Vector_t<double, 3> E  = {0.0, 0.0, 0.0};
        Vector_t<double, 3> Bf = {0.0, 0.0, 0.0}, Br_ = {0.0, 0.0, 0.0};

        ASSERT_FALSE(forward->getFieldstrength({0.010, 0.0, z}, E, Bf));
        ASSERT_FALSE(reverse->getFieldstrength({0.010, 0.0, -z}, E, Br_));

        EXPECT_NEAR(Br_[2], -Bf[2], 1e-12) << "Bz not negated at iz=" << iz;
        EXPECT_NEAR(Br_[0], Bf[0], 1e-12) << "Br changed at iz=" << iz;
    }
}

// ===========================================================================
// Test: the same file cannot be cached in two orientations
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, ZReverseCacheConflictThrows) {
    std::string fname = writeCylinderMap(
            tmpFile("conflict.g4blmap"), 0.0, 4, 25.0, 3, 10.0, constantField(1.0),
            constantField(0.0));

    Fieldmap* first = Fieldmap::getFieldmap(fname, false, false);
    ASSERT_NE(first, nullptr);
    EXPECT_THROW(Fieldmap::getFieldmap(fname, false, true), GeneralOpalException);
    // Same setting is fine
    EXPECT_EQ(Fieldmap::getFieldmap(fname, false, false), first);
}

// ===========================================================================
// Test: ZREVERSE is rejected for map types that do not support it
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, ZReverseRejectedForOtherMapTypes) {
    const std::string path = tmpFile("plain2d.T7");
    std::ofstream f(path);
    f << "2DMagnetoStatic XZ\n";
    f << "0.0 10.0 4\n";
    f << "0.0 5.0 2\n";
    for (int jr = 0; jr <= 2; ++jr) {
        for (int iz = 0; iz <= 4; ++iz) {
            f << "1.0 0.0\n";
        }
    }
    f.close();

    EXPECT_THROW(Fieldmap::getFieldmap(path, false, true), GeneralOpalException);
}

// ===========================================================================
// Test: rows longer than READ_BUFFER_LENGTH (256 chars)
//
// Regression guard.  Fieldmap::getLine() reads through a 256 character buffer,
// so a row of 60 values -- like the 51-value rows of the real muE4 WSX map --
// would be silently truncated if this reader ever went back to the inherited
// line helpers.
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, RowsLongerThanReadBuffer) {
    const int nZ = 3, nR = 60;

    const std::string path = tmpFile("longrows.g4blmap");
    {
        std::ofstream f(path);
        f << "param normB=1. current=1.\n";
        f << "cylinder Z0=0.0 nZ=" << nZ << " dZ=10.0 nR=" << nR << " dR=1.0\n";
        f << "Bz\n";
        for (int iz = 0; iz < nZ; ++iz) {
            for (int ir = 0; ir < nR; ++ir) {
                f << "1.234567e-01" << (ir + 1 == nR ? '\n' : '\t');
            }
        }
        f << "Br\n";
        for (int iz = 0; iz < nZ; ++iz) {
            for (int ir = 0; ir < nR; ++ir) {
                // Br grows with r so a truncated read shows up as a wrong value
                f << (ir * 1.0e-03) << (ir + 1 == nR ? '\n' : '\t');
            }
        }
        f.close();
    }

    // A row is 60 * 13 characters, comfortably past the 256 byte buffer
    {
        std::ifstream check(path);
        std::string line;
        std::getline(check, line);
        std::getline(check, line);
        std::getline(check, line);
        std::getline(check, line);
        ASSERT_GT(line.size(), 256u);
    }

    Fieldmap* fm = Fieldmap::getFieldmap(path);
    ASSERT_NE(fm, nullptr);
    Fieldmap::readMap(path);

    // Sample near the outer edge of the grid, which only a complete row reaches
    const double r        = 55.0 * Units::mm2m;
    Vector_t<double, 3> R = {r, 0.0, 0.010};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    EXPECT_FALSE(fm->getFieldstrength(R, E, B));
    EXPECT_NEAR(B[2], 1.234567e-01, 1e-12);
    EXPECT_NEAR(B[0], 55.0e-03, 1e-12);
}

// ===========================================================================
// Test: a missing Br block is zero-filled
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, MissingBrIsZero) {
    std::string fname = writeCylinderMap(
            tmpFile("bzonly.g4blmap"), 0.0, 5, 50.0, 4, 10.0, constantField(1.0),
            /*br=*/nullptr);

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.010, 0.0, 0.100};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    EXPECT_FALSE(fm->getFieldstrength(R, E, B));
    EXPECT_NEAR(B[2], 1.0, 1e-12);
    EXPECT_NEAR(B[0], 0.0, 1e-12);
    EXPECT_NEAR(B[1], 0.0, 1e-12);
}

// ===========================================================================
// Test: unsupported G4beamline constructs are rejected with a clear error
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, UnsupportedConstructsThrow) {
    // 3D cartesian `grid` map
    {
        const std::string path = tmpFile("grid.g4blmap");
        std::ofstream f(path);
        f << "param current=1.\n";
        f << "grid X0=-10.0 Y0=-10.0 Z0=-10.0 nX=2 nY=2 nZ=2 dX=10.0 dY=10.0 dZ=10.0\n";
        f << "data\n";
        f << "-1.0e+01\t-1.0e+01\t-1.0e+01\t0.0\t0.0\t0.0\n";
        f.close();
        EXPECT_THROW(Fieldmap::getFieldmap(path), GeneralOpalException);
    }

    // one-point-per-row `data` block inside a cylinder map
    {
        const std::string path = tmpFile("cyldata.g4blmap");
        std::ofstream f(path);
        f << "cylinder Z0=0.0 nZ=3 dZ=10.0 nR=3 dR=10.0\n";
        f << "data\n";
        f << "0.0 0.0 0.0 0.0\n";
        f.close();
        EXPECT_THROW(Fieldmap::getFieldmap(path), GeneralOpalException);
    }

    // symmetry declaration
    {
        const std::string path = tmpFile("extend.g4blmap");
        std::ofstream f(path);
        f << "cylinder Z0=0.0 nZ=3 dZ=10.0 nR=3 dR=10.0\n";
        f << "extendZ flip=Br,Bz\n";
        f.close();
        EXPECT_THROW(Fieldmap::getFieldmap(path), GeneralOpalException);
    }

    // electric field block
    {
        const std::string path = tmpFile("efield.g4blmap");
        std::ofstream f(path);
        f << "cylinder Z0=0.0 nZ=2 dZ=10.0 nR=2 dR=10.0\n";
        f << "Bz\n1.0\t1.0\n1.0\t1.0\n";
        f << "Ez\n1.0\t1.0\n1.0\t1.0\n";
        f.close();
        EXPECT_THROW(Fieldmap::getFieldmap(path), GeneralOpalException);
    }

    // no Bz block at all
    {
        const std::string path = tmpFile("nobz.g4blmap");
        std::ofstream f(path);
        f << "cylinder Z0=0.0 nZ=2 dZ=10.0 nR=2 dR=10.0\n";
        f << "Br\n0.0\t1.0\n0.0\t1.0\n";
        f.close();
        EXPECT_THROW(Fieldmap::getFieldmap(path), GeneralOpalException);
    }
}

// ===========================================================================
// Test: malformed rows and headers are rejected
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, MalformedFileThrows) {
    // too few values on a row
    {
        const std::string path = tmpFile("short.g4blmap");
        std::ofstream f(path);
        f << "cylinder Z0=0.0 nZ=2 dZ=10.0 nR=3 dR=10.0\n";
        f << "Bz\n1.0\t1.0\n1.0\t1.0\t1.0\n";
        f.close();
        EXPECT_THROW(Fieldmap::getFieldmap(path), GeneralOpalException);
    }

    // too many values on a row
    {
        const std::string path = tmpFile("long.g4blmap");
        std::ofstream f(path);
        f << "cylinder Z0=0.0 nZ=2 dZ=10.0 nR=2 dR=10.0\n";
        f << "Bz\n1.0\t1.0\t1.0\n1.0\t1.0\n";
        f.close();
        EXPECT_THROW(Fieldmap::getFieldmap(path), GeneralOpalException);
    }

    // block ends early
    {
        const std::string path = tmpFile("truncated.g4blmap");
        std::ofstream f(path);
        f << "cylinder Z0=0.0 nZ=3 dZ=10.0 nR=2 dR=10.0\n";
        f << "Bz\n1.0\t1.0\n1.0\t1.0\n";
        f.close();
        EXPECT_THROW(Fieldmap::getFieldmap(path), GeneralOpalException);
    }

    // cylinder line missing a required key
    {
        const std::string path = tmpFile("nokey.g4blmap");
        std::ofstream f(path);
        f << "cylinder Z0=0.0 nZ=2 dZ=10.0 dR=10.0\n";
        f << "Bz\n1.0\t1.0\n1.0\t1.0\n";
        f.close();
        EXPECT_THROW(Fieldmap::getFieldmap(path), GeneralOpalException);
    }

    // too few grid points to interpolate on
    {
        const std::string path = tmpFile("tiny.g4blmap");
        std::ofstream f(path);
        f << "cylinder Z0=0.0 nZ=1 dZ=10.0 nR=2 dR=10.0\n";
        f << "Bz\n1.0\t1.0\n";
        f.close();
        EXPECT_THROW(Fieldmap::getFieldmap(path), GeneralOpalException);
    }

    // missing file
    EXPECT_THROW(Fieldmap::getFieldmap(tmpFile("nonexistent.g4blmap")), GeneralOpalException);
}

// ===========================================================================
// Test: comments and blank lines are ignored
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, CommentsAndBlankLinesIgnored) {
    const std::string path = tmpFile("comments.g4blmap");
    std::ofstream f(path);
    f << "# a leading comment\n";
    f << "\n";
    f << "param normB=1. current=1.\n";
    f << "cylinder Z0=0.0 nZ=3 dZ=10.0 nR=2 dR=10.0   # inline comment\n";
    f << "Bz\n";
    f << "2.0\t2.0\n";
    f << "\n";
    f << "2.0\t2.0\n";
    f << "# between rows\n";
    f << "2.0\t2.0\n";
    f.close();

    Fieldmap* fm = Fieldmap::getFieldmap(path);
    ASSERT_NE(fm, nullptr);
    Fieldmap::readMap(path);

    Vector_t<double, 3> R = {0.0, 0.0, 0.010};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};
    EXPECT_FALSE(fm->getFieldstrength(R, E, B));
    EXPECT_NEAR(B[2], 2.0, 1e-12);
}

// ===========================================================================
// Test: dictionary caching and the readMap / freeMap cycle
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, CachingAndReadFreeCycle) {
    std::string fname = writeCylinderMap(
            tmpFile("cycle.g4blmap"), 0.0, 5, 50.0, 4, 10.0, constantField(1.0),
            constantField(0.0));

    Fieldmap* fm = Fieldmap::getFieldmap(fname);
    ASSERT_NE(dynamic_cast<G4BL2DMagnetoStatic*>(fm), nullptr);
    EXPECT_EQ(Fieldmap::getFieldmap(fname), fm);

    Fieldmap::readMap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.100};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};
    fm->getFieldstrength(R, E, B);
    EXPECT_NEAR(B[2], 1.0, 1e-12);

    fm->freeMap();
    fm->readMap();

    B = {0.0, 0.0, 0.0};
    fm->getFieldstrength(R, E, B);
    EXPECT_NEAR(B[2], 1.0, 1e-12);
}

// ===========================================================================
// Test: unimplemented interface methods throw, getInfo does not crash
// ===========================================================================
TEST_F(G4BL2DMagnetoStaticTest, UnimplementedMethods) {
    std::string fname = writeCylinderMap(
            tmpFile("misc.g4blmap"), 0.0, 4, 25.0, 3, 10.0, constantField(1.0), constantField(0.0));

    Fieldmap* fm = Fieldmap::getFieldmap(fname);

    Vector_t<double, 3> R = {0.0, 0.0, 0.010};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    EXPECT_THROW(fm->getFieldDerivative(R, E, B, DX), GeneralOpalException);
    EXPECT_THROW(fm->getFrequency(), GeneralOpalException);
    EXPECT_THROW(fm->setFrequency(100.0), GeneralOpalException);

    EXPECT_NO_THROW(fm->swap());

    Inform msg("test");
    EXPECT_NO_THROW(fm->getInfo(&msg));
}
