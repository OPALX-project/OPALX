/**
 * \file TestG4BL3DGrid.cpp
 * \brief Unit tests for the G4beamline `grid` fieldmap reader
 *
 * Tests cover:
 * - Header parsing (param / grid, optional param line, missing keys)
 * - Field dimensions, mm -> m conversion, isInside()
 * - Trilinear interpolation, exact on a field that is linear in each coordinate
 * - Absolute Tesla: the map is NOT normalized
 * - normB and current from the param line
 * - Rows placed by their own coordinates: repeated, missing and off-grid rows are errors
 * - Rows in a shuffled order still land in the right place
 * - Rejection of extend* / points sections
 * - ZREVERSE rejected (cylinder maps only)
 * - getFieldDerivative / getFrequency / setFrequency throw
 */

#include "Fields/Fieldmap.h"
#include "Fields/G4BL3DGrid.h"
#include "Ippl.h"
#include "Utilities/GeneralOpalException.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <random>
#include <string>
#include <vector>

namespace {

    /// B as a function of position in mm, returning Tesla.
    using FieldFn = std::function<Vector_t<double, 3>(double, double, double)>;

    struct GridSpec {
        double x0 = -10.0, y0 = -20.0, z0 = 0.0;  // mm
        int nx = 3, ny = 4, nz = 5;
        double dx = 10.0, dy = 10.0, dz = 5.0;  // mm
    };

    /// Write a G4beamline `grid` map. Rows come out z fastest, then y, then x -- the order
    /// real files use -- unless `shuffle` is set.
    std::string writeGridMap(
            const std::string& path, const GridSpec& g, const FieldFn& field,
            const std::string& paramLine = "param current=1.", bool shuffle = false) {
        std::vector<std::string> rows;
        rows.reserve(static_cast<size_t>(g.nx) * g.ny * g.nz);

        for (int ix = 0; ix < g.nx; ++ix) {
            for (int iy = 0; iy < g.ny; ++iy) {
                for (int iz = 0; iz < g.nz; ++iz) {
                    const double x = g.x0 + ix * g.dx;
                    const double y = g.y0 + iy * g.dy;
                    const double z = g.z0 + iz * g.dz;
                    const auto B   = field(x, y, z);
                    std::ostringstream row;
                    row << std::setprecision(17) << x << " " << y << " " << z << " " << B(0) << " "
                        << B(1) << " " << B(2);
                    rows.push_back(row.str());
                }
            }
        }

        if (shuffle) {
            std::mt19937 rng(12345);
            std::shuffle(rows.begin(), rows.end(), rng);
        }

        std::ofstream f(path);
        if (!paramLine.empty()) {
            f << paramLine << "\n";
        }
        f << "grid X0=" << g.x0 << " Y0=" << g.y0 << " Z0=" << g.z0 << " nX=" << g.nx
          << " nY=" << g.ny << " nZ=" << g.nz << " dX=" << g.dx << " dY=" << g.dy << " dZ=" << g.dz
          << "\n";
        f << "data\n";
        for (const auto& row : rows) {
            f << row << "\n";
        }
        return path;
    }

    /// Write the nine-column form, x y z Bx By Bz Ex Ey Ez. Magnetic values are Tesla and
    /// electric ones MV/m, which is what G4beamline writes.
    std::string writeGridMapWithE(
            const std::string& path, const GridSpec& g, const FieldFn& bfield,
            const FieldFn& efield, const std::string& paramLine = "param current=1.") {
        std::ofstream f(path);
        if (!paramLine.empty()) {
            f << paramLine << "\n";
        }
        f << "grid X0=" << g.x0 << " Y0=" << g.y0 << " Z0=" << g.z0 << " nX=" << g.nx
          << " nY=" << g.ny << " nZ=" << g.nz << " dX=" << g.dx << " dY=" << g.dy << " dZ=" << g.dz
          << "\n";
        f << "data\n";
        for (int ix = 0; ix < g.nx; ++ix) {
            for (int iy = 0; iy < g.ny; ++iy) {
                for (int iz = 0; iz < g.nz; ++iz) {
                    const double x = g.x0 + ix * g.dx;
                    const double y = g.y0 + iy * g.dy;
                    const double z = g.z0 + iz * g.dz;
                    const auto B   = bfield(x, y, z);
                    const auto E   = efield(x, y, z);
                    f << std::setprecision(17) << x << " " << y << " " << z << " " << B(0) << " "
                      << B(1) << " " << B(2) << " " << E(0) << " " << E(1) << " " << E(2) << "\n";
                }
            }
        }
        return path;
    }

    FieldFn uniformField(double bx, double by, double bz) {
        return [bx, by, bz](double, double, double) {
            return Vector_t<double, 3>(bx, by, bz);
        };
    }

    /// Linear in each coordinate, so trilinear interpolation reproduces it exactly.
    FieldFn linearField() {
        return [](double x, double y, double z) {
            return Vector_t<double, 3>(0.001 * x, 0.002 * y, 0.003 * z);
        };
    }

}  // namespace

class G4BL3DGridTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        std::filesystem::create_directories("data");
    }

    static void TearDownTestSuite() { ippl::finalize(); }

    void SetUp() override {
        const auto* info = ::testing::UnitTest::GetInstance()->current_test_info();
        stem_            = std::string("data/g4bl3d_") + info->name();
    }

    void TearDown() override {
        for (const auto& file : files_) {
            Fieldmap::freeMap(file);
            std::filesystem::remove(file);
        }
        files_.clear();
    }

    std::string path(const std::string& suffix = ".g4blmap") {
        files_.push_back(stem_ + suffix);
        return files_.back();
    }

    /// Load and read a map through the factory, the way an element does.
    Fieldmap* load(const std::string& file) {
        Fieldmap* map = Fieldmap::getFieldmap(file, false, false);
        Fieldmap::readMap(file);
        return map;
    }

    std::string stem_;
    std::vector<std::string> files_;
};

TEST_F(G4BL3DGridTest, DetectedAsAGridMapAndReportsItsBox) {
    const auto file = path();
    GridSpec g;
    writeGridMap(file, g, uniformField(0.0, 0.0, 0.5));

    EXPECT_EQ(Fieldmap::readHeader(file), TG4BL3DGrid);

    Fieldmap* map = load(file);
    ASSERT_NE(map, nullptr);

    double zBegin = 0.0, zEnd = 0.0;
    map->getFieldDimensions(zBegin, zEnd);
    EXPECT_NEAR(zBegin, 0.0, 1e-12);
    EXPECT_NEAR(zEnd, 0.020, 1e-12);  // (5 - 1) * 5 mm

    double x0 = 0, x1 = 0, y0 = 0, y1 = 0, z0 = 0, z1 = 0;
    map->getFieldDimensions(x0, x1, y0, y1, z0, z1);
    EXPECT_NEAR(x0, -0.010, 1e-12);
    EXPECT_NEAR(x1, 0.010, 1e-12);  // -10 + (3 - 1) * 10 mm
    EXPECT_NEAR(y0, -0.020, 1e-12);
    EXPECT_NEAR(y1, 0.010, 1e-12);  // -20 + (4 - 1) * 10 mm
}

// Absolute Tesla: the tabulated value is what comes back, with no rescaling of the peak.
TEST_F(G4BL3DGridTest, ValuesAreAbsoluteTesla) {
    const auto file = path();
    writeGridMap(file, GridSpec{}, uniformField(0.1, -0.2, 0.3));
    Fieldmap* map = load(file);

    Vector_t<double, 3> E(0.0), B(0.0);
    EXPECT_FALSE(map->getFieldstrength(Vector_t<double, 3>(0.0, -0.005, 0.010), E, B));
    EXPECT_NEAR(B(0), 0.1, 1e-12);
    EXPECT_NEAR(B(1), -0.2, 1e-12);
    EXPECT_NEAR(B(2), 0.3, 1e-12);
}

// Trilinear interpolation is exact for a field linear in each coordinate, so the grid
// spacing contributes no error at all. Sampling deliberately off the grid nodes.
TEST_F(G4BL3DGridTest, InterpolationIsExactForALinearField) {
    const auto file = path();
    GridSpec g;
    writeGridMap(file, g, linearField());
    Fieldmap* map = load(file);

    for (double fx = 0.05; fx < 1.0; fx += 0.31) {
        for (double fy = 0.05; fy < 1.0; fy += 0.29) {
            for (double fz = 0.05; fz < 1.0; fz += 0.27) {
                // Somewhere inside the box, in mm then converted.
                const double xmm = g.x0 + fx * (g.nx - 1) * g.dx;
                const double ymm = g.y0 + fy * (g.ny - 1) * g.dy;
                const double zmm = g.z0 + fz * (g.nz - 1) * g.dz;

                Vector_t<double, 3> E(0.0), B(0.0);
                ASSERT_FALSE(map->getFieldstrength(
                        Vector_t<double, 3>(xmm * 1e-3, ymm * 1e-3, zmm * 1e-3), E, B));

                EXPECT_NEAR(B(0), 0.001 * xmm, 1e-12);
                EXPECT_NEAR(B(1), 0.002 * ymm, 1e-12);
                EXPECT_NEAR(B(2), 0.003 * zmm, 1e-12);
            }
        }
    }
}

// Rows are placed by the coordinates they carry, not by where they sit in the file.
TEST_F(G4BL3DGridTest, RowOrderDoesNotMatter) {
    const auto ordered  = path("_ordered.g4blmap");
    const auto shuffled = path("_shuffled.g4blmap");
    GridSpec g;
    writeGridMap(ordered, g, linearField());
    writeGridMap(shuffled, g, linearField(), "param current=1.", /*shuffle=*/true);

    Fieldmap* a = load(ordered);
    Fieldmap* b = load(shuffled);

    const Vector_t<double, 3> R(0.003, -0.007, 0.011);
    Vector_t<double, 3> E(0.0), Ba(0.0), Bb(0.0);
    ASSERT_FALSE(a->getFieldstrength(R, E, Ba));
    ASSERT_FALSE(b->getFieldstrength(R, E, Bb));
    EXPECT_NEAR(Ba(0), Bb(0), 1e-15);
    EXPECT_NEAR(Ba(1), Bb(1), 1e-15);
    EXPECT_NEAR(Ba(2), Bb(2), 1e-15);
}

// normB and current from the file's own param line are folded in on load, so the element's
// scale is exactly G4beamline's deck-side current=.
TEST_F(G4BL3DGridTest, ParamLineScalesTheField) {
    const auto file = path();
    writeGridMap(file, GridSpec{}, uniformField(0.0, 0.0, 1.0), "param normB=3. current=2.");
    Fieldmap* map = load(file);

    Vector_t<double, 3> E(0.0), B(0.0);
    ASSERT_FALSE(map->getFieldstrength(Vector_t<double, 3>(0.0, -0.005, 0.010), E, B));
    EXPECT_NEAR(B(2), 1.5, 1e-12);
}

TEST_F(G4BL3DGridTest, ParamLineIsOptional) {
    const auto file = path();
    writeGridMap(file, GridSpec{}, uniformField(0.0, 0.0, 0.7), "");
    Fieldmap* map = load(file);

    Vector_t<double, 3> E(0.0), B(0.0);
    ASSERT_FALSE(map->getFieldstrength(Vector_t<double, 3>(0.0, -0.005, 0.010), E, B));
    EXPECT_NEAR(B(2), 0.7, 1e-12);
}

// The upper face in each axis has no cell to interpolate in, so it is outside.
TEST_F(G4BL3DGridTest, IsInsideExcludesTheUpperFaces) {
    const auto file = path();
    GridSpec g;
    writeGridMap(file, g, uniformField(0.0, 0.0, 0.5));
    Fieldmap* map = load(file);

    Vector_t<double, 3> E(0.0), B(0.0);
    EXPECT_FALSE(map->getFieldstrength(Vector_t<double, 3>(-0.010, -0.020, 0.0), E, B));
    EXPECT_TRUE(map->getFieldstrength(Vector_t<double, 3>(0.010, 0.0, 0.010), E, B));   // x face
    EXPECT_TRUE(map->getFieldstrength(Vector_t<double, 3>(0.0, 0.010, 0.010), E, B));   // y face
    EXPECT_TRUE(map->getFieldstrength(Vector_t<double, 3>(0.0, 0.0, 0.020), E, B));     // z face
    EXPECT_TRUE(map->getFieldstrength(Vector_t<double, 3>(-0.011, 0.0, 0.010), E, B));  // below
}

TEST_F(G4BL3DGridTest, RepeatedRowIsRejected) {
    const auto file = path();
    writeGridMap(file, GridSpec{}, uniformField(0.0, 0.0, 0.5));
    {
        std::ofstream f(file, std::ios::app);
        f << "-10 -20 0 0 0 0.5\n";
    }
    Fieldmap::getFieldmap(file, false, false);
    EXPECT_THROW(Fieldmap::readMap(file), GeneralOpalException);
}

TEST_F(G4BL3DGridTest, MissingRowIsRejected) {
    const auto file = path();
    writeGridMap(file, GridSpec{}, uniformField(0.0, 0.0, 0.5));

    // Drop the last data row.
    std::vector<std::string> lines;
    {
        std::ifstream in(file);
        std::string line;
        while (std::getline(in, line)) {
            lines.push_back(line);
        }
    }
    lines.pop_back();
    {
        std::ofstream out(file);
        for (const auto& line : lines) {
            out << line << "\n";
        }
    }

    Fieldmap::getFieldmap(file, false, false);
    EXPECT_THROW(Fieldmap::readMap(file), GeneralOpalException);
}

TEST_F(G4BL3DGridTest, OffGridRowIsRejected) {
    const auto file = path();
    writeGridMap(file, GridSpec{}, uniformField(0.0, 0.0, 0.5));
    {
        std::ofstream f(file, std::ios::app);
        f << "-7.3 -20 0 0 0 0.5\n";  // x is not on the 10 mm grid
    }
    Fieldmap::getFieldmap(file, false, false);
    EXPECT_THROW(Fieldmap::readMap(file), GeneralOpalException);
}

TEST_F(G4BL3DGridTest, ShortRowIsRejected) {
    const auto file = path();
    writeGridMap(file, GridSpec{}, uniformField(0.0, 0.0, 0.5));
    {
        std::ofstream f(file, std::ios::app);
        f << "0 0 0 1.0\n";
    }
    Fieldmap::getFieldmap(file, false, false);
    EXPECT_THROW(Fieldmap::readMap(file), GeneralOpalException);
}

TEST_F(G4BL3DGridTest, MissingGridKeyIsRejected) {
    const auto file = path();
    {
        std::ofstream f(file);
        f << "grid X0=0 Y0=0 nX=2 nY=2 nZ=2 dX=1 dY=1 dZ=1\n";  // no Z0
        f << "data\n";
    }
    EXPECT_THROW(Fieldmap::getFieldmap(file, false, false), GeneralOpalException);
}

TEST_F(G4BL3DGridTest, ExtendSectionIsRejectedWithANamedReason) {
    const auto file = path();
    {
        std::ofstream f(file);
        f << "grid X0=0 Y0=0 Z0=0 nX=2 nY=2 nZ=2 dX=1 dY=1 dZ=1\n";
        f << "extendX factor=1.0\n";
    }
    try {
        Fieldmap::getFieldmap(file, false, false);
        FAIL() << "expected a throw";
    } catch (const GeneralOpalException& e) {
        EXPECT_NE(std::string(e.what()).find("extend"), std::string::npos);
    }
}

TEST_F(G4BL3DGridTest, PointsSectionIsRejected) {
    const auto file = path();
    {
        std::ofstream f(file);
        f << "grid X0=0 Y0=0 Z0=0 nX=2 nY=2 nZ=2 dX=1 dY=1 dZ=1\n";
        f << "points\n";
    }
    EXPECT_THROW(Fieldmap::getFieldmap(file, false, false), GeneralOpalException);
}

// ZREVERSE mirrors in z and negates the longitudinal component, which is well defined for an
// axisymmetric map and not for a cartesian one (it would need an x choice too).
TEST_F(G4BL3DGridTest, ZReverseIsRejected) {
    const auto file = path();
    writeGridMap(file, GridSpec{}, uniformField(0.0, 0.0, 0.5));
    EXPECT_THROW(Fieldmap::getFieldmap(file, false, true), GeneralOpalException);
}

TEST_F(G4BL3DGridTest, UnsupportedQueriesThrow) {
    const auto file = path();
    writeGridMap(file, GridSpec{}, uniformField(0.0, 0.0, 0.5));
    Fieldmap* map = load(file);

    Vector_t<double, 3> E(0.0), B(0.0);
    EXPECT_THROW(map->getFieldDerivative(Vector_t<double, 3>(0.0), E, B, DX), GeneralOpalException);
    EXPECT_THROW(map->getFrequency(), GeneralOpalException);
    EXPECT_THROW(map->setFrequency(1.0), GeneralOpalException);
}

// G4beamline also writes a nine-column form with Ex Ey Ez -- the muE4 quadrupole map
// qsm01a_210_track.g4blmap is one, with all-zero E. Accept it, but only when E really is zero.
TEST_F(G4BL3DGridTest, NineColumnRowsWithZeroEfieldAreAccepted) {
    const auto file = path();
    GridSpec g;
    writeGridMap(file, g, uniformField(0.0, 0.0, 0.5));

    // Rewrite every data row with three trailing zero E columns.
    std::vector<std::string> lines;
    {
        std::ifstream in(file);
        std::string line;
        while (std::getline(in, line)) {
            lines.push_back(line);
        }
    }
    {
        std::ofstream out(file);
        for (size_t i = 0; i < lines.size(); ++i) {
            out << lines[i] << (i >= 3 ? " 0.0 0.0 0.0" : "") << "\n";
        }
    }

    Fieldmap* map = load(file);
    Vector_t<double, 3> E(0.0), B(0.0);
    ASSERT_FALSE(map->getFieldstrength(Vector_t<double, 3>(0.0, -0.005, 0.010), E, B));
    EXPECT_NEAR(B(2), 0.5, 1e-12);
}

// Nine columns of zeros are the common case -- most maps written that way have no electric
// field at all -- and no storage should be allocated for them.
TEST_F(G4BL3DGridTest, AllZeroElectricColumnsAllocateNoStorage) {
    const auto file = path();
    writeGridMapWithE(file, GridSpec{}, uniformField(0.0, 0.0, 0.5), uniformField(0.0, 0.0, 0.0));

    Fieldmap* map = load(file);
    auto* grid    = dynamic_cast<G4BL3DGrid*>(map);
    ASSERT_NE(grid, nullptr);
    EXPECT_FALSE(grid->hasEField());

    Vector_t<double, 3> E(0.0), B(0.0);
    ASSERT_FALSE(map->getFieldstrength(Vector_t<double, 3>(0.0, -0.005, 0.010), E, B));
    EXPECT_NEAR(B(2), 0.5, 1e-12);
    EXPECT_EQ(E(0), 0.0);
    EXPECT_EQ(E(1), 0.0);
    EXPECT_EQ(E(2), 0.0);
}

// The electric columns are MV/m in the file and V/m once loaded.
TEST_F(G4BL3DGridTest, ElectricFieldIsReadAndConvertedToVoltsPerMetre) {
    const auto file = path();
    writeGridMapWithE(file, GridSpec{}, uniformField(0.0, 0.0, 0.5), uniformField(0.0, 0.0, 2.5));

    Fieldmap* map = load(file);
    auto* grid    = dynamic_cast<G4BL3DGrid*>(map);
    ASSERT_NE(grid, nullptr);
    EXPECT_TRUE(grid->hasEField());

    Vector_t<double, 3> E(0.0), B(0.0);
    ASSERT_FALSE(map->getFieldstrength(Vector_t<double, 3>(0.0, -0.005, 0.010), E, B));
    EXPECT_NEAR(B(2), 0.5, 1e-12);
    EXPECT_NEAR(E(2), 2.5e6, 1e-6);
}

// Linear in each coordinate, so trilinear interpolation has to reproduce it exactly.
TEST_F(G4BL3DGridTest, ElectricFieldInterpolationIsExactForALinearField) {
    const auto file = path();
    writeGridMapWithE(file, GridSpec{}, uniformField(0.0, 0.0, 0.0), linearField());

    Fieldmap* map = load(file);
    // Mid-cell in every axis: x = -5 mm, y = -15 mm, z = 2.5 mm.
    const Vector_t<double, 3> R(-0.005, -0.015, 0.0025);
    Vector_t<double, 3> E(0.0), B(0.0);
    ASSERT_FALSE(map->getFieldstrength(R, E, B));
    EXPECT_NEAR(E(0), 0.001 * -5.0 * 1e6, 1e-6);
    EXPECT_NEAR(E(1), 0.002 * -15.0 * 1e6, 1e-6);
    EXPECT_NEAR(E(2), 0.003 * 2.5 * 1e6, 1e-6);
}

// G4beamline scales the two fields by separate pairs of keys: normB/current for the
// magnetic one and normE/gradient for the electric one. Neither pair may touch the other.
TEST_F(G4BL3DGridTest, TheTwoFieldsScaleIndependently) {
    const auto file = path();
    writeGridMapWithE(
            file, GridSpec{}, uniformField(0.0, 0.0, 0.5), uniformField(0.0, 0.0, 2.0),
            "param current=2. normB=6. gradient=4. normE=1.");

    Fieldmap* map = load(file);
    Vector_t<double, 3> E(0.0), B(0.0);
    ASSERT_FALSE(map->getFieldstrength(Vector_t<double, 3>(0.0, -0.005, 0.010), E, B));
    EXPECT_NEAR(B(2), 0.5 * 6.0 / 2.0, 1e-12);   // normB / current = 3
    EXPECT_NEAR(E(2), 2.0e6 * 1.0 / 4.0, 1e-6);  // normE / gradient = 1/4
}

// gradient=0 is the electric counterpart of current=0: G4beamline would divide by it.
TEST_F(G4BL3DGridTest, ZeroGradientIsRejected) {
    const auto file = path();
    writeGridMapWithE(
            file, GridSpec{}, uniformField(0.0, 0.0, 0.5), uniformField(0.0, 0.0, 1.0),
            "param current=1. gradient=0.");
    EXPECT_THROW(Fieldmap::getFieldmap(file, false, false), GeneralOpalException);
}

// Seven or eight numbers is a malformed row either way.
TEST_F(G4BL3DGridTest, SevenColumnRowIsRejected) {
    const auto file = path();
    writeGridMap(file, GridSpec{}, uniformField(0.0, 0.0, 0.5));
    {
        std::ofstream f(file, std::ios::app);
        f << "-10 -20 0 0 0 0.5 0.0\n";
    }
    Fieldmap::getFieldmap(file, false, false);
    EXPECT_THROW(Fieldmap::readMap(file), GeneralOpalException);
}
