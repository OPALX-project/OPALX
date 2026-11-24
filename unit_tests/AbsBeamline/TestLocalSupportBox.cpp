/**
 * \file TestLocalSupportBox.cpp
 * \brief Unit tests for ElementBase::getLocalSupportBox() and its overrides.
 *
 * The box is the element-local region in which an element may act on a particle. It must
 * contain the element's selection support (isInside) and the region where its kernels act,
 * and its transverse half-width is max(aperture, reach), where an aperture of 1e5 m or more
 * counts as none. Curved supports are checked by sampling points through the inverse chart.
 */

#include <gtest/gtest.h>

#include <unistd.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>

#include "Ippl.h"

#include "AbstractObjects/OpalData.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/FieldmapElementRep.h"
#include "BeamlineCore/MonitorRep.h"
#include "BeamlineCore/MultipoleTRep.h"
#include "BeamlineCore/SBendRep.h"
#include "BeamlineCore/VariableRFCavityRep.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"

#define private public
#include "AbsBeamline/TravelingWave.h"
#undef private

extern Inform* gmsg;

namespace {
    using Vector = Vector_t<double, 3>;

    constexpr double reach = 0.5;

    bool contains(const BoundingBox& box, const Vector& p, double tolerance = 1e-12) {
        const auto [lower, upper] = box.getCorners();
        for (unsigned d = 0; d < 3; ++d) {
            if (p[d] < lower[d] - tolerance || p[d] > upper[d] + tolerance) return false;
        }
        return true;
    }

    /// Entrance-frame point with arc coordinates (x, y, s) on a sector arc of curvature h:
    /// the inverse of GeometryHelper::toBendArcCoords.
    Vector arcPoint(double h, double length, double x, double y, double s) {
        if (s <= 0.0) return Vector({x, y, s});
        const double phi = h * std::min(s, length);
        Vector point({(std::cos(phi) - 1.0) / h, y, std::sin(phi) / h});
        const Vector normal({std::cos(phi), 0.0, std::sin(phi)});
        const Vector tangent({-std::sin(phi), 0.0, std::cos(phi)});
        if (s > length) point += (s - length) * tangent;
        return point + x * normal;
    }

    /// Expect that every point of the band |x| <= halfX, |y| <= halfY, s in [sBegin, sEnd]
    /// around the arc lies in the box.
    void expectArcBandInside(
            const BoundingBox& box, double h, double length, double halfX, double halfY,
            double sBegin, double sEnd) {
        const int n = 60;
        for (int i = 0; i <= n; ++i) {
            const double s = sBegin + (sEnd - sBegin) * i / n;
            for (double x : {-halfX, 0.0, halfX}) {
                for (double y : {-halfY, halfY}) {
                    const Vector p = arcPoint(h, length, x, y, s);
                    EXPECT_TRUE(contains(box, p, 1e-9))
                            << "s=" << s << " x=" << x << " y=" << y << " point " << p;
                }
            }
        }
    }

    class DriftWithExtent : public DriftRep {
    public:
        DriftWithExtent(double begin, double end)
            : DriftRep("D_EXTENT"), begin_m(begin), end_m(end) {}
        void getFieldExtent(double& zBegin, double& zEnd) const override {
            zBegin = begin_m;
            zEnd   = end_m;
        }

    private:
        double begin_m;
        double end_m;
    };

    class BoxTravelingWave : public TravelingWave {
    public:
        BoxTravelingWave() : TravelingWave("TW") {}

        void setFieldWindow(double start, double startExit, double end, double period) {
            startField_m     = start;
            startExitField_m = startExit;
            endField_m       = end;
            periodLength_m   = period;
        }
    };

    std::string temporaryPath(const std::string& suffix) {
        return (std::filesystem::temp_directory_path()
                / ("opalx-support-box-" + std::to_string(getpid()) + suffix))
                .string();
    }
}  // namespace

class LocalSupportBoxTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        OpalData::getInstance()->storeInputFn("unit_test.opal");
        if (gmsg == nullptr) gmsg = new Inform(nullptr, -1);
        Options::enableHDF5 = false;
    }

    static void TearDownTestSuite() {
        if (gmsg != nullptr) {
            delete gmsg;
            gmsg = nullptr;
        }
        ippl::finalize();
    }
};

TEST_F(LocalSupportBoxTest, ElementWithoutApertureUsesReach) {
    DriftRep drift("D");
    drift.getGeometry().setElementLength(1.0);

    const auto [lower, upper] = drift.getLocalSupportBox(reach).getCorners();
    EXPECT_DOUBLE_EQ(lower[0], -reach);
    EXPECT_DOUBLE_EQ(upper[0], reach);
    EXPECT_DOUBLE_EQ(lower[1], -reach);
    EXPECT_DOUBLE_EQ(upper[1], reach);
    EXPECT_DOUBLE_EQ(lower[2], 0.0);
    EXPECT_DOUBLE_EQ(upper[2], 1.0);
}

TEST_F(LocalSupportBoxTest, ApertureWiderThanReachSetsHalfWidth) {
    DriftRep drift("D");
    drift.getGeometry().setElementLength(1.0);
    drift.setAperture(ApertureType::RECTANGULAR, {0.8, 0.3});

    const auto [lower, upper] = drift.getLocalSupportBox(reach).getCorners();
    EXPECT_DOUBLE_EQ(upper[0], 0.8);
    EXPECT_DOUBLE_EQ(lower[0], -0.8);
    EXPECT_DOUBLE_EQ(upper[1], reach);
}

TEST_F(LocalSupportBoxTest, ApertureAbove1e5CountsAsNone) {
    DriftRep drift("D");
    drift.getGeometry().setElementLength(1.0);
    drift.setAperture(ApertureType::ELLIPTICAL, {2e5, 2e5});

    const auto [lower, upper] = drift.getLocalSupportBox(reach).getCorners();
    EXPECT_DOUBLE_EQ(upper[0], reach);
    EXPECT_DOUBLE_EQ(upper[1], reach);
}

TEST_F(LocalSupportBoxTest, ZRangeJoinsFieldExtentAndBody) {
    DriftWithExtent wide(-0.2, 1.3);
    wide.getGeometry().setElementLength(1.0);
    auto [lowerWide, upperWide] = wide.getLocalSupportBox(reach).getCorners();
    EXPECT_DOUBLE_EQ(lowerWide[2], -0.2);
    EXPECT_DOUBLE_EQ(upperWide[2], 1.3);

    DriftWithExtent narrow(0.2, 0.5);
    narrow.getGeometry().setElementLength(1.0);
    auto [lowerNarrow, upperNarrow] = narrow.getLocalSupportBox(reach).getCorners();
    EXPECT_DOUBLE_EQ(lowerNarrow[2], 0.0);
    EXPECT_DOUBLE_EQ(upperNarrow[2], 1.0);
}

TEST_F(LocalSupportBoxTest, ReachAttributeIsUnsetByDefaultAndCopied) {
    DriftRep drift("D");
    EXPECT_FALSE(drift.getReach().has_value());

    drift.setReach(0.25);
    const DriftRep copy(drift);
    ASSERT_TRUE(copy.getReach().has_value());
    EXPECT_DOUBLE_EQ(*copy.getReach(), 0.25);
}

TEST_F(LocalSupportBoxTest, SectorBendBoxContainsArcBandAndFringes) {
    for (double h : {0.8, -0.8}) {
        SBendRep bend("B");
        bend.getGeometry() = Geometry::makeSBend(1.5, h);
        bend.setFullGap(0.04);  // fringe reach 5 * gap = 0.2 m per side
        bend.setAperture(ApertureType::RECTANGULAR, {0.05, 0.05});

        const BoundingBox box = bend.getLocalSupportBox(reach);
        expectArcBandInside(box, h, 1.5, reach, reach, -0.2, 1.7);

        // Sanity check of the inverse chart: points near the design arc are selected.
        for (double s : {-0.1, 0.4, 1.2, 1.6})
            EXPECT_TRUE(bend.isInside(arcPoint(h, 1.5, 0.01, 0.01, s))) << "s=" << s;

        // The entrance fringe is straight along z: the reach must not pad it upstream. Only
        // the sub-micrometre margin for the arc sampling may.
        EXPECT_NEAR(box.getCorners().first[2], -0.2, 1e-6);
    }
}

TEST_F(LocalSupportBoxTest, StraightMultipoleTUsesItsHalfApertures) {
    MultipoleTRep magnet("MT");
    magnet.setBendAngle(0.0, false);
    magnet.setElementLength(0.5);
    magnet.setAperture(1.2, 2.0);  // full vertical, full horizontal
    magnet.setFringeField(0.25, 0.05, 0.05);
    magnet.setMaxOrder(5, 10);
    magnet.setTransProfile({1});

    const auto [lower, upper] = magnet.getLocalSupportBox(reach).getCorners();
    EXPECT_NEAR(upper[0], 1.0, 1e-12);
    EXPECT_NEAR(lower[0], -1.0, 1e-12);
    EXPECT_NEAR(upper[1], 0.6, 1e-12);
    EXPECT_NEAR(lower[2], -0.5, 1e-12);
    EXPECT_NEAR(upper[2], 1.0, 1e-12);
}

TEST_F(LocalSupportBoxTest, RotatedMultipoleTCoversItsRotatedAperture) {
    MultipoleTRep magnet("MT");
    magnet.setBendAngle(0.0, false);
    magnet.setElementLength(0.5);
    magnet.setAperture(0.4, 2.0);
    magnet.setFringeField(0.25, 0.05, 0.05);
    magnet.setRotation(Physics::pi / 4.0);
    magnet.setMaxOrder(5, 10);
    magnet.setTransProfile({1});

    const double expected     = std::cos(Physics::pi / 4.0) * (1.0 + 0.2);
    const auto [lower, upper] = magnet.getLocalSupportBox(reach).getCorners();
    EXPECT_GE(upper[0], expected - 1e-12);
    EXPECT_GE(upper[1], expected - 1e-12);
    EXPECT_LE(lower[0], -expected + 1e-12);
}

TEST_F(LocalSupportBoxTest, CurvedMultipoleTBoxContainsArcBand) {
    MultipoleTRep magnet("MT");
    magnet.setBendAngle(0.6, false);
    magnet.setElementLength(1.0);
    magnet.setAperture(0.1, 0.1);
    magnet.setFringeField(0.5, 0.05, 0.05);
    magnet.setMaxOrder(5, 10);
    magnet.setTransProfile({1});

    const double h = magnet.getGeometry().getCurvature();
    ASSERT_NE(h, 0.0);
    expectArcBandInside(magnet.getLocalSupportBox(reach), h, 1.0, reach, reach, -0.3, 1.3);
}

TEST_F(LocalSupportBoxTest, TravelingWaveCoversItsWholeFieldWindow) {
    // The particle kernel acts up to endField_m, past where isInside() stops (startExitField_m).
    BoxTravelingWave cavity;
    cavity.getGeometry().setElementLength(4.0);
    cavity.setFieldWindow(-1.0, 3.0, 4.5, 2.0);

    const auto [lower, upper] = cavity.getLocalSupportBox(reach).getCorners();
    EXPECT_DOUBLE_EQ(lower[2], -1.0);
    EXPECT_DOUBLE_EQ(upper[2], 4.5);
    EXPECT_DOUBLE_EQ(upper[0], reach);
}

TEST_F(LocalSupportBoxTest, MonitorCoversItsPlaneAndItsBody) {
    MonitorRep monitor("M");
    monitor.getGeometry().setElementLength(0.1);

    const auto [lower, upper] = monitor.getLocalSupportBox(reach).getCorners();
    EXPECT_DOUBLE_EQ(lower[2], -0.05);
    EXPECT_DOUBLE_EQ(upper[2], 0.1);
}

TEST_F(LocalSupportBoxTest, VariableRFCavityUsesHalfWidthAndHeight) {
    VariableRfCavityRep cavity("V");
    cavity.getGeometry().setElementLength(1.0);
    cavity.setWidth(1.6);
    cavity.setHeight(0.4);

    const auto [lower, upper] = cavity.getLocalSupportBox(reach).getCorners();
    EXPECT_DOUBLE_EQ(upper[0], 0.8);
    EXPECT_DOUBLE_EQ(upper[1], reach);
    EXPECT_DOUBLE_EQ(upper[2], 1.0);
}

TEST_F(LocalSupportBoxTest, FieldmapBoxCoversTheMapCrossSection) {
    // G4beamline cylinder: z = 100 .. 600 mm, r = 0 .. 800 mm.
    const std::string file = temporaryPath("_cylinder.g4blmap");
    {
        std::ofstream f(file);
        f << "param normB=1. current=1.\n";
        f << "cylinder Z0=100 nZ=11 dZ=50 nR=5 dR=200\n";
        for (const char* component : {"Bz", "Br"}) {
            f << component << "\n";
            for (int iz = 0; iz < 11; ++iz)
                for (int ir = 0; ir < 5; ++ir)
                    f << 0.1 << (ir + 1 == 5 ? '\n' : '\t');
        }
    }
    FieldmapElementRep element("FM");
    element.setFieldMapFN(file);
    element.initialise(nullptr);

    // Without an APERTURE the box is the map itself, wider than the reach.
    const auto [lower, upper] = element.getLocalSupportBox(reach).getCorners();
    EXPECT_NEAR(lower[0], -0.8, 1e-12);
    EXPECT_NEAR(upper[0], 0.8, 1e-12);
    EXPECT_NEAR(lower[2], 0.1, 1e-12);
    EXPECT_NEAR(upper[2], 0.6, 1e-12);
    std::filesystem::remove(file);
}

TEST_F(LocalSupportBoxTest, SectorBendSegmentBoxesCoverTheBandButNotTheTangentLine) {
    for (double h : {0.8, -0.8}) {
        SBendRep bend("B");
        bend.getGeometry() = Geometry::makeSBend(1.5, h);
        bend.setFullGap(0.04);
        bend.setAperture(ApertureType::RECTANGULAR, {0.05, 0.05});
        const auto boxes = bend.getLocalSupportBoxes(0.1);
        ASSERT_GT(boxes.size(), 3u);

        const auto inAny = [&](const Vector& p) {
            return std::any_of(boxes.begin(), boxes.end(), [&](const auto& box) {
                return box.contains(p, 1e-9);
            });
        };
        for (int i = 0; i <= 200; ++i) {
            const double s = -0.2 + 1.9 * i / 200;
            for (double x : {-0.1, 0.0, 0.1})
                for (double y : {-0.1, 0.1})
                    EXPECT_TRUE(inAny(arcPoint(h, 1.5, x, y, s))) << "h=" << h << " s=" << s;
        }
        // The entrance tangent at z = 1 m is about 0.27 m off the arc.
        EXPECT_FALSE(inAny(Vector({0.0, 0.0, 1.0}))) << "h=" << h;
    }
}

TEST_F(LocalSupportBoxTest, StraightElementHasOneAxisAlignedSegmentBox) {
    DriftRep drift("D");
    drift.getGeometry().setElementLength(1.0);
    const auto boxes = drift.getLocalSupportBoxes(reach);
    ASSERT_EQ(boxes.size(), 1u);
    EXPECT_TRUE(boxes.front().contains(Vector({0.49, -0.49, 0.99})));
    EXPECT_FALSE(boxes.front().contains(Vector({0.51, 0.0, 0.5})));
}
