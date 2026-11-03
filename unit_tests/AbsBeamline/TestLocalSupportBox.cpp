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

#include "AbsBeamline/CyclotronSector.h"
#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "AbsBeamline/ScalingFFAMagnet.h"
#include "AbstractObjects/OpalData.h"
#include "BeamlineCore/BoxRep.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/FieldmapElementRep.h"
#include "BeamlineCore/MonitorRep.h"
#include "BeamlineCore/MultipoleTRep.h"
#include "BeamlineCore/RBendRep.h"
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
        ElementBase* clone() const override { return new DriftWithExtent(*this); }
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
        double getAmplitude() const override { return 0.0; }
        double getFrequency() const override { return 0.0; }
        double getPhase() const override { return 0.0; }
        ElementBase* clone() const override { return new BoxTravelingWave(*this); }
        Geometry& getGeometry() override { return geometry_m; }
        const Geometry& getGeometry() const override { return geometry_m; }

        void setFieldWindow(double start, double startExit, double end, double period) {
            startField_m     = start;
            startExitField_m = startExit;
            endField_m       = end;
            periodLength_m   = period;
        }

        Geometry geometry_m;
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
    std::unique_ptr<ElementBase> copy(drift.clone());
    ASSERT_TRUE(copy->getReach().has_value());
    EXPECT_DOUBLE_EQ(*copy->getReach(), 0.25);
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

TEST_F(LocalSupportBoxTest, RectangularBendWithoutPoleFacesIncludesSagitta) {
    RBendRep bend("R");
    bend.getGeometry() = Geometry::makeRBend(1.0, 0.6);
    bend.setFullGap(0.04);

    const double fringe       = 0.2;
    const double sagitta      = 0.5 * std::tan(0.6 / 4.0);
    const double needed       = reach + sagitta + fringe * std::tan(0.3);
    const auto [lower, upper] = bend.getLocalSupportBox(reach).getCorners();
    EXPECT_GE(upper[0], needed - 1e-12);
    EXPECT_LE(lower[0], -needed + 1e-12);
    EXPECT_DOUBLE_EQ(lower[2], -fringe);
    EXPECT_DOUBLE_EQ(upper[2], 1.0 + fringe);
}

TEST_F(LocalSupportBoxTest, RectangularBendWithPoleFacesContainsDesignPath) {
    RBendRep bend("R");
    Geometry geometry = Geometry::makeRBend(1.0, 0.6);
    geometry.setPoleFaceAnglesEnabled(true);
    geometry.setEntranceAngle(0.3);
    geometry.setExitAngle(0.3);
    bend.getGeometry() = geometry;
    bend.setFullGap(0.04);

    const BoundingBox box = bend.getLocalSupportBox(reach);
    for (const Vector& p : bend.getGeometry().getDesignPath()) {
        for (double dx : {-reach, reach}) {
            EXPECT_TRUE(contains(box, p + Vector({dx, 0.0, 0.0}), 1e-9)) << p;
        }
    }
    // Upstream along the entrance tangent to the start of the fringe at z = -0.2.
    const Vector upstream({-0.2 * std::tan(0.3), 0.0, -0.2});
    EXPECT_TRUE(contains(box, upstream + Vector({-reach, 0.0, 0.0}), 1e-9));
    EXPECT_TRUE(contains(box, upstream + Vector({reach, 0.0, 0.0}), 1e-9));
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
    EXPECT_NEAR(lower[2], -0.3, 1e-12);
    EXPECT_NEAR(upper[2], 0.8, 1e-12);
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

TEST_F(LocalSupportBoxTest, TravelingWaveStopsWhereSelectionStops) {
    BoxTravelingWave cavity;
    cavity.geometry_m.setElementLength(4.0);
    cavity.setFieldWindow(-1.0, 3.0, 4.0, 2.0);

    const auto [lower, upper] = cavity.getLocalSupportBox(reach).getCorners();
    EXPECT_DOUBLE_EQ(lower[2], -1.0);
    EXPECT_DOUBLE_EQ(upper[2], 3.0);
    EXPECT_DOUBLE_EQ(upper[0], reach);
}

TEST_F(LocalSupportBoxTest, MonitorIsCentredOnItsPlane) {
    MonitorRep monitor("M");
    monitor.getGeometry().setElementLength(0.1);

    const auto [lower, upper] = monitor.getLocalSupportBox(reach).getCorners();
    EXPECT_DOUBLE_EQ(lower[2], -0.05);
    EXPECT_DOUBLE_EQ(upper[2], 0.05);
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

TEST_F(LocalSupportBoxTest, BoxElementUsesItsBlock) {
    BoxRep block("BLOCK");
    block.getGeometry().setElementLength(1.0);
    block.setHalfSize(0.9, 0.1);

    const auto [lower, upper] = block.getLocalSupportBox(reach).getCorners();
    EXPECT_DOUBLE_EQ(upper[0], 0.9);
    EXPECT_DOUBLE_EQ(lower[0], -0.9);
    EXPECT_DOUBLE_EQ(upper[1], reach);
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

    const auto [lower, upper] = element.getLocalSupportBox(reach).getCorners();
    EXPECT_LE(lower[0], -0.8 + 1e-12);
    EXPECT_GE(upper[0], 0.8 - 1e-12);
    EXPECT_NEAR(lower[2], 0.1, 1e-12);
    EXPECT_NEAR(upper[2], 0.6, 1e-12);
    std::filesystem::remove(file);
}

TEST_F(LocalSupportBoxTest, ScalingFFABoxContainsSupportForBothSignsOfR0) {
    for (double r0 : {3.0, -3.0}) {
        ScalingFFAMagnet magnet("FFA");
        magnet.setEndField(std::make_shared<endfieldmodel::Tanh>(0.1, 0.02, 20));
        magnet.setTanDelta(std::tan(Physics::pi / 4.0));
        magnet.setR0(r0);
        magnet.setRMin(2.5);
        magnet.setRMax(3.5);
        magnet.setDipoleConstant(1.0);
        magnet.setFieldIndex(15);
        magnet.setMaxOrder(10);
        magnet.setPhiStart(0.2);
        magnet.setAzimuthalExtent(0.3);
        magnet.setVerticalExtent(0.1);
        magnet.initialise();

        const BoundingBox box = magnet.getLocalSupportBox(reach);
        const double sign     = r0 > 0.0 ? 1.0 : -1.0;
        const int n           = 40;
        for (int i = 0; i <= n; ++i) {
            const double r = 2.5 + 1.0 * i / n;
            for (int j = 0; j <= n; ++j) {
                const double psi = -0.3 + 0.6 * j / n;
                const double phi = psi + std::log(r / std::abs(r0)) + 0.2;
                for (double y : {-0.1, 0.1}) {
                    const Vector p({sign * r * std::cos(phi) - r0, y, r * std::sin(phi)});
                    if (i > 0 && i < n && j > 0 && j < n)
                        EXPECT_TRUE(magnet.isInside(p))
                                << "r0=" << r0 << " r=" << r << " psi=" << psi;
                    EXPECT_TRUE(contains(box, p, 1e-9))
                            << "r0=" << r0 << " r=" << r << " psi=" << psi << " point " << p;
                }
            }
        }
    }
}

TEST_F(LocalSupportBoxTest, CyclotronSectorBoxContainsAnnularSector) {
    // Synthetic PSI layout from TestCyclotronSector: r = 1.0 .. 1.4 m, one 45-degree sector.
    const std::string file = temporaryPath("_cyclotron.map");
    {
        std::ofstream out(file);
        out << "1000 100 0 5.625\n";
        for (int i = 0; i < 13; ++i)
            out << "header ";
        out << "5 8 a b c d e 2 a b c d 0.100000000E+01-0.100000000E+01 a b c d e f LREC= a b "
               "c d e\n";
        for (int i = 0; i < 5; ++i) {
            if (i) out << "a b c d e f\n";
            for (int channel = 0; channel < 4; ++channel)
                for (int j = 0; j < 8; ++j)
                    out << (channel == 0 ? 10 * (1 + 0.1 * i) : 0) << ' ';
            out << '\n';
        }
    }
    CyclotronSector sector("SECTOR");
    sector.configure(CyclotronSectorFieldMap::read(file), 8, -0.05, 0.05, 1, {});

    const BoundingBox box = sector.getLocalSupportBox(reach);
    const double centre   = 1.2;
    const int n           = 30;
    for (int i = 0; i <= n; ++i) {
        const double r = 1.0 + 0.4 * i / n;
        for (int j = 0; j < n; ++j) {
            const double theta = 0.25 * Physics::pi * j / n;
            for (double y : {-0.05, 0.05}) {
                const Vector p({r * std::cos(theta) - centre, y, r * std::sin(theta)});
                // The radial limits are inclusive, so round-off decides isInside() there.
                if (i > 0 && i < n)
                    EXPECT_TRUE(sector.isInside(p)) << "r=" << r << " theta=" << theta;
                EXPECT_TRUE(contains(box, p, 1e-9)) << "r=" << r << " theta=" << theta;
            }
        }
    }
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

TEST_F(LocalSupportBoxTest, ScalingFFAPatchBoxesCoverTheSupportTightly) {
    for (double r0 : {3.0, -3.0}) {
        ScalingFFAMagnet magnet("FFA");
        magnet.setEndField(std::make_shared<endfieldmodel::Tanh>(0.1, 0.02, 20));
        magnet.setTanDelta(std::tan(Physics::pi / 4.0));
        magnet.setR0(r0);
        magnet.setRMin(2.5);
        magnet.setRMax(3.5);
        magnet.setDipoleConstant(1.0);
        magnet.setFieldIndex(15);
        magnet.setMaxOrder(10);
        magnet.setPhiStart(0.2);
        magnet.setAzimuthalExtent(0.3);
        magnet.setVerticalExtent(0.1);
        magnet.initialise();

        const auto boxes = magnet.getLocalSupportBoxes(reach);
        ASSERT_GT(boxes.size(), 1u);
        const auto inAny = [&](const Vector& p) {
            return std::any_of(boxes.begin(), boxes.end(), [&](const auto& box) {
                return box.contains(p, 1e-9);
            });
        };
        const double sign = r0 > 0.0 ? 1.0 : -1.0;
        const auto point  = [&](double r, double psi, double y) {
            const double phi = psi + std::log(r / std::abs(r0)) + 0.2;
            return Vector({sign * r * std::cos(phi) - r0, y, r * std::sin(phi)});
        };
        for (int i = 0; i <= 40; ++i)
            for (int j = 0; j <= 40; ++j)
                for (double y : {-0.1, 0.1})
                    EXPECT_TRUE(inAny(point(2.5 + 0.025 * i, -0.3 + 0.015 * j, y)))
                            << "r0=" << r0 << " i=" << i << " j=" << j;

        // On a grid over the bounding box, the patches add little beyond the support itself.
        const auto [lower, upper] = magnet.getLocalSupportBox(reach).getCorners();
        int insideBox = 0, insideSupport = 0, insidePatches = 0;
        for (int i = 0; i <= 40; ++i) {
            for (int k = 0; k <= 40; ++k) {
                const Vector p(
                        {lower(0) + (upper(0) - lower(0)) * i / 40.0, 0.0,
                         lower(2) + (upper(2) - lower(2)) * k / 40.0});
                ++insideBox;
                insideSupport += magnet.isInside(p) ? 1 : 0;
                insidePatches += inAny(p) ? 1 : 0;
            }
        }
        EXPECT_LT(insidePatches, 1.25 * insideSupport) << "r0=" << r0 << " box " << insideBox;
        EXPECT_LT(insidePatches, insideBox) << "r0=" << r0;
    }
}
