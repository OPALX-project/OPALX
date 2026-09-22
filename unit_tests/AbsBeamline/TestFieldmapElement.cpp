/**
 * \file TestFieldmapElement.cpp
 * \brief Unit tests for the FIELDMAP element
 *
 * Tests cover:
 * - The element's box comes from the map's extent, not from the deck
 * - The field is the map's field, times SCALE
 * - Leaving the map transversely means no field, and is NOT a loss
 * - A genuine aperture hit still is a loss
 * - markOutsideAperture gates on the field window, not on [0, L]
 * - getSupportEnvelope: aperture first, map extent second, false for a 1D map
 * - Non-magnetostatic maps are rejected
 * - Type name and field extent
 */

#include "AbsBeamline/FieldmapElement.h"
#include "AbstractObjects/OpalData.h"
#include "BeamlineCore/FieldmapElementRep.h"
#include "BeamlineGeometry/Geometry.h"
#include "Fields/Fieldmap.h"
#include "Ippl.h"
#include "PartBunch/BunchStateHandler.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Utilities/GeneralOpalException.h"
#include "Utilities/Options.h"

#include "gtest/gtest.h"

extern Inform* gmsg;

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace {

    /// Uniform Bz over a cylinder, in the G4beamline `cylinder` format. That reader stores
    /// absolute Tesla and does not normalise, so the tabulated value is what comes back and
    /// the expected numbers below are exact.
    std::string writeUniformCylinderMap(
            const std::string& path, double z0_mm, int nZ, double dZ_mm, int nR, double dR_mm,
            double bz, double br = 0.0) {
        std::ofstream f(path);
        f << "param normB=1. current=1.\n";
        f << "cylinder Z0=" << z0_mm << " nZ=" << nZ << " dZ=" << dZ_mm << " nR=" << nR
          << " dR=" << dR_mm << "\n";
        f << "Bz\n";
        for (int iz = 0; iz < nZ; ++iz) {
            for (int ir = 0; ir < nR; ++ir) {
                f << bz << (ir + 1 == nR ? '\n' : '\t');
            }
        }
        f << "Br\n";
        for (int iz = 0; iz < nZ; ++iz) {
            for (int ir = 0; ir < nR; ++ir) {
                f << br << (ir + 1 == nR ? '\n' : '\t');
            }
        }
        return path;
    }

    /// A one-dimensional on-axis map. It has no transverse extent at all: its off-axis field
    /// is an expansion about the axis rather than a tabulated box.
    std::string writeAstra1DMap(
            const std::string& path, const std::vector<double>& z_m,
            const std::vector<double>& bz) {
        std::ofstream f(path);
        f << "AstraMagnetoStatic 8 TRUE\n";
        for (std::size_t i = 0; i < z_m.size(); ++i) {
            f << z_m[i] << " " << bz[i] << "\n";
        }
        return path;
    }

    /// A time-dependent map, which the element must refuse.
    std::string writeFM2DDynamicMap(const std::string& path) {
        std::ofstream f(path);
        f << "2DDynamic XZ\n";
        f << "0.0 10.0 4\n";  // z from/to [cm], intervals
        f << "100.0\n";       // frequency [MHz]
        f << "0.0 3.0 3\n";   // r from/to [cm], intervals
        for (int i = 0; i < 5 * 4; ++i) {
            f << "1.0 0.0 0.0 0.0\n";
        }
        return path;
    }

}  // namespace

class FieldmapElementTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        gmsg                = new Inform(nullptr, -1);
        Options::enableHDF5 = false;
        std::filesystem::create_directories("data");
    }

    static void TearDownTestSuite() {
        delete gmsg;
        gmsg = nullptr;
        ippl::finalize();
    }

    void SetUp() override {
        const auto* info = ::testing::UnitTest::GetInstance()->current_test_info();
        testStem_        = std::string("fieldmapelement_") + info->name();
        OpalData::getInstance()->storeInputFn(testStem_ + ".opal");
        OpalData::getInstance()->setOpenMode(OpalData::OpenMode::WRITE);
    }

    void TearDown() override {
        for (const auto& file : mapFiles_) {
            Fieldmap::freeMap(file);
            std::filesystem::remove(file);
        }
        mapFiles_.clear();
    }

    /// Register a map file so TearDown releases and removes it.
    std::string track(const std::string& path) {
        mapFiles_.push_back(path);
        return path;
    }

    std::string mapPath(const std::string& suffix) const {
        return (std::filesystem::path("data") / (testStem_ + suffix)).string();
    }

    std::string testStem_;
    std::vector<std::string> mapFiles_;
};

// The element's box is the map's box. Nothing about it comes from the deck.
TEST_F(FieldmapElementTest, BoxComesFromTheMap) {
    // z = 100 .. 600 mm, r = 0 .. 40 mm.
    const auto file = track(mapPath("_box.g4blmap"));
    writeUniformCylinderMap(file, 100.0, 11, 50.0, 5, 10.0, 0.25);

    FieldmapElementRep element("FM_BOX");
    element.setFieldMapFN(file);
    ASSERT_EQ(element.getGeometry().getElementLength(), 0.0);  // nothing set at parse time

    element.initialise(nullptr);

    double zBegin = 0.0, zEnd = 0.0;
    element.getFieldExtent(zBegin, zEnd);
    EXPECT_NEAR(zBegin, 0.100, 1e-12);
    EXPECT_NEAR(zEnd, 0.600, 1e-12);
    EXPECT_NEAR(element.getGeometry().getElementLength(), 0.500, 1e-12);

    EXPECT_TRUE(element.hasTransverseExtent());
    double h = 0.0, v = 0.0;
    ASSERT_TRUE(element.getSupportEnvelope(h, v));
    EXPECT_NEAR(h, 0.040, 1e-12);
    EXPECT_NEAR(v, 0.040, 1e-12);
}

// The field is the map's field, and SCALE is a plain multiplier on it -- not a normalised
// strength. The G4beamline reader stores absolute Tesla, so 1.0 gives the tabulated value back.
TEST_F(FieldmapElementTest, FieldIsTheMapTimesScale) {
    const auto file = track(mapPath("_scale.g4blmap"));
    writeUniformCylinderMap(file, 0.0, 11, 100.0, 5, 10.0, 0.25);

    FieldmapElementRep element("FM_SCALE");
    element.setFieldMapFN(file);
    element.setScale(2.0);
    element.initialise(nullptr);
    Fieldmap::readMap(file);

    const Vector_t<double, 3> R(0.0, 0.0, 0.5);
    const Vector_t<double, 3> P(0.0);
    Vector_t<double, 3> E(0.0), B(0.0);
    element.apply(R, P, 0.0, E, B);

    EXPECT_NEAR(B(0), 0.0, 1e-12);
    EXPECT_NEAR(B(1), 0.0, 1e-12);
    EXPECT_NEAR(B(2), 0.5, 1e-12);  // 0.25 T x SCALE 2

    // apply() accumulates, like every other element.
    element.apply(R, P, 0.0, E, B);
    EXPECT_NEAR(B(2), 1.0, 1e-12);
}

// Leaving the map transversely gives no field and is NOT a loss. Reporting it as one would
// make OrbitThreader treat a benign miss as hitting material and add a metre to the path
// length -- which is the bug Solenoid still has.
TEST_F(FieldmapElementTest, OutsideTheMapMeansNoFieldAndNoLoss) {
    const auto file = track(mapPath("_oob.g4blmap"));
    writeUniformCylinderMap(file, 0.0, 11, 100.0, 5, 10.0, 0.25);  // r up to 40 mm

    FieldmapElementRep element("FM_OOB");
    element.setFieldMapFN(file);
    element.initialise(nullptr);
    Fieldmap::readMap(file);

    const Vector_t<double, 3> P(0.0);

    // Well inside the z range, but past the map's radius.
    const Vector_t<double, 3> outside(0.5, 0.0, 0.5);
    Vector_t<double, 3> E(0.0), B(0.0);
    element.apply(outside, P, 0.0, E, B);
    EXPECT_NEAR(B(2), 0.0, 1e-12);

    B = 0.0;
    EXPECT_FALSE(element.applyToReferenceParticle(outside, P, 0.0, E, B));
    EXPECT_NEAR(B(2), 0.0, 1e-12);

    // Inside, the same call does give the field and still reports no loss.
    const Vector_t<double, 3> inside(0.0, 0.0, 0.5);
    B = 0.0;
    EXPECT_FALSE(element.applyToReferenceParticle(inside, P, 0.0, E, B));
    EXPECT_NEAR(B(2), 0.25, 1e-12);
}

// A real aperture hit is still a loss, so the benign case above is not just "never true".
TEST_F(FieldmapElementTest, ApertureHitIsStillALoss) {
    const auto file = track(mapPath("_apert.g4blmap"));
    writeUniformCylinderMap(file, 0.0, 11, 100.0, 5, 10.0, 0.25);

    FieldmapElementRep element("FM_APERT");
    element.setFieldMapFN(file);
    element.setAperture(ApertureType::ELLIPTICAL, {0.01, 0.01});
    element.initialise(nullptr);
    Fieldmap::readMap(file);

    const Vector_t<double, 3> P(0.0);
    Vector_t<double, 3> E(0.0), B(0.0);

    // Inside the field window in z, outside the 10 mm aperture.
    EXPECT_TRUE(
            element.applyToReferenceParticle(Vector_t<double, 3>(0.02, 0.0, 0.5), P, 0.0, E, B));

    // Outside the field window in z, the aperture is not enforced.
    EXPECT_FALSE(
            element.applyToReferenceParticle(Vector_t<double, 3>(0.02, 0.0, 5.0), P, 0.0, E, B));
}

// A configured aperture wins over the map extent, because it is the physical bore.
TEST_F(FieldmapElementTest, ApertureWinsOverTheMapExtentForTheEnvelope) {
    const auto file = track(mapPath("_env.g4blmap"));
    writeUniformCylinderMap(file, 0.0, 11, 100.0, 9, 10.0, 0.25);  // r up to 80 mm

    FieldmapElementRep element("FM_ENV");
    element.setFieldMapFN(file);
    element.setAperture(ApertureType::ELLIPTICAL, {0.03, 0.02});
    element.initialise(nullptr);

    double h = 0.0, v = 0.0;
    ASSERT_TRUE(element.getSupportEnvelope(h, v));
    EXPECT_NEAR(h, 0.03, 1e-12);
    EXPECT_NEAR(v, 0.02, 1e-12);
}

// A one-dimensional map declares no transverse extent. That is not an error: the element
// loads, takes its length from the map, and simply has no envelope to report.
TEST_F(FieldmapElementTest, OneDimensionalMapHasNoTransverseExtent) {
    const auto file = track(mapPath("_1d.astra"));
    writeAstra1DMap(file, {0.0, 0.25, 0.5, 0.75, 1.0}, {1.0, 1.0, 1.0, 1.0, 1.0});

    FieldmapElementRep element("FM_1D");
    element.setFieldMapFN(file);
    element.initialise(nullptr);

    EXPECT_FALSE(element.hasTransverseExtent());
    double h = 0.0, v = 0.0;
    EXPECT_FALSE(element.getSupportEnvelope(h, v));

    // The length still comes from the map.
    EXPECT_NEAR(element.getGeometry().getElementLength(), 1.0, 1e-9);
}

// Time-dependent maps are refused. Their readers ignore the scale argument and their phase
// lives on the element rather than in the map, so RFCAVITY is the element for those.
TEST_F(FieldmapElementTest, DynamicMapIsRejected) {
    const auto file = track(mapPath("_dyn.map"));
    writeFM2DDynamicMap(file);

    FieldmapElementRep element("FM_DYN");
    element.setFieldMapFN(file);
    EXPECT_THROW(element.initialise(nullptr), GeneralOpalException);
}

TEST_F(FieldmapElementTest, MissingFieldMapIsRejected) {
    FieldmapElementRep element("FM_NONE");
    EXPECT_THROW(element.initialise(nullptr), GeneralOpalException);
}

TEST_F(FieldmapElementTest, ReportsItsOwnType) {
    FieldmapElementRep element("FM_TYPE");
    EXPECT_EQ(element.getType(), ElementType::FIELDMAP);
    EXPECT_EQ(element.getTypeString(), "Fieldmap");
}

// isInside uses the field window, so an element whose map does not start at z = 0 is
// selected over the map's own interval rather than over [0, L].
TEST_F(FieldmapElementTest, SelectionFollowsTheOffsetFieldWindow) {
    const auto file = track(mapPath("_offset.g4blmap"));
    writeUniformCylinderMap(file, 200.0, 11, 50.0, 5, 10.0, 0.25);  // z = 200 .. 700 mm

    FieldmapElementRep element("FM_OFFSET");
    element.setFieldMapFN(file);
    element.initialise(nullptr);

    // Inside [0, L) = [0, 0.5) but before the field starts: not selected.
    EXPECT_FALSE(element.isInside(Vector_t<double, 3>(0.0, 0.0, 0.1)));
    // Inside the field window: selected.
    EXPECT_TRUE(element.isInside(Vector_t<double, 3>(0.0, 0.0, 0.3)));
    // Past the field window: not selected.
    EXPECT_FALSE(element.isInside(Vector_t<double, 3>(0.0, 0.0, 0.8)));
}

// markOutsideAperture gates on the field window [zBegin, zEnd), not on the body window [0, L]
// that ElementBase assumes. The two differ whenever the map's z range does not start at zero,
// and this map starts at 200 mm.
TEST_F(FieldmapElementTest, ScrapingFollowsTheOffsetFieldWindow) {
    const auto file = track(mapPath("_scrape.g4blmap"));
    writeUniformCylinderMap(file, 200.0, 11, 50.0, 5, 10.0, 0.25);  // z = 200 .. 700 mm

    FieldmapElementRep element("FM_SCRAPE");
    element.setFieldMapFN(file);
    element.setAperture(ApertureType::ELLIPTICAL, {0.01, 0.01});
    element.setFlagDeleteOnTransverseExit(true);
    element.initialise(nullptr);

    const Vector_t<int, 3> nr(8);
    const Vector_t<double, 3> rmin(-1.0);
    const Vector_t<double, 3> hr(0.25);
    ippl::NDIndex<3> domain;
    for (unsigned d = 0; d < 3; ++d) {
        domain[d] = ippl::Index(nr[d]);
    }
    std::array<bool, 3> decomp{true, true, true};
    ippl::UniformCartesian<double, 3> mesh(domain, hr, rmin);
    ippl::FieldLayout<3> layout(MPI_COMM_WORLD, domain, decomp, false);

    auto pc = std::make_shared<ParticleContainer_t>(mesh, layout);
    pc->setBunchStateHandler(std::make_shared<BunchStateHandler>());
    pc->allocateParticles(3);
    pc->createParticles(3);

    auto R     = pc->R.getView();
    auto hostR = Kokkos::create_mirror_view(R);
    hostR(0)   = Vector_t<double, 3>(0.05, 0.0, 0.100);  // outside the aperture, BEFORE the field
    hostR(1)   = Vector_t<double, 3>(0.05, 0.0, 0.400);  // outside the aperture, inside the field
    hostR(2)   = Vector_t<double, 3>(0.00, 0.0, 0.400);  // inside both
    Kokkos::deep_copy(R, hostR);

    // Only the middle one is inside the field window, so only it is scraped. Gating on [0, L)
    // instead would have caught the first and missed the second.
    EXPECT_EQ(element.markOutsideAperture(pc), 1u);

    auto invalid     = pc->InvalidMask.getView();
    auto hostInvalid = Kokkos::create_mirror_view(invalid);
    Kokkos::deep_copy(hostInvalid, invalid);
    EXPECT_FALSE(hostInvalid(0));
    EXPECT_TRUE(hostInvalid(1));
    EXPECT_FALSE(hostInvalid(2));
}
