#include <gtest/gtest.h>

#include "BeamlineCore/RFCavityRep.h"
#include "Fields/FM3DDynamic.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Units.h"
#include "Structure/Beam.h"
#include "Utilities/GeneralOpalException.h"

#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <sstream>

namespace {
    using Vector  = Vector_t<double, 3>;
    using Samples = std::array<double, 6>;

    Samples analyticField(double x, double y, double z) {
        return {1 + 2 * x - 3 * y + 4 * z + x * y * z,
                -2 + x + y * z,
                2 + 10 * z + 100 * x + 30 * y + x * y,
                4 + x * y - 2 * z,
                -5 + y + x * z,
                6 - x + 3 * y * z};
    }
}  // namespace

class FM3DDynamicTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        gmsg = new Inform(nullptr, -1);
    }

    static void TearDownTestSuite() {
        Fieldmap::clearDictionary();
        delete gmsg;
        gmsg = nullptr;
        ippl::finalize();
    }

    void SetUp() override {
        const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        directory        = std::filesystem::temp_directory_path()
                    / ("opalx_fm3dynamic_" + std::to_string(stamp) + "_"
                       + std::to_string(ippl::Comm->rank()));
        std::filesystem::create_directories(directory);
    }

    void TearDown() override {
        Fieldmap::clearDictionary();
        std::filesystem::remove_all(directory);
    }

    std::string write(const std::string& name, const std::string& contents) {
        const auto path = (directory / name).string();
        std::ofstream file(path);
        file << contents;
        return path;
    }

    std::string writeMap(
            const std::string& name = "field.map", const std::string& flag = "FALSE",
            const std::function<Samples(double, double, double)>& sample = analyticField,
            double xBegin = -0.02, double xEnd = 0.04, int xIntervals = 2) {
        std::ostringstream file;
        file << std::setprecision(17) << "3DDynamic " << flag << " # amplitudes\n\n"
             << "# frequency in MHz\n125\n"
             << xBegin * 100 << ' ' << xEnd * 100 << ' ' << xIntervals << "\n"
             << "-3 3 3\n-5 15 4\n";
        for (int indexX = 0; indexX <= xIntervals; ++indexX) {
            const double x = xBegin + (xEnd - xBegin) * indexX / xIntervals;
            for (int indexY = 0; indexY <= 3; ++indexY) {
                const double y = -0.03 + 0.02 * indexY;
                for (int indexZ = 0; indexZ <= 4; ++indexZ) {
                    const double z    = -0.05 + 0.05 * indexZ;
                    const auto values = sample(x, y, z);
                    for (const double value : values)
                        file << value << ' ';
                    file << "# grid point\n";
                }
            }
        }
        file << "# trailing comment without a newline";
        return write(name, file.str());
    }

    Fieldmap* load(const std::string& path) {
        auto* map = Fieldmap::getFieldmap(path);
        Fieldmap::readMap(path);
        return map;
    }

    void expectField(Fieldmap* map, const Vector& position, double normalization = 1.0) {
        Vector electric(7.0), magnetic(-0.25);
        ASSERT_FALSE(map->getFieldstrength(position, electric, magnetic));
        const auto expected = analyticField(position(0), position(1), position(2));
        for (unsigned component = 0; component < 3; ++component) {
            EXPECT_NEAR(electric(component), 7.0 + expected[component] * 1e6 / normalization, 2e-8);
            EXPECT_NEAR(
                    magnetic(component),
                    -0.25 + expected[component + 3] * Physics::mu_0 / normalization, 1e-15);
        }
    }

    std::filesystem::path directory;
};

TEST_F(FM3DDynamicTest, FactoryMetadataAndUnloadedAccess) {
    const auto path = writeMap();
    auto* map       = Fieldmap::getFieldmap(path);
    ASSERT_NE(dynamic_cast<FM3DDynamic*>(map), nullptr);
    EXPECT_EQ(map, Fieldmap::getFieldmap(path));
    EXPECT_DOUBLE_EQ(map->getFrequency(), 125 * Physics::two_pi * Units::MHz2Hz);
    double xBegin, xEnd, yBegin, yEnd, zBegin, zEnd;
    map->getFieldDimensions(xBegin, xEnd, yBegin, yEnd, zBegin, zEnd);
    EXPECT_DOUBLE_EQ(xBegin, -0.02);
    EXPECT_DOUBLE_EQ(xEnd, 0.04);
    EXPECT_DOUBLE_EQ(yBegin, -0.03);
    EXPECT_DOUBLE_EQ(yEnd, 0.03);
    EXPECT_DOUBLE_EQ(zBegin, -0.05);
    EXPECT_DOUBLE_EQ(zEnd, 0.15);
    Vector electric(0.0), magnetic(0.0);
    EXPECT_THROW(map->getFieldstrength(Vector(0.0), electric, magnetic), GeneralOpalException);
    Fieldmap::readMap(path);
    Fieldmap::readMap(path);
    expectField(map, Vector(0.0));
}

TEST_F(FM3DDynamicTest, AllComponentsAndTrilinearInterpolationIncludingLastCell) {
    auto* map = load(writeMap());
    for (const Vector& position :
         {Vector(-0.02, -0.03, -0.05), Vector(0.01, -0.01, 0.05), Vector(-0.011, 0.007, -0.013),
          Vector(0.039, 0.029, 0.149)}) {
        expectField(map, position);
    }
}

TEST_F(FM3DDynamicTest, PositionsImmediatelyInsideUpperFacesStillInterpolate) {
    auto* map = load(writeMap());
    expectField(map, Vector(std::nextafter(0.04, 0.0), 0, 0));
    expectField(map, Vector(0, std::nextafter(0.03, 0.0), 0));
    expectField(map, Vector(0, 0, std::nextafter(0.15, 0.0)));
}

TEST_F(FM3DDynamicTest, OutsideAndUpperFacesLeaveFieldsUnchanged) {
    auto* map        = load(writeMap());
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const Vector& position :
         {Vector(-0.021, 0, 0), Vector(0.04, 0, 0), Vector(0, -0.031, 0), Vector(0, 0.03, 0),
          Vector(0, 0, -0.051), Vector(0, 0, 0.15), Vector(nan, 0, 0),
          Vector(0, 0, std::numeric_limits<double>::infinity())}) {
        Vector electric(7.0), magnetic(-0.25);
        EXPECT_FALSE(map->isInside(position));
        EXPECT_TRUE(map->getFieldstrength(position, electric, magnetic));
        for (unsigned component = 0; component < 3; ++component) {
            EXPECT_DOUBLE_EQ(electric(component), 7.0);
            EXPECT_DOUBLE_EQ(magnetic(component), -0.25);
        }
    }
}

TEST_F(FM3DDynamicTest, NormalizesToInterpolatedOnaxisPeakAndIncludesFinalPlane) {
    for (const std::string flag : {"", "TRUE", "true"}) {
        auto* map = load(writeMap("normalized" + flag + ".map", flag));
        expectField(map, Vector(0.039, 0.029, 0.149), 3.5);
        std::vector<std::pair<double, double>> profile;
        map->getOnaxisEz(profile);
        ASSERT_EQ(profile.size(), 5u);
        for (unsigned index = 0; index < profile.size(); ++index) {
            EXPECT_NEAR(profile[index].first, 0.05 * index, 1e-15);
            EXPECT_NEAR(profile[index].second, (1.5 + 0.5 * index) / 3.5, 1e-14);
        }
    }
}

TEST_F(FM3DDynamicTest, OnaxisSamplingAlsoWorksAtGridNodesAndTransverseUpperFace) {
    for (const auto bounds : {std::pair{-0.03, 0.03}, std::pair{-0.06, 0.0}}) {
        auto* map = load(writeMap(
                std::to_string(bounds.first) + ".map", "TRUE", analyticField, bounds.first,
                bounds.second));
        std::vector<std::pair<double, double>> profile;
        map->getOnaxisEz(profile);
        ASSERT_EQ(profile.size(), 5u);
        EXPECT_NEAR(profile.back().second, 1.0, 1e-14);
    }
}

TEST_F(FM3DDynamicTest, NormalizationUsesAbsolutePeak) {
    auto* map = load(writeMap("negative.map", "TRUE", [](double, double, double z) {
        return Samples{2, 0, -4 + 10 * z, 3, 0, 0};
    }));
    Vector electric(0.0), magnetic(0.0);
    ASSERT_FALSE(map->getFieldstrength(Vector(0, 0, -0.05), electric, magnetic));
    EXPECT_NEAR(electric(2), -1e6, 1e-9);
    EXPECT_NEAR(electric(0), 2e6 / 4.5, 1e-9);
    EXPECT_NEAR(magnetic(0), 3 * Physics::mu_0 / 4.5, 1e-18);
}

TEST_F(FM3DDynamicTest, ZeroLongitudinalFieldRequiresDisabledNormalization) {
    const auto transverse = [](double, double, double) {
        return Samples{2, 3, 0, 4, 5, 6};
    };
    const auto invalid = writeMap("zero.map", "TRUE", transverse);
    Fieldmap::getFieldmap(invalid);
    EXPECT_THROW(Fieldmap::readMap(invalid), GeneralOpalException);
    auto* map = load(writeMap("transverse.map", "FALSE", transverse));
    Vector electric(0.0), magnetic(0.0);
    ASSERT_FALSE(map->getFieldstrength(Vector(0.0), electric, magnetic));
    EXPECT_DOUBLE_EQ(electric(2), 0.0);
    EXPECT_NEAR(electric(0), 2e6, 1e-9);
    EXPECT_NEAR(magnetic(2), 6 * Physics::mu_0, 1e-18);
}

TEST_F(FM3DDynamicTest, MapsAwayFromAxisWorkWithoutNormalization) {
    const auto invalid = writeMap("offaxis.map", "TRUE", analyticField, 0.01, 0.04);
    Fieldmap::getFieldmap(invalid);
    EXPECT_THROW(Fieldmap::readMap(invalid), GeneralOpalException);
    auto* map = load(writeMap("offaxis_false.map", "FALSE", analyticField, 0.01, 0.04));
    expectField(map, Vector(0.02, 0.0, 0.0));
    std::vector<std::pair<double, double>> profile;
    EXPECT_THROW(map->getOnaxisEz(profile), GeneralOpalException);
}

TEST_F(FM3DDynamicTest, RejectsInvalidHeadersBeforeAllocating) {
    const std::vector<std::string> headers{
            "3DDynamic MAYBE\n125\n-1 1 1\n-1 1 1\n0 1 1\n",
            "3DDynamic FALSE extra\n125\n-1 1 1\n-1 1 1\n0 1 1\n",
            "3DDynamic\n-125\n-1 1 1\n-1 1 1\n0 1 1\n",
            "3DDynamic\n1e309\n-1 1 1\n-1 1 1\n0 1 1\n",
            "3DDynamic\n125\n1 1 1\n-1 1 1\n0 1 1\n",
            "3DDynamic\n125\n-1 1 0\n-1 1 1\n0 1 1\n",
            "3DDynamic\n125\n-1 1 -1\n-1 1 1\n0 1 1\n",
            "3DDynamic\n125\n-1 1 1.5\n-1 1 1\n0 1 1\n",
            "3DDynamic\n125\n-1 1 2147483647\n-1 1 1\n0 1 1\n",
            "3DDynamic\n125\n-1 1 2147483646\n-1 1 2147483646\n0 1 2147483646\n",
            "3DDynamic\n125\n-1 1 1\n1 -1 1\n0 1 1\n",
            "3DDynamic\n125\n-1 1 1\n-1 1 1\n0 1 0\n"};
    for (std::size_t index = 0; index < headers.size(); ++index) {
        SCOPED_TRACE(index);
        const auto path = write(std::to_string(index) + ".map", headers[index]);
        EXPECT_THROW(Fieldmap::getFieldmap(path), GeneralOpalException);
    }
}

TEST_F(FM3DDynamicTest, RejectsBadDataAndAllowsRetryAfterRepair) {
    const std::string header = "3DDynamic FALSE\n125\n-1 1 1\n-1 1 1\n0 1 1\n";
    std::string valid;
    for (int point = 0; point < 8; ++point)
        valid += "1 2 3 4 5 6\n";
    const std::vector<std::string> invalid{
            "",
            valid + "1 2 3 4 5 6\n",
            "1 2 3 4 5\n" + valid,
            "1 2 3 4 5 6 7\n" + valid,
            "nan 2 3 4 5 6\n" + valid,
            "1 2 inf 4 5 6\n" + valid,
            "1e308 2 3 4 5 6\n" + valid.substr(12)};
    for (std::size_t index = 0; index < invalid.size(); ++index) {
        SCOPED_TRACE(index);
        const auto name = std::to_string(index) + ".map";
        const auto path = write(name, header + invalid[index]);
        auto* map       = Fieldmap::getFieldmap(path);
        EXPECT_THROW(Fieldmap::readMap(path), GeneralOpalException);
        write(name, header + valid);
        ASSERT_NO_THROW(Fieldmap::readMap(path));
        Vector electric(0.0), magnetic(0.0);
        ASSERT_FALSE(map->getFieldstrength(Vector(0, 0, 0.005), electric, magnetic));
        EXPECT_NEAR(electric(2), 3e6, 1e-9);
    }
}

TEST_F(FM3DDynamicTest, ParticleKernelAccumulatesWithIndependentRFScalesAndBounds) {
    const auto path = writeMap();
    auto* map       = dynamic_cast<FM3DDynamic*>(load(path));
    ASSERT_NE(map, nullptr);
    Beam beam;
    PartBunch_t bunch(
            {1.}, {1.}, {&beam}, {0}, 1., "LF2", opalx::spacecharge::CartesianDomainConfig3D{});
    auto particles = bunch.getParticleContainer();
    const std::vector<Vector> positions{Vector(0.005, 0.005, 0.025), Vector(0.039, 0.029, 0.149),
                                        Vector(-0.02, -0.03, -0.05), Vector(0, 0, -0.025),
                                        Vector(0.04, 0, 0.05),       Vector(0, 0, 0.15)};
    particles->createParticles(positions.size());
    auto hostPositions = Kokkos::create_mirror_view(particles->R.getView());
    auto hostElectric  = Kokkos::create_mirror_view(particles->E.getView());
    auto hostMagnetic  = Kokkos::create_mirror_view(particles->B.getView());
    for (std::size_t index = 0; index < positions.size(); ++index)
        hostPositions(index) = positions[index];
    Kokkos::deep_copy(particles->R.getView(), hostPositions);
    for (const bool rf : {false, true}) {
        Kokkos::deep_copy(particles->E.getView(), Vector(7.0));
        Kokkos::deep_copy(particles->B.getView(), Vector(-0.25));
        if (rf)
            map->applyRFField(particles, 2.5, -1.75, -0.025, 0.149);
        else
            map->applyField(particles, 1.0, 1.0);
        Kokkos::deep_copy(hostElectric, particles->E.getView());
        Kokkos::deep_copy(hostMagnetic, particles->B.getView());
        for (std::size_t index = 0; index < positions.size(); ++index) {
            Vector electric(0.0), magnetic(0.0);
            if (!rf || (positions[index](2) >= -0.025 && positions[index](2) < 0.149))
                map->getFieldstrength(positions[index], electric, magnetic);
            for (unsigned component = 0; component < 3; ++component) {
                EXPECT_NEAR(
                        hostElectric(index)(component),
                        7.0 + (rf ? 2.5 : 1.0) * electric(component), 2e-8);
                EXPECT_NEAR(
                        hostMagnetic(index)(component),
                        -0.25 + (rf ? -1.75 : 1.0) * magnetic(component), 1e-15);
            }
        }
    }
}

TEST_F(FM3DDynamicTest, RFCavityParticleAndPointEvaluationAgreeAtMidstepTime) {
    const auto path = writeMap();
    Beam beam;
    PartBunch_t bunch(
            {1.}, {1.}, {&beam}, {0}, 1., "LF2", opalx::spacecharge::CartesianDomainConfig3D{});
    bunch.setT(1.3e-9);
    bunch.setdT(0.4e-9);
    auto particles = bunch.getParticleContainer();
    particles->createParticles(1);
    const Vector position(0.015, -0.005, 0.09);
    Kokkos::deep_copy(particles->R.getView(), position);
    RFCavityRep cavity("3D cavity");
    cavity.setFieldMapFN(path);
    cavity.setFrequencym(125 * Physics::two_pi * Units::MHz2Hz);
    cavity.setAmplitudem(2.0);
    cavity.setAmplitudeError(0.3);
    cavity.setPhaseError(-0.15);
    cavity.initialise(&bunch);
    cavity.goOnline(0.0);
    for (const double phase : {0.0, 0.7, Physics::pi / 2, Physics::pi}) {
        cavity.setPhasem(phase);
        Kokkos::deep_copy(particles->E.getView(), Vector(7.0));
        Kokkos::deep_copy(particles->B.getView(), Vector(-0.25));
        cavity.apply(particles);
        const auto electric =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles->E.getView());
        const auto magnetic =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles->B.getView());
        Vector expectedElectric(7.0), expectedMagnetic(-0.25);
        cavity.apply(
                position, Vector(0.0), bunch.getT() + 0.5 * bunch.getdT(), expectedElectric,
                expectedMagnetic);
        for (unsigned component = 0; component < 3; ++component) {
            EXPECT_NEAR(electric(0)(component), expectedElectric(component), 2e-8);
            EXPECT_NEAR(magnetic(0)(component), expectedMagnetic(component), 1e-15);
        }
    }
    cavity.goOffline();
}
