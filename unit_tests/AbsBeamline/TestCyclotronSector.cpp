#include <gtest/gtest.h>
#include <unistd.h>
#include <filesystem>
#include <fstream>
#include "AbsBeamline/CyclotronSector.h"
#include "Algorithms/DirectedTurnCounter.h"
#include "PartBunch/BunchStateHandler.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Utilities/OpalException.h"

class CyclotronTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }
    static void TearDownTestSuite() { ippl::finalize(); }
    std::string filename;
    void SetUp() override {
        filename = (std::filesystem::temp_directory_path()
                    / ("opalx-cyclotron-test-" + std::to_string(getpid()) + ".map"))
                           .string();
        // Synthetic PSI layout: 5 radii, 8 angular points, one 45-degree sector.
        std::ofstream out(filename);
        out << "1000 100 0 5.625\n";
        for (int i = 0; i < 13; ++i)
            out << "header ";
        out << "5 8 a b c d e 2 a b c d 0.100000000E+01-0.100000000E+01 a b c d e f LREC= a b c d "
               "e\n";
        for (int i = 0; i < 5; ++i) {
            if (i) out << "a b c d e f\n";
            for (int channel = 0; channel < 4; ++channel)
                for (int j = 0; j < 8; ++j)
                    out << (channel == 0 ? 10 * (1 + 0.1 * i) : 0) << ' ';
            out << '\n';
        }
    }
    void TearDown() override { std::filesystem::remove(filename); }
};
TEST_F(CyclotronTest, MapUnitsDerivativesSeamAndSharedStorage) {
    auto map = CyclotronSectorFieldMap::read(filename);
    EXPECT_EQ(map.get(), CyclotronSectorFieldMap::read(filename).get());
    EXPECT_NEAR(map->rmin, 1, 1e-15);
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j <= 8; ++j) {
            EXPECT_NEAR(map->host(i, j, 0), 1 + 0.1 * i, 1e-14);
            EXPECT_NEAR(map->host(i, j, 1), 1, 1e-13);
            EXPECT_NEAR(map->host(i, j, 2), 0, 1e-13);
        }
    EXPECT_NEAR(CyclotronSectorFieldMap::interpolate(map->host, 0, 4, 8, 5, 8), 1.4, 1e-14);
}
TEST_F(CyclotronTest, WedgeBoundsAndOffPlaneField) {
    CyclotronSector sector("test");
    sector.configure(CyclotronSectorFieldMap::read(filename), 8, -0.05, 0.05, 1, {});
    const double angle = 0.2, radius = 1.25;
    Vector_t<double, 3> r(radius * std::cos(angle) - 1.2, 0.01, radius * std::sin(angle));
    Vector_t<double, 3> e(0), b(0), p(0);
    ASSERT_TRUE(sector.isInside(r));
    sector.apply(r, p, 0, e, b);
    EXPECT_NEAR(b[0], 0.01 * std::cos(angle), 1e-14);
    EXPECT_NEAR(b[1], 1.25, 1e-14);
    EXPECT_NEAR(b[2], 0.01 * std::sin(angle), 1e-14);
    EXPECT_FALSE(sector.isInside(Vector_t<double, 3>(0, 0, -1e-8)));
    EXPECT_FALSE(sector.isInside(Vector_t<double, 3>(0, 0.051, 0)));
    EXPECT_TRUE(sector.isInside(Vector_t<double, 3>(0, 0, 0)));
}
TEST_F(CyclotronTest, RejectsTruncatedMap) {
    {
        std::ofstream out(filename);
        out << "1000 100 0 5.625";
    }
    EXPECT_THROW(CyclotronSectorFieldMap::read(filename), OpalException);
}
TEST_F(CyclotronTest, EightRotatedSectorsOwnEverySeamExactlyOnce) {
    CyclotronSector sector("sector");
    sector.configure(CyclotronSectorFieldMap::read(filename), 8, -0.05, 0.05, 1, {});
    for (int seam = 0; seam < 8; ++seam)
        for (double offset : {-1e-8, 0.0, 1e-8}) {
            const double theta = seam * std::acos(-1.0) / 4 + offset;
            int owners         = 0;
            double total       = 0;
            for (int k = 0; k < 8; ++k) {
                const double a = theta - k * std::acos(-1.0) / 4;
                Vector_t<double, 3> r(1.25 * std::cos(a) - 1.2, 0, 1.25 * std::sin(a));
                owners += sector.isInsideBody(r);
                Vector_t<double, 3> p(0), e(0), b(0);
                sector.apply(r, p, 0, e, b);
                total += b[1];
            }
            EXPECT_EQ(owners, 1) << "seam=" << seam << " offset=" << offset;
            EXPECT_NEAR(total, 1.25, 1e-13);
        }
}
TEST_F(CyclotronTest, DeviceParticleFieldMatchesHostIncludingTrim) {
    ippl::NDIndex<3> domain;
    for (int d = 0; d < 3; ++d)
        domain[d] = ippl::Index(8);
    ippl::UniformCartesian<double, 3> mesh(
            domain, Vector_t<double, 3>(0.25), Vector_t<double, 3>(-1));
    std::array<bool, 3> decomp{true, true, true};
    ippl::FieldLayout<3> layout(MPI_COMM_WORLD, domain, decomp, false);
    auto pc = std::make_shared<ParticleContainer_t>(mesh, layout);
    pc->setBunchStateHandler(std::make_shared<BunchStateHandler>());
    pc->allocateParticles(2);
    pc->createParticles(2);
    auto r = Kokkos::create_mirror_view(pc->R.getView());
    r(0)   = Vector_t<double, 3>(0, 0.01, 0.1);
    r(1)   = Vector_t<double, 3>(0, 0, -0.1);
    Kokkos::deep_copy(pc->R.getView(), r);
    Kokkos::deep_copy(pc->B.getView(), Vector_t<double, 3>(0));
    CyclotronSector sector("test");
    sector.configure(
            CyclotronSectorFieldMap::read(filename), 8, -0.05, 0.05, 1, {{1.1, 1.3, 0.0014, 600}});
    sector.apply(pc);
    auto b = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->B.getView());
    for (int i = 0; i < 2; ++i) {
        Vector_t<double, 3> expected(0), e(0), p(0);
        sector.apply(r(i), p, 0, e, expected);
        for (int d = 0; d < 3; ++d)
            EXPECT_NEAR(b(i)[d], expected[d], 1e-14);
    }
}
TEST_F(CyclotronTest, MirroredCoilHasTailsAndZeroStrengthIsIdentity) {
    CyclotronTrimCoil coil{4.35, 4.47, 0.0014, 600};
    double br = 0, by = 0;
    coil.add(2.1314, 0, br, by);
    EXPECT_DOUBLE_EQ(br, 0);
    EXPECT_NE(by, 0);  // Shape radii are not hard cutoffs.
    double left = 0, right = 0, unused = 0;
    coil.add(4.38, 0, unused, left);
    coil.add(4.44, 0, unused, right);
    EXPECT_NEAR(left, right, 1e-15);
    coil.bmax = 0;
    br        = 2;
    by        = 3;
    coil.add(4.41, 0.01, br, by);
    EXPECT_DOUBLE_EQ(br, 2);
    EXPECT_DOUBLE_EQ(by, 3);
}
TEST_F(CyclotronTest, CountsDirectedReturnsWithoutInitialOrOppositeCrossing) {
    Vector_t<double, 3> origin(1, 0, 0), momentum(0, 0, 1);
    DirectedTurnCounter counter(origin, momentum);
    EXPECT_FALSE(counter.update(origin, momentum));
    EXPECT_FALSE(counter.update(Vector_t<double, 3>(0, 0, 1), Vector_t<double, 3>(-1, 0, 0)));
    EXPECT_FALSE(counter.update(Vector_t<double, 3>(-1, 0, -1e-10), Vector_t<double, 3>(0, 0, -1)));
    EXPECT_FALSE(counter.update(Vector_t<double, 3>(0, 0, -1), Vector_t<double, 3>(1, 0, 0)));
    EXPECT_TRUE(counter.update(origin, momentum));
    EXPECT_FALSE(counter.update(origin, momentum));
    EXPECT_EQ(counter.count(), 1);
}
