#include <gtest/gtest.h>
#include <unistd.h>
#include <filesystem>
#include <fstream>
#include "AbsBeamline/CyclotronSector.h"
#include "Algorithms/DeviceExternalField.h"
#include "Algorithms/DirectedTurnCounter.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/RFCavityRep.h"
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

TEST_F(CyclotronTest, SpatialWedgeCandidatesAreConservativeAcrossFramesAndSeams) {
    CyclotronSector sector("candidate");
    sector.configure(CyclotronSectorFieldMap::read(filename), 8, -0.05, 0.05, 1, {});
    const CoordinateSystemTrafo frame(
            Vector_t<double, 3>(2, -3, 4),
            Quaternion(std::cos(.3), Vector_t<double, 3>(0, std::sin(.3), 0)));
    const auto& geometry = sector.getGeometry();
    for (int i = 0; i <= 20; ++i) {
        const double angle = (std::acos(-1.0) / 4) * i / 20;
        const auto point   = frame.transformFrom(
                Vector_t<double, 3>(1.25 * std::cos(angle) - 1.2, 0, 1.25 * std::sin(angle)));
        EXPECT_TRUE(device_external::cyclotronSectorMayIntersect(geometry, frame, point, point))
                << i;
    }
    const auto behind = frame.transformFrom(Vector_t<double, 3>(0, 0, -.1));
    EXPECT_FALSE(device_external::cyclotronSectorMayIntersect(geometry, frame, behind, behind));
    const auto beyond = frame.transformFrom(Vector_t<double, 3>(-.9, 0, 1));
    EXPECT_FALSE(device_external::cyclotronSectorMayIntersect(geometry, frame, beyond, beyond));
    // A finite bunch straddling the entrance must not be rejected by its centre.
    EXPECT_TRUE(
            device_external::cyclotronSectorMayIntersect(
                    geometry, CoordinateSystemTrafo(), Vector_t<double, 3>(0, -.01, -.1),
                    Vector_t<double, 3>(.1, .01, .1)));
    EXPECT_TRUE(
            device_external::cyclotronSectorMayIntersect(
                    geometry, CoordinateSystemTrafo(), Vector_t<double, 3>(0, -100, 0),
                    Vector_t<double, 3>(100, 100, 100)));
    auto invalid = Vector_t<double, 3>(std::numeric_limits<double>::quiet_NaN());
    EXPECT_TRUE(device_external::cyclotronSectorMayIntersect(geometry, frame, invalid, invalid));
    // SYM=2 permits the complete upper half-plane, including negative local X.
    const auto half = Geometry::makeSBend(1.2 * std::acos(-1.0), 1 / 1.2);
    const Vector_t<double, 3> left(-2.4, 0, .1);
    EXPECT_TRUE(
            device_external::cyclotronSectorMayIntersect(
                    half, CoordinateSystemTrafo(), left, left));
}

TEST_F(CyclotronTest, SpatialSelectionRestoresOmittedDeviceFieldAndPreservesOtherCandidates) {
    auto sector = std::make_shared<CyclotronSector>("active");
    sector->configure(CyclotronSectorFieldMap::read(filename), 8, -.05, .05, 1, {});
    sector->setCSTrafoGlobal2Local(CoordinateSystemTrafo(
            Vector_t<double, 3>(.3, -.2, .4),
            Quaternion(std::cos(.15), Vector_t<double, 3>(0, std::sin(.15), 0))));
    sector->setMisalignment(CoordinateSystemTrafo(Vector_t<double, 3>(.01, 0, .02), Quaternion()));
    std::shared_ptr<ElementBase> staleSector(sector->clone());
    staleSector->setCSTrafoGlobal2Local(
            CoordinateSystemTrafo(Vector_t<double, 3>(0, 0, 10), Quaternion()));
    auto drift = std::make_shared<DriftRep>("preserved");
    const CoordinateSystemTrafo toLab(
            Vector_t<double, 3>(-.1, .2, -.3),
            Quaternion(std::cos(.1), Vector_t<double, 3>(0, std::sin(.1), 0)));
    const auto frame = sector->getMisalignment() * (sector->getCSTrafoGlobal2Local() * toLab);
    const Vector_t<double, 3> local(1.25 * std::cos(.2) - 1.2, .01, 1.25 * std::sin(.2));
    const auto point = frame.transformFrom(local);
    // Model a stale periodic index selecting the neighbouring sector instead.
    auto selected = device_external::selectCyclotronCandidates(
            {sector, staleSector}, toLab, point, point, {staleSector, drift});
    ASSERT_EQ(selected.size(), 2u);
    EXPECT_EQ(selected.count(sector), 1u);
    EXPECT_EQ(selected.count(staleSector), 0u);
    EXPECT_EQ(selected.count(drift), 1u);

    ippl::NDIndex<3> domain;
    for (unsigned d = 0; d < 3; ++d)
        domain[d] = ippl::Index(8);
    Mesh_t<3> mesh(domain, Vector_t<double, 3>(1), Vector_t<double, 3>(-4));
    std::array<bool, 3> decomp{true, true, true};
    FieldLayout_t<3> layout(ippl::Comm->getCommunicator(), domain, decomp, true);
    auto pc = std::make_shared<ParticleContainer_t>(mesh, layout);
    pc->setBunchStateHandler(std::make_shared<BunchStateHandler>());
    pc->createParticles(ippl::Comm->rank() == 0 ? 1 : 0);
    Kokkos::deep_copy(pc->R.getView(), point);
    const Vector_t<double, 3> momentum(.1, .2, .3);
    Kokkos::deep_copy(pc->P.getView(), momentum);
    Kokkos::deep_copy(pc->B.getView(), Vector_t<double, 3>(0));
    Kokkos::deep_copy(pc->E.getView(), Vector_t<double, 3>(0));
    pc->transformBunch(frame);
    // Exercise the existing Kokkos kernel, not a host replacement field.
    sector->apply(pc);
    pc->transformBunch(frame.inverted());
    const auto b = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->B.getView());
    const auto r = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
    const auto p = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
    Vector_t<double, 3> expected(0), e(0);
    sector->apply(local, Vector_t<double, 3>(0), 0, e, expected);
    expected = frame.rotateFrom(expected);
    for (std::size_t i = 0; i < pc->getLocalNum(); ++i)
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_NEAR(b(i)(d), expected(d), 1e-14);
            EXPECT_NEAR(r(i)(d), point(d), 1e-14);
            EXPECT_NEAR(p(i)(d), momentum(d), 1e-14);
        }
}

namespace {
    // NVCC extended lambdas require an accessible enclosing function; GoogleTest's
    // generated TestBody is private. Keep the real device kernels in a free helper.
    void checkDeviceBunchTransport(const std::string& filename) {
        const auto rfFile = filename + ".rf";
        {
            std::ofstream out(rfFile);
            out << "2\n0 1 0\n1 1 0\n";
        }
        CyclotronSector sector("zero");
        sector.configure(CyclotronSectorFieldMap::read(filename), 8, -.05, .05, 0, {});
        RFCavityRep cavity("gap");
        const double mass = 938272088.16, dt = 1e-11, time = 2e-9;
        cavity.configureCyclotronGap(
                rfFile, .001, 1, CyclotronRFKick{847000, Physics::two_pi * 50.65e6, 0, 0, .999});
        device_external::CyclotronLattice lattice;
        lattice.sectors = decltype(lattice.sectors)("sectors", 1);
        lattice.gaps    = decltype(lattice.gaps)("gaps", 1);
        auto sectors    = Kokkos::create_mirror_view(lattice.sectors);
        auto gaps       = Kokkos::create_mirror_view(lattice.gaps);
        const CoordinateSystemTrafo identity;
        const device_external::Rigid frame{identity.getOrigin(), identity.getRotationMatrix()};
        const Vector_t<double, 3> start(.02, .001, .25), momentum(.01, .002, .4);
        const CoordinateSystemTrafo plane(Vector_t<double, 3>(0, 0, .2502), Quaternion());
        sectors(0) = {frame, sector.deviceField()};
        gaps(0)    = {{plane.getOrigin(), plane.getRotationMatrix()}, cavity.deviceGap()};
        Kokkos::deep_copy(lattice.sectors, sectors);
        Kokkos::deep_copy(lattice.gaps, gaps);
        Kokkos::View<Vector_t<double, 3>*> result("result", 2);
        Kokkos::View<int*> accepted("accepted", 2);
        const size_t count = ippl::Comm->rank() == 0 ? 2 : 0;
        Kokkos::parallel_for(
                "bunch event test", count, KOKKOS_LAMBDA(size_t i) {
                    auto r = start, p = momentum;
                    r(2) -= i * .0001;
                    accepted(i) = lattice.advance(r, p, time, dt, mass);
                    result(i)   = p;
                });
        auto host          = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), result);
        auto ok            = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), accepted);
        const double gamma = std::sqrt(1 + dot(momentum, momentum));
        for (size_t i = 0; i < count; ++i) {
            ASSERT_EQ(ok(i), 1);
            const double arrival  = (.0002 + i * .0001) * gamma / (Physics::c * momentum(2));
            const double expected = 847000 * std::cos(Physics::two_pi * 50.65e6 * (time + arrival));
            EXPECT_NEAR((std::sqrt(1 + dot(host(i), host(i))) - gamma) * mass, expected, 5e-7);
            EXPECT_EQ(host(i)(1), momentum(1));
            EXPECT_EQ(host(i)(0), momentum(0));
        }
        // Frozen uniform self-field split equals a full momentum impulse without RF.
        lattice.gaps = decltype(lattice.gaps)("no gaps", 0);
        const Vector_t<double, 3> electric(0, 1000, 0);
        Kokkos::parallel_for(
                "self field test", count, KOKKOS_LAMBDA(size_t i) {
                    auto r = start, p = momentum;
                    accepted(i) = lattice.advanceWithSelfFields(
                            r, p, electric, Vector_t<double, 3>(0), time, dt, mass);
                    result(i) = p;
                });
        Kokkos::deep_copy(host, result);
        Kokkos::deep_copy(ok, accepted);
        for (size_t i = 0; i < count; ++i) {
            ASSERT_EQ(ok(i), 1);
            EXPECT_NEAR(host(i)(1), momentum(1) + dt * Physics::c / mass * 1000, 1e-16);
        }
        // Failures are transactional and distinguish support, input and RF causes.
        Kokkos::parallel_for(
                "support diagnostics", count, KOKKOS_LAMBDA(size_t i) {
                    auto r = start, p = momentum;
                    r(1) = 1;
                    rf_event::Diagnostic diagnostic;
                    const bool ok = lattice.advance(r, p, time, dt, mass, &diagnostic);
                    accepted(i)   = ok;
                    result(i)     = Vector_t<double, 3>(int(diagnostic.failure), r(1), p(1));
                });
        Kokkos::deep_copy(host, result);
        Kokkos::deep_copy(ok, accepted);
        for (size_t i = 0; i < count; ++i) {
            EXPECT_EQ(ok(i), 0);
            EXPECT_EQ(host(i)(0), 1);
            EXPECT_EQ(host(i)(1), 1);
            EXPECT_EQ(host(i)(2), momentum(1));
        }
        lattice.gaps               = decltype(lattice.gaps)("invalid kick gap", 1);
        gaps(0).field.kick.voltage = -1e12;
        Kokkos::deep_copy(lattice.gaps, gaps);
        Kokkos::parallel_for(
                "kick diagnostics", count, KOKKOS_LAMBDA(size_t i) {
                    auto r = start, p = momentum;
                    rf_event::Diagnostic diagnostic;
                    accepted(i) = lattice.advance(r, p, time, dt, mass, &diagnostic);
                    result(i)   = Vector_t<double, 3>(int(diagnostic.failure), diagnostic.gap, 0);
                });
        Kokkos::deep_copy(host, result);
        Kokkos::deep_copy(ok, accepted);
        for (size_t i = 0; i < count; ++i) {
            EXPECT_EQ(ok(i), 0);
            EXPECT_EQ(host(i)(0), 2);
            EXPECT_EQ(host(i)(1), 0);
        }
        std::filesystem::remove(rfFile);
    }
}  // namespace

TEST_F(CyclotronTest, DeviceBunchTransportUsesIndividualArrivalTimesAndThreeDimensionalEnergy) {
    checkDeviceBunchTransport(filename);
}
