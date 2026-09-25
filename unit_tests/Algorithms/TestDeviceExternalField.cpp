#include <unistd.h>
#include <array>
#include <filesystem>
#include <fstream>
#include "AbsBeamline/CyclotronSector.h"
#include "AbsBeamline/VariableRFCavity.h"
#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/DeviceExternalField.h"
#include "Algorithms/ExternalFieldRayTracker.h"
#include "Algorithms/PartData.h"
#include "BeamlineCore/ConstantFocusingRep.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/MarkerRep.h"
#include "BeamlineCore/MultipoleRep.h"
#include "BeamlineCore/RBendRep.h"
#include "BeamlineCore/SBendRep.h"
#include "Beamlines/FlaggedBeamline.h"
#include "Elements/OpalBeamline.h"
#include "PartBunch/PartBunch.h"
#include "Structure/Beam.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

using Vector = device_external::Vector;

namespace {
    struct MagneticFieldKernel {
        device_external::Lattice device;
        Kokkos::View<Vector*> points;
        Kokkos::View<Vector*> fields;

        KOKKOS_INLINE_FUNCTION void operator()(int i) const {
            fields(i) = device.magnetic(points(i));
        }
    };

    struct SupportKernel {
        Kokkos::View<const device_external::Element*> elements;
        Kokkos::View<Vector*> points;
        Kokkos::View<int*> flags;

        KOKKOS_INLINE_FUNCTION void operator()(int i) const {
            flags(i) = elements(0).contains(points(i));
        }
    };
}  // namespace

class DeviceExternalFieldTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        gmsg = new Inform(nullptr, -1);
    }
    static void TearDownTestSuite() {
        delete gmsg;
        gmsg = nullptr;
        ippl::finalize();
    }
};

TEST_F(DeviceExternalFieldTest, OrderedCandidatesPreserveOccurrencesAndRejectUnknowns) {
    Beam beam;
    PartBunch_t bunch(
            {1.}, {1.}, {&beam}, {0}, 1., "LF2", opalx::spacecharge::CartesianDomainConfig3D{});
    FlaggedBeamline line;
    DefaultVisitor visitor(line, false, false);
    OpalBeamline lattice;
    for (const auto* name : {"shared", "other", "shared"}) {
        DriftRep drift(name);
        drift.getGeometry().setElementLength(1.);
        drift.setCSTrafoGlobal2Local(CoordinateSystemTrafo());
        drift.fixPosition();
        lattice.visit(drift, visitor, bunch);
    }
    const auto declared = lattice.getElementByType(ElementType::ANY);
    const auto first = declared.front(), last = declared.back();
    ASSERT_NE(first, last);
    ASSERT_EQ(first->getName(), last->getName());
    lattice.prepareSections();
    std::set<std::shared_ptr<ElementBase>> selected;
    selected.insert(last);
    selected.insert(first);
    const std::vector<std::shared_ptr<ElementBase>> expected{first, last};
    EXPECT_EQ(device_external::orderedCandidates(lattice, selected), expected);
    EXPECT_TRUE(device_external::orderedCandidates(lattice, {}).empty());

    // A merge can repeat an existing pointer. Retain set selection semantics,
    // without accidentally merging distinct occurrences that share a name.
    OpalBeamline duplicate = lattice;
    lattice.merge(duplicate);
    lattice.prepareSections();
    EXPECT_EQ(device_external::orderedCandidates(lattice, selected), expected);
    EXPECT_EQ(device_external::orderedCandidates(lattice, lattice.getElements()).size(), 3u);
    selected.insert(std::make_shared<DriftRep>("shared"));
    EXPECT_THROW(device_external::orderedCandidates(lattice, selected), OpalException);
    EXPECT_THROW(device_external::orderedCandidates(lattice, {nullptr}), OpalException);
}

TEST_F(DeviceExternalFieldTest, PreparedOrderMakesDeviceFrameAndFieldWorkIndependentOfAddresses) {
    Beam beam;
    PartBunch_t bunch(
            {1.}, {1.}, {&beam}, {0}, 1., "LF2", opalx::spacecharge::CartesianDomainConfig3D{});
    const auto particles = bunch.getParticleContainer();
    particles->createParticles(3);
    FlaggedBeamline line;
    DefaultVisitor visitor(line, false, false);
    std::array<std::array<double, 36>, 2> snapshots{};
    const std::array<std::string, 3> names{"first", "second", "third"};
    for (unsigned arrangement = 0; arrangement < 2; ++arrangement) {
        OpalBeamline lattice;
        for (unsigned i = 0; i < 3; ++i) {
            MultipoleRep magnet("unassigned");
            magnet.getGeometry().setElementLength(3.);
            lattice.visit(magnet, visitor, bunch);
        }
        // Assign physical roles AFTER allocating the occurrences. This forces
        // opposite pointer-to-physics order, independent of allocator behavior.
        const auto allocated = lattice.getElements();
        unsigned index       = 0;
        for (const auto& occurrence : allocated) {
            const unsigned role = arrangement == 0 ? index : 2 - index;
            auto& magnet        = dynamic_cast<Multipole&>(*occurrence);
            magnet.setName(names[role]);
            magnet.setNormalComponent(0, 0.4 * (role + 1));
            magnet.setNormalComponent(1, 0.2 * (role + 1));
            magnet.setSkewComponent(1, 0.07 * (role + 1));
            magnet.setAperture(ApertureType::RECTANGULAR, {2., 2.});
            const double angle = 0.11 + 0.19 * role;
            magnet.setCSTrafoGlobal2Local(CoordinateSystemTrafo(
                    Vector(0.03 * (role + 1), -0.05 * (role + 1), -0.7),
                    Quaternion(std::cos(angle / 2), 0., std::sin(angle / 2), 0.)));
            magnet.fixPosition();
            ++index;
        }
        lattice.prepareSections();
        ASSERT_EQ((*allocated.begin())->getName(), names[arrangement == 0 ? 0 : 2]);
        std::set<std::shared_ptr<ElementBase>> selected;
        for (auto it = allocated.rbegin(); it != allocated.rend(); ++it)
            selected.insert(*it);
        const auto ordered = device_external::orderedCandidates(lattice, selected);
        ASSERT_EQ(ordered.size(), names.size());
        for (unsigned i = 0; i < names.size(); ++i)
            EXPECT_EQ(ordered[i]->getName(), names[i]);

        auto r = Kokkos::create_mirror(particles->R.getView());
        auto p = Kokkos::create_mirror(particles->P.getView());
        for (unsigned i = 0; i < 3; ++i) {
            r(i) = Vector(0.0123456789 * (i + 1), -0.0234567891 * (i + 1), 0.3123456789);
            p(i) = Vector(0.0034567891 * (i + 1), -0.0045678912 * (i + 1), 0.0734567891);
        }
        Kokkos::deep_copy(particles->R.getView(), r);
        Kokkos::deep_copy(particles->P.getView(), p);
        Kokkos::deep_copy(particles->E.getView(), Vector(0.123456789, -0.234567891, 0.345678912));
        Kokkos::deep_copy(particles->B.getView(), Vector(0.));
        // Exercise the actual device element kernel and the same R/P/E/B
        // round trips used by production, including overlapping field supports.
        for (const auto& element : ordered) {
            const auto transform = lattice.getCSTrafoLab2Local(element);
            particles->transformBunch(transform);
            element->apply(particles);
            particles->transformBunch(transform.inverted());
        }
        const auto outR =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles->R.getView());
        const auto outP =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles->P.getView());
        const auto outE =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles->E.getView());
        const auto outB =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles->B.getView());
        for (unsigned i = 0; i < 3; ++i) {
            EXPECT_GT(dot(outB(i), outB(i)), 0.1);
            for (unsigned d = 0; d < 3; ++d) {
                snapshots[arrangement][12 * i + d]     = outR(i)(d);
                snapshots[arrangement][12 * i + 3 + d] = outP(i)(d);
                snapshots[arrangement][12 * i + 6 + d] = outE(i)(d);
                snapshots[arrangement][12 * i + 9 + d] = outB(i)(d);
            }
        }
    }
    for (unsigned i = 0; i < snapshots[0].size(); ++i) {
        EXPECT_TRUE(std::isfinite(snapshots[0][i]));
        EXPECT_EQ(snapshots[0][i], snapshots[1][i]) << "phase/field component " << i;
    }
}

TEST_F(DeviceExternalFieldTest, ThinMagnetDeviceImpulseEnergyAndSpectators) {
    Beam beam;
    PartBunch_t bunch(
            {1.}, {1.}, {&beam}, {0}, 1., "LF2", opalx::spacecharge::CartesianDomainConfig3D{});
    FlaggedBeamline line;
    DefaultVisitor visitor(line, false, false);
    OpalBeamline lattice;
    MultipoleRep magnet("thin");
    magnet.getGeometry().setElementLength(0.003);
    // Multipole's existing n=0 setter stores half its argument.
    magnet.setNormalComponent(0, 1.6);
    ASSERT_DOUBLE_EQ(magnet.getNormalComponent(0), 0.8);
    magnet.setAperture(ApertureType::RECTANGULAR, {0.02, 0.02});
    magnet.setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector(0, 0, 0.01), Quaternion()));
    magnet.fixPosition();
    lattice.visit(magnet, visitor, bunch);
    lattice.prepareSections();
    ASSERT_TRUE(device_external::Builder::supports(lattice));
    const auto device = device_external::Builder::build(lattice);
    const double mass = 9.382720813e8, dt = 0.1 * std::sqrt(2.) / Physics::c;
    Kokkos::View<Vector*> r("r", 3), p("p", 3), e("e", 3), b("b", 3);
    auto rh = Kokkos::create_mirror_view(r);
    rh(0) = rh(2) = Vector(0);
    rh(1)         = Vector(0.1, 0, 0);
    Kokkos::deep_copy(r, rh);
    Kokkos::deep_copy(p, Vector(0, 0, 1));
    ASSERT_EQ(device.transport(r, p, e, b, 3, {}, dt, mass, 1), 0);
    const auto result    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), p);
    const auto positions = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), r);
    EXPECT_NEAR(result(0)(0), -Physics::c * 0.8 * 0.003 / mass, 2e-13);
    EXPECT_NEAR(dot(result(0), result(0)), 1., 2e-13);
    EXPECT_DOUBLE_EQ(result(1)(0), 0.);  // Exact spatial support, not an infinite slab.
    for (unsigned d = 0; d < 3; ++d) {
        EXPECT_DOUBLE_EQ(result(0)(d), result(2)(d));
        EXPECT_DOUBLE_EQ(positions(0)(d), positions(2)(d));
    }
    const PartData reference(1., mass, 1e6);
    ExternalFieldRayTracker host(lattice, reference);
    ExternalFieldRayTracker::State start;
    start.momentum      = Vector(0, 0, 1);
    const auto expected = host.advance(start, dt);
    for (unsigned d = 0; d < 3; ++d) {
        EXPECT_NEAR(positions(0)(d), expected.position(d), 2e-13);
        EXPECT_NEAR(result(0)(d), expected.momentum(d), 2e-13);
    }
    EXPECT_EQ(device.transport(r, p, e, b, 3, {}, 0, mass, 1), 3);
}

TEST_F(DeviceExternalFieldTest, DeviceFieldsMatchPlacedHostBendsWithFringesAndMultipole) {
    Beam beam;
    PartBunch_t bunch(
            {1.}, {1.}, {&beam}, {0}, 1., "LF2", opalx::spacecharge::CartesianDomainConfig3D{});
    for (int kind = 0; kind < 3; ++kind) {
        FlaggedBeamline line;
        DefaultVisitor visitor(line, false, false);
        OpalBeamline lattice;
        SBendRep sector("sector");
        RBendRep rectangular("rectangular");
        MultipoleRep quad("quad");
        sector.getGeometry() = Geometry::makeSBend(1., 0.5);
        rectangular.getGeometry().setElementLength(1.);
        quad.getGeometry().setElementLength(1.);
        sector.setFieldComponents({0.8, 0.3}, {0.1, -0.2});
        rectangular.setFieldComponents({0.8, 0.3}, {0.1, -0.2});
        sector.setFullGap(0.03);
        rectangular.setFullGap(0.03);
        sector.setFringeIntegral(0.4);
        rectangular.setFringeIntegral(0.4);
        quad.setNormalComponent(0, 0.8);
        quad.setNormalComponent(1, 0.3);
        quad.setSkewComponent(0, 0.1);
        quad.setSkewComponent(1, -0.2);
        ElementBase& element = kind == 0   ? static_cast<ElementBase&>(sector)
                               : kind == 1 ? static_cast<ElementBase&>(rectangular)
                                           : static_cast<ElementBase&>(quad);
        const CoordinateSystemTrafo frame(
                Vector(1., -2., 0.3), Quaternion(std::cos(0.2), 0., std::sin(0.2), 0.));
        element.setCSTrafoGlobal2Local(frame);
        element.fixPosition();
        element.setAperture(ApertureType::ELLIPTICAL, {0.2, 0.1});
        lattice.visit(element, visitor, bunch);
        lattice.prepareSections();
        ASSERT_TRUE(device_external::Builder::supports(lattice));
        const auto device = device_external::Builder::build(lattice);
        Kokkos::View<Vector*> points("points", 7), fields("fields", 7);
        auto input           = Kokkos::create_mirror_view(points);
        const Vector local[] = {Vector(0.01, 0.02, -0.02), Vector(0.01, 0.02, 0.02),
                                Vector(-0.06, 0.01, 0.49), Vector(-0.23, 0.01, 0.97),
                                Vector(0, 0, 1.04),        Vector(0.5, 0, 0.5),
                                Vector(0, 0.2, 0.5)};
        for (unsigned i = 0; i < 7; ++i)
            input(i) = frame.transformFrom(local[i]);
        Kokkos::deep_copy(points, input);
        Kokkos::parallel_for("test spatial field", 7, MagneticFieldKernel{device, points, fields});
        const auto output = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), fields);
        for (unsigned i = 0; i < 7; ++i) {
            Vector electric(0), expected(0);
            lattice.getFieldAt(input(i), Vector(0, 0, 1), 0., electric, expected);
            for (unsigned d = 0; d < 3; ++d)
                EXPECT_NEAR(output(i)(d), expected(d), 2e-13);
        }
    }
}

TEST_F(DeviceExternalFieldTest, FrameTransformAndReverseDriftStayOnDevice) {
    device_external::Lattice lattice;
    const CoordinateSystemTrafo frame(
            Vector(1, 2, 3), Quaternion(std::cos(0.3), 0., std::sin(0.3), 0.));
    const device_external::Rigid deviceFrame{frame.getOrigin(), frame.getRotationMatrix()};
    Kokkos::View<Vector*> r("r", 1), p("p", 1), e("e", 1), b("b", 1);
    const Vector start(0.01, 0.02, 0.03), momentum(0.1, 0.2, 1.);
    Kokkos::deep_copy(r, start);
    Kokkos::deep_copy(p, momentum);
    ASSERT_EQ(lattice.transport(r, p, e, b, 1, deviceFrame, 1e-10, 1e9, -1.), 0);
    auto forward = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), r);
    for (unsigned d = 0; d < 3; ++d)
        EXPECT_NEAR(
                forward(0)(d),
                start(d)
                        + Physics::c * 1e-10 * momentum(d) / std::sqrt(1 + dot(momentum, momentum)),
                2e-15);
    ASSERT_EQ(lattice.transport(r, p, e, b, 1, deviceFrame, -1e-10, 1e9, -1.), 0);
    const auto backward = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), r);
    for (unsigned d = 0; d < 3; ++d)
        EXPECT_NEAR(backward(0)(d), start(d), 2e-15);
    EXPECT_EQ(lattice.transport(r, p, e, b, 1, deviceFrame, 1e-10, 0, 1), 1);
}

TEST_F(DeviceExternalFieldTest, UniformRFCavityGeometryUsesItsWidthAndRejectsMagneticOnlyOracle) {
    Beam beam;
    PartBunch_t bunch(
            {1.}, {1.}, {&beam}, {0}, 1., "LF2", opalx::spacecharge::CartesianDomainConfig3D{});
    FlaggedBeamline line;
    DefaultVisitor visitor(line, false, false);
    OpalBeamline lattice;
    VariableRFCavity cavity("buncher");
    cavity.setLength(0.1);
    cavity.setWidth(0.4);
    cavity.setHeight(0.2);
    // Generic aperture loss is distinct from RF WIDTH/HEIGHT field support.
    cavity.setAperture(ApertureType::RECTANGULAR, {1e-5, 1e-5});
    cavity.setCSTrafoGlobal2Local(CoordinateSystemTrafo());
    cavity.fixPosition();
    lattice.visit(cavity, visitor, bunch);
    lattice.prepareSections();
    ASSERT_TRUE(device_external::Builder::supports(lattice));
    const auto device = device_external::Builder::build(lattice);
    EXPECT_FALSE(device.magneticOnly);
    EXPECT_DOUBLE_EQ(device.maximumStep, 0.1 / (4 * Physics::c));
    Kokkos::View<Vector*> points("RF support points", 8);
    auto input = Kokkos::create_mirror_view(points);
    input(0)   = Vector(0, 0, 0);
    input(1)   = Vector(0, 0, 0.1);
    input(2)   = Vector(0.2, 0, 0.05);
    input(3)   = Vector(-0.2, 0.1, 0.05);
    input(4)   = Vector(0.201, 0, 0.05);
    input(5)   = Vector(0, -0.101, 0.05);
    input(6)   = Vector(0, 0, -1e-8);
    input(7)   = Vector(0.1, 0.01, 0.05);
    Kokkos::deep_copy(points, input);
    Kokkos::View<int*> flags("RF support flags", 8);
    Kokkos::parallel_for("RF support check", 8, SupportKernel{device.elements, points, flags});
    const auto actual = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), flags);
    for (int i = 0; i < 8; ++i)
        EXPECT_EQ(bool(actual(i)), cavity.isInside(input(i)));
    device_external::State state;
    state.momentum = Vector(0, 0, 1);
    EXPECT_FALSE(device.advance(state, 1e-12, 9.38e8, 1));
    auto stored = std::dynamic_pointer_cast<VariableRFCavity>(*lattice.getElements().begin());
    ASSERT_TRUE(stored);
    stored->setHeight(0);
    EXPECT_FALSE(device_external::Builder::supports(lattice));
    EXPECT_THROW(device_external::Builder::build(lattice), OpalException);
}

TEST_F(DeviceExternalFieldTest, RejectsUnimplementedDeviceFieldInsteadOfHostFallback) {
    Beam beam;
    PartBunch_t bunch(
            {1.}, {1.}, {&beam}, {0}, 1., "LF2", opalx::spacecharge::CartesianDomainConfig3D{});
    FlaggedBeamline line;
    DefaultVisitor visitor(line, false, false);
    OpalBeamline lattice;
    MultipoleRep sextupole("sextupole");
    sextupole.getGeometry().setElementLength(1);
    sextupole.setNormalComponent(2, 0.5);
    sextupole.setCSTrafoGlobal2Local(CoordinateSystemTrafo());
    sextupole.fixPosition();
    lattice.visit(sextupole, visitor, bunch);
    lattice.prepareSections();
    EXPECT_FALSE(device_external::Builder::supports(lattice));
    EXPECT_THROW(device_external::Builder::build(lattice), OpalException);
    auto stored = *lattice.getElements().begin();
    dynamic_cast<Multipole&>(*stored).setNormalComponent(2, 0.);
    EXPECT_TRUE(device_external::Builder::supports(lattice));
    dynamic_cast<Multipole&>(*stored).setSkewComponent(3, 0.25);
    EXPECT_FALSE(device_external::Builder::supports(lattice));
    EXPECT_THROW(device_external::Builder::build(lattice), OpalException);
    dynamic_cast<Multipole&>(*stored).setSkewComponent(3, 0.);
    EXPECT_TRUE(device_external::Builder::supports(lattice));
    stored->setMisalignment(CoordinateSystemTrafo(Vector(1e-4, 0, 0), Quaternion()));
    EXPECT_FALSE(device_external::Builder::supports(lattice));
    EXPECT_THROW(device_external::Builder::build(lattice), OpalException);
    stored->setMisalignment(
            CoordinateSystemTrafo(Vector(0), Quaternion(std::cos(0.1), 0., std::sin(0.1), 0.)));
    EXPECT_FALSE(device_external::Builder::supports(lattice));
    EXPECT_THROW(device_external::Builder::build(lattice), OpalException);
    stored->setMisalignment(CoordinateSystemTrafo());
    EXPECT_TRUE(device_external::Builder::supports(lattice));
}

TEST_F(DeviceExternalFieldTest,
       EligibilityAcceptsPassivesAndRejectsOtherFieldModelsWithoutMutation) {
    Beam beam;
    PartBunch_t bunch(
            {1.}, {1.}, {&beam}, {0}, 1., "LF2", opalx::spacecharge::CartesianDomainConfig3D{});
    FlaggedBeamline line;
    DefaultVisitor visitor(line, false, false);
    OpalBeamline lattice;
    DriftRep drift("drift");
    drift.getGeometry().setElementLength(1.);
    drift.setCSTrafoGlobal2Local(CoordinateSystemTrafo());
    drift.fixPosition();
    MarkerRep marker("marker");
    marker.setCSTrafoGlobal2Local(CoordinateSystemTrafo());
    marker.fixPosition();
    lattice.visit(drift, visitor, bunch);
    lattice.visit(marker, visitor, bunch);
    lattice.prepareSections();
    const auto before = lattice.getElements();
    ASSERT_EQ(before.size(), 2u);
    EXPECT_TRUE(device_external::Builder::supports(lattice));
    EXPECT_TRUE(device_external::Builder::supports(lattice));
    EXPECT_EQ(lattice.getElements(), before);
    EXPECT_NO_THROW(device_external::Builder::build(lattice));

    ConstantFocusingRep other("other-field-model");
    other.getGeometry().setElementLength(1.);
    other.setCSTrafoGlobal2Local(CoordinateSystemTrafo());
    other.fixPosition();
    lattice.visit(other, visitor, bunch);
    lattice.prepareSections();
    const auto unsupported = lattice.getElements();
    EXPECT_FALSE(device_external::Builder::supports(lattice));
    EXPECT_EQ(lattice.getElements(), unsupported);
    EXPECT_THROW(device_external::Builder::build(lattice), OpalException);
}

TEST_F(DeviceExternalFieldTest, CyclotronSectorRetainsItsMappedTrackingPath) {
    // A real mapped sector must remain on its established device field path.
    // Use the same synthetic PSI-map layout as TestCyclotronSector.
    const auto filename = std::filesystem::temp_directory_path()
                          / ("opalx-ring-eligibility-" + std::to_string(getpid()) + ".map");
    {
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
                    out << (channel == 0 ? 10. : 0.) << ' ';
            out << '\n';
        }
    }
    const auto map = CyclotronSectorFieldMap::read(filename.string());
    std::filesystem::remove(filename);
    Beam beam;
    PartBunch_t bunch(
            {1.}, {1.}, {&beam}, {0}, 1., "LF2", opalx::spacecharge::CartesianDomainConfig3D{});
    FlaggedBeamline line;
    DefaultVisitor visitor(line, false, false);
    OpalBeamline lattice;
    CyclotronSector sector("mapped-cyclotron-sector");
    sector.configure(map, 8, -0.05, 0.05, 1., {});
    sector.setCSTrafoGlobal2Local(CoordinateSystemTrafo());
    sector.fixPosition();
    lattice.visit(sector, visitor, bunch);
    lattice.prepareSections();
    const auto before = lattice.getElements();
    ASSERT_EQ(before.size(), 1u);
    ASSERT_EQ((*before.begin())->getType(), ElementType::CYCLOTRONSECTOR);
    EXPECT_FALSE(device_external::Builder::supports(lattice));
    EXPECT_THROW(device_external::Builder::build(lattice), OpalException);
    EXPECT_EQ(lattice.getElements(), before);
    Vector electric(0.), magnetic(0.);
    (*before.begin())->apply(Vector(0., 0., 0.1), Vector(0., 0., 0.4), 0., electric, magnetic);
    EXPECT_NEAR(magnetic(1), 1., 1e-14);
}

TEST_F(DeviceExternalFieldTest, EligibilityRejectsHigherBendComponentsButAllowsZeroPadding) {
    Beam beam;
    PartBunch_t bunch(
            {1.}, {1.}, {&beam}, {0}, 1., "LF2", opalx::spacecharge::CartesianDomainConfig3D{});
    for (unsigned kind = 0; kind < 2; ++kind) {
        FlaggedBeamline line;
        DefaultVisitor visitor(line, false, false);
        OpalBeamline lattice;
        SBendRep sector("sector");
        RBendRep rectangular("rectangular");
        sector.getGeometry() = Geometry::makeSBend(1., 0.2);
        rectangular.getGeometry().setElementLength(1.);
        sector.setFieldComponents({1., 0.2, 0.}, {0.3, 0.4, 0.});
        rectangular.setFieldComponents({1., 0.2, 0.}, {0.3, 0.4, 0.});
        ElementBase& element = kind == 0 ? static_cast<ElementBase&>(sector)
                                         : static_cast<ElementBase&>(rectangular);
        element.setCSTrafoGlobal2Local(CoordinateSystemTrafo());
        element.fixPosition();
        lattice.visit(element, visitor, bunch);
        lattice.prepareSections();
        ASSERT_TRUE(device_external::Builder::supports(lattice));
        EXPECT_NO_THROW(device_external::Builder::build(lattice));
        const auto stored = *lattice.getElements().begin();
        for (bool normal : {true, false}) {
            const std::vector<double> normalCoefficients = {1., 0.2, normal ? 0.1 : 0.};
            const std::vector<double> skewCoefficients   = {0.3, 0.4, normal ? 0. : 0.1};
            if (kind == 0)
                dynamic_cast<SBend&>(*stored).setFieldComponents(
                        normalCoefficients, skewCoefficients);
            else
                dynamic_cast<RBend&>(*stored).setFieldComponents(
                        normalCoefficients, skewCoefficients);
            EXPECT_FALSE(device_external::Builder::supports(lattice));
            EXPECT_THROW(device_external::Builder::build(lattice), OpalException);
        }
    }
}
