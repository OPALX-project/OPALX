//
// Tests of the OpalBeamline dispatch to the elements near a container: the external fields in
// each element's frame and the aperture loss marks.
//
#include <gtest/gtest.h>

#include "Ippl.h"

#include "../KeywordObjects.h"
#include "BeamlineCore/MultipoleRep.h"
#include "Elements/OpalBeamline.h"
#include "PartBunch/PartBunch.h"
#include "Structure/Beam.h"

#include <cmath>
#include <memory>
#include <vector>

extern Inform* gmsg;

namespace {
    using Vector = Vector_t<double, 3>;

    class OpalBeamlineDispatchTest : public ::testing::Test {
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

        void SetUp() override {
            const auto beam = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
            bunch_m         = std::make_unique<PartBunch_t>(
                    std::vector<double>{1.0}, std::vector<double>{1.0},
                    std::vector{opalx::run::describeBeamForBunch(*beam)}, std::vector<size_t>{8},
                    1.0, "LF2", opalx::spacecharge::CartesianDomainConfig3D{});
            pc_m = bunch_m->getParticleContainer(0);
            // The reference frame sits 0.3 m downstream of the lab origin.
            pc_m->setToLabTrafo(toLab_m);
            pc_m->createParticles(positions_m.size());
            auto R = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc_m->R.getView());
            auto P = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc_m->P.getView());
            for (size_t i = 0; i < positions_m.size(); ++i) {
                R(i) = positions_m[i];
                P(i) = Vector({0.0, 0.0, 1.0});
            }
            Kokkos::deep_copy(pc_m->R.getView(), R);
            Kokkos::deep_copy(pc_m->P.getView(), P);
            bunch_m->resetPcActive();
            bunch_m->setTimeStep(1e-12);

            // A quadrupole over lab z in [0.2, 0.5) m with a 5 mm bore.
            quadrupole_m->getGeometry().setElementLength(0.3);
            quadrupole_m->setNormalComponent(1, 2.0);
            quadrupole_m->setAperture(ApertureType::RECTANGULAR, {5e-3, 5e-3});
            quadrupole_m->setCSTrafoGlobal2Local(
                    CoordinateSystemTrafo(Vector({0.0, 0.0, 0.2}), Quaternion()));
            quadrupole_m->fixPosition();
            beamline_m.addElement(quadrupole_m, *bunch_m);
            beamline_m.prepareSections();
            beamline_m.buildElementMapping(defaultElementMappingGrid, defaultElementReach);
        }

        template <class Attribute>
        static auto host(const Attribute& attribute) {
            return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), attribute.getView());
        }

        const CoordinateSystemTrafo toLab_m{Vector({0.0, 0.0, -0.3}), Quaternion()};
        // Reference frame; local z of the quadrupole is z + 0.1. The first three are in its
        // body, the third outside its bore; the last two are in the drift behind it.
        const std::vector<Vector> positions_m{
                Vector({1e-3, 2e-3, 0.0}), Vector({-2e-3, 1e-3, 0.1}), Vector({8e-3, 0.0, 0.05}),
                Vector({1e-3, -1e-3, 0.4}), Vector({9e-3, 0.0, 0.45})};
        std::unique_ptr<PartBunch_t> bunch_m;
        std::shared_ptr<ParticleContainer_t> pc_m;
        std::shared_ptr<MultipoleRep> quadrupole_m = std::make_shared<MultipoleRep>("Q");
        OpalBeamline beamline_m;
    };

    TEST_F(OpalBeamlineDispatchTest, ExternalFieldsAreTheElementFieldInItsFrame) {
        beamline_m.applyExternalFields(*bunch_m);
        EXPECT_TRUE(pc_m->frame().isReference());

        auto R = host(pc_m->R), E = host(pc_m->E), B = host(pc_m->B);
        for (size_t i = 0; i < positions_m.size(); ++i) {
            const Vector local = quadrupole_m->getCSTrafoGlobal2Local().transformTo(
                    toLab_m.transformTo(positions_m[i]));
            Vector expectedE(0.0), expectedB(0.0);
            quadrupole_m->apply(local, Vector({0.0, 0.0, 1.0}), 0.0, expectedE, expectedB);
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_NEAR(R(i)(d), positions_m[i](d), 1e-15) << i;
                EXPECT_EQ(E(i)(d), 0.0) << i;
                EXPECT_NEAR(B(i)(d), expectedB(d), 1e-15) << i;
                if (i >= 2) EXPECT_EQ(B(i)(d), 0.0) << i;
            }
        }
        // The particles in the bore feel the gradient.
        EXPECT_GT(std::abs(B(0)(0)), 1e-3);
        EXPECT_GT(std::abs(B(1)(1)), 1e-3);
    }

    TEST_F(OpalBeamlineDispatchTest, AperturesMarkOnlyParticlesOutsideTheBoreOfTheirElement) {
        const auto before = host(pc_m->R);
        EXPECT_EQ(beamline_m.applyElementApertures(*bunch_m), 1u);
        // Marks are counted once.
        EXPECT_EQ(beamline_m.applyElementApertures(*bunch_m), 0u);
        EXPECT_TRUE(pc_m->frame().isReference());

        const auto after   = host(pc_m->R);
        const auto invalid = host(pc_m->InvalidMask);
        for (size_t i = 0; i < positions_m.size(); ++i) {
            EXPECT_EQ(invalid(i), i == 2) << i;
            for (unsigned d = 0; d < 3; ++d)
                EXPECT_EQ(after(i)(d), before(i)(d)) << i;
        }
    }
}  // namespace
