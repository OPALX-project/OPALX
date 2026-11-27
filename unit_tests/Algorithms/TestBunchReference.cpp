//
// Tests of BunchReference: where each beam's reference particle starts and how it advances in a
// field-free line.
//
#include <gtest/gtest.h>

#include "Ippl.h"

#include "../BunchBeam.h"
#include "../KeywordObjects.h"
#include "AbsBeamline/Drift.h"
#include "Algorithms/BunchReference.h"
#include "Distribution/EmittingSamplers.h"
#include "Distribution/SamplingBase.hpp"
#include "Elements/OpalBeamline.h"
#include "Moments/BeamMoments.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"
#include "SpaceCharge/SpaceChargeAlgorithm.h"
#include "SpaceCharge/SpaceChargeSolver.h"
#include "Steppers/NumericalIntegrator.h"
#include "Structure/Beam.h"

#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <vector>

extern Inform* gmsg;

namespace {
    using Vector = Vector_t<double, 3>;

    /// Reports the emission offsets R0 and P0 of a source.
    class OffsetSampler : public SamplingBase {
    public:
        using SamplingBase::SamplingBase;
        bool hasInitialReferencePosition() const override { return true; }
        Vector getInitialReferencePosition() const override { return Vector({1e-3, 0.0, 0.02}); }
        bool hasInitialReferenceMomentum() const override { return true; }
        Vector getInitialReferenceMomentum() const override { return Vector({0.0, 0.01, 0.2}); }
    };

    /// A space-charge algorithm without self field; the references see none anyway.
    class NoSelfField final : public opalx::spacecharge::SpaceChargeAlgorithm {
    public:
        opalx::spacecharge::SpaceChargeSolveResult solve(
                const opalx::spacecharge::SpaceChargeSolveContext&) override {
            return {};
        }
        opalx::spacecharge::SolveFrame solveFrame() const override {
            return opalx::spacecharge::SolveFrame::Reference;
        }
    };

    class BunchReferenceTest : public ::testing::Test {
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

        // Three beams of protons at PC = 0.1 GeV; only the first holds particles.
        void SetUp() override {
            const auto beam = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
            Attributes::setPredefinedString(*beam->findAttribute("PARTICLE"), "PROTON");
            Attributes::setReal(*beam->findAttribute("PC"), 0.1);
            bunch_m = std::make_unique<PartBunch_t>(
                    std::vector(3, opalx::test::bunchBeam(*beam, 8)), 1.0,
                    opalx::spacecharge::CartesianDomainConfig3D{});
            auto pc = bunch_m->getParticleContainer(0);
            pc->createParticles(3);
            auto R = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
            auto P = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
            for (size_t i = 0; i < 3; ++i) {
                R(i) = Vector({1e-3 * i, -2e-3, 0.01 + 1e-3 * i});
                P(i) = Vector({0.001 * i, 0.0, 0.1 + 0.01 * i});
            }
            Kokkos::deep_copy(pc->R.getView(), R);
            Kokkos::deep_copy(pc->P.getView(), P);
            bunch_m->resetPcActive();
        }

        /// The integration scheme @p name, whose reference step BunchReference uses.
        std::unique_ptr<NumericalIntegrator> scheme(const std::string& name) {
            return makeNumericalIntegrator(
                    name, {*bunch_m, beamline_m, solver_m, samplers_m, {}, 0, false});
        }

        static void expectNear(const Vector& actual, const Vector& expected, double tolerance) {
            for (unsigned d = 0; d < 3; ++d)
                EXPECT_NEAR(actual(d), expected(d), tolerance) << d;
        }

        std::unique_ptr<PartBunch_t> bunch_m;
        OpalBeamline beamline_m;
        opalx::spacecharge::SpaceChargeSolver solver_m{std::make_unique<NoSelfField>(), 3};
        EmittingSamplers samplers_m;
    };

    TEST_F(BunchReferenceTest, InitialPoseFromMeanThenSamplerThenDesign) {
        const auto mean = opalx::moments::compute(
                *bunch_m->getParticleContainer(0),
                opalx::moments::Moment::MeanPosition | opalx::moments::Moment::MeanMomentum);
        const auto sampler = std::make_shared<OffsetSampler>(bunch_m->getParticleContainer(1));
        const auto boris   = scheme("BORIS");
        BunchReference references(*bunch_m, beamline_m, *boris, {}, 0, {});
        references.setInitialPoses(EmittingSamplers({{}, {sampler}, {}}));

        const auto& withParticles = *bunch_m->getParticleContainer(0);
        expectNear(withParticles.getRefPartR(), mean.meanR(), 1e-18);
        expectNear(withParticles.getRefPartP(), mean.meanP(), 1e-16);

        const auto& emitting = *bunch_m->getParticleContainer(1);
        expectNear(emitting.getRefPartR(), sampler->getInitialReferencePosition(), 0.0);
        expectNear(emitting.getRefPartP(), sampler->getInitialReferenceMomentum(), 0.0);

        const auto& design        = *bunch_m->getParticleContainer(2);
        const PartData& reference = *design.getReference();
        expectNear(design.getRefPartR(), Vector(0.0), 0.0);
        expectNear(
                design.getRefPartP(), Vector({0.0, 0.0, reference.getP() / reference.getM()}), 0.0);
        EXPECT_NEAR(design.getRefPartP()(2), 0.1 / 0.93827208816, 1e-6);
    }

    TEST_F(BunchReferenceTest, BeamReferenceMomentumOverridesTheMean) {
        const auto mean = opalx::moments::compute(
                *bunch_m->getParticleContainer(0), opalx::moments::Moment::MeanPosition);
        const auto boris = scheme("BORIS");
        BunchReference references(*bunch_m, beamline_m, *boris, {Vector({0.0, 0.0, 0.5})}, 0, {});
        references.setInitialPoses(EmittingSamplers());

        const auto& pc = *bunch_m->getParticleContainer(0);
        expectNear(pc.getRefPartR(), mean.meanR(), 1e-18);
        expectNear(pc.getRefPartP(), Vector({0.0, 0.0, 0.5}), 0.0);
    }

    TEST_F(BunchReferenceTest, AdvanceInADriftMovesByBetaCDt) {
        const auto drift = std::make_shared<Drift>("D");
        drift->getGeometry().setElementLength(2.0);
        drift->setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector(0.0), Quaternion()));
        drift->fixPosition();
        beamline_m.addElement(drift, *bunch_m);
        beamline_m.prepareSections();
        beamline_m.buildElementMapping(defaultElementMappingGrid, defaultElementReach);

        // The reference sits at the origin of its frame and moves obliquely, with the reference
        // step of each scheme.
        auto pc         = bunch_m->getParticleContainer(0);
        const Vector p0 = Vector({0.03, 0.0, 0.12});
        const double dt = 1e-11;
        bunch_m->setTimeStep(dt);
        for (const std::string name : {"BORIS", "RK4"}) {
            SCOPED_TRACE(name);
            pc->setToLabTrafo(CoordinateSystemTrafo());
            pc->setRefPartR(Vector(0.0));
            pc->setRefPartP(p0);
            pc->set_sPos(0.3);
            const auto integrator = scheme(name);
            BunchReference references(*bunch_m, beamline_m, *integrator, {}, 0, {});

            // Nothing moves before the bunch time is positive.
            bunch_m->setT(0.0);
            references.advance();
            expectNear(pc->getRefPartR(), Vector(0.0), 0.0);
            EXPECT_EQ(pc->get_sPos(), 0.3);

            bunch_m->setT(dt);
            references.advance();
            const double norm = std::sqrt(dot(p0, p0));
            const double ds   = norm / std::sqrt(1.0 + norm * norm) * Physics::c * dt;
            expectNear(pc->getRefPartR(), ds / norm * p0, 1e-16);
            expectNear(pc->getRefPartP(), p0, 1e-16);
            EXPECT_NEAR(pc->get_sPos(), 0.3 + ds, 1e-15);
            EXPECT_FALSE(references.hitElement());
        }
    }
}  // namespace
