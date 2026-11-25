//
// Tests of EmittingSamplers: the containers emission keeps active, the emission time step and
// progress the tracker reads, and emit() with samplers that append particles.
//
#include <gtest/gtest.h>

#include "Ippl.h"

#include "../BunchBeam.h"
#include "../KeywordObjects.h"
#include "Distribution/EmittingSamplers.h"
#include "Distribution/SamplingBase.hpp"
#include "Elements/OpalBeamline.h"
#include "PartBunch/PartBunch.h"
#include "Structure/Beam.h"

#include <memory>
#include <utility>
#include <vector>

extern Inform* gmsg;

namespace {
    using Vector = Vector_t<double, 3>;

    /// Emits until @p endTime with a fixed emitted fraction; each call appends one particle at
    /// `position` in the beamline frame.
    class StubSampler : public SamplingBase {
    public:
        StubSampler(
                std::shared_ptr<ParticleContainer_t> pc, double fraction, double endTime,
                double timeStep = 0.0)
            : SamplingBase(std::move(pc)),
              fraction_m(fraction),
              endTime_m(endTime),
              timeStep_m(timeStep) {}

        void emitParticles(double, double) override {
            ++calls;
            const size_t index = pc_m->getLocalNum();
            pc_m->createParticles(1);
            auto R   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc_m->R.getView());
            auto P   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc_m->P.getView());
            R(index) = position;
            P(index) = Vector({0.0, 0.0, 1.0});
            Kokkos::deep_copy(pc_m->R.getView(), R);
            Kokkos::deep_copy(pc_m->P.getView(), P);
        }
        bool isEmissionDone(double t) const override { return t >= endTime_m; }
        double getEmittedFraction() const override { return fraction_m; }
        double getEmissionTimeStep() const override { return timeStep_m; }

        const Vector position = Vector({0.01, -0.02, 0.003});
        int calls             = 0;

    private:
        double fraction_m, endTime_m, timeStep_m;
    };

    class EmittingSamplersTest : public ::testing::Test {
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
                    std::vector(2, opalx::test::bunchBeam(*beam, 8)), 1.0,
                    opalx::spacecharge::CartesianDomainConfig3D{});
        }

        std::shared_ptr<StubSampler> stub(
                size_t container, double fraction, double endTime, double timeStep = 0.0) {
            return std::make_shared<StubSampler>(
                    bunch_m->getParticleContainer(container), fraction, endTime, timeStep);
        }

        std::unique_ptr<PartBunch_t> bunch_m;
    };

    TEST_F(EmittingSamplersTest, ContainerIsActiveWhileItsSamplerEmits) {
        EmittingSamplers samplers({{}, {stub(1, 0.5, 1e-9)}});
        bunch_m->resetPcActive();
        samplers.activateEmittingContainers(*bunch_m, 0.0);
        EXPECT_FALSE(bunch_m->isPcActive(0));
        EXPECT_TRUE(bunch_m->isPcActive(1));

        bunch_m->resetPcActive();
        samplers.activateEmittingContainers(*bunch_m, 1e-9);
        EXPECT_FALSE(bunch_m->isPcActive(1));
    }

    TEST_F(EmittingSamplersTest, EmissionTimeStepIsTheSmallestOfTheEmittingSamplers) {
        EXPECT_FALSE(EmittingSamplers().emissionTimeStep(0.0));
        // A sampler without a time step of its own leaves the schedule's.
        EXPECT_FALSE(EmittingSamplers({{stub(0, 0.5, 2e-9)}}).emissionTimeStep(0.0));

        EmittingSamplers samplers(
                {{stub(0, 0.5, 2e-9, 3e-12), stub(0, 0.5, 2e-9)}, {stub(1, 0.5, 1e-9, 2e-12)}});
        EXPECT_EQ(samplers.emissionTimeStep(0.0), 2e-12);
        EXPECT_EQ(samplers.emissionTimeStep(1.5e-9), 3e-12);
        EXPECT_FALSE(samplers.emissionTimeStep(2e-9));
    }

    TEST_F(EmittingSamplersTest, SpaceChargeProgressIsTheLeastEmittedFraction) {
        EmittingSamplers samplers({{stub(0, 0.3, 1e-9)}, {stub(1, 0.6, 2e-9), stub(1, 1.0, 3e-9)}});
        for (const bool sync : {false, true}) {
            auto progress = samplers.spaceChargeEmissionProgress(0.0, sync);
            EXPECT_TRUE(progress.active);
            EXPECT_EQ(progress.fraction, 0.3);
            progress = samplers.spaceChargeEmissionProgress(1.5e-9, sync);
            EXPECT_TRUE(progress.active);
            EXPECT_EQ(progress.fraction, 0.6);
            // Only a sampler that emitted everything is left.
            progress = samplers.spaceChargeEmissionProgress(2.5e-9, sync);
            EXPECT_FALSE(progress.active);
            EXPECT_EQ(progress.fraction, 1.0);
        }
    }

    TEST_F(EmittingSamplersTest, EmitAsksEverySamplerOnceAndLeavesTrackedParticles) {
        auto pc = bunch_m->getParticleContainer(0);
        // The reference frame is away from the beamline frame the samplers write in.
        const CoordinateSystemTrafo toLab(
                Vector({0.1, -0.2, 0.5}),
                getQuaternion(Vector({0.3, 0.0, 1.0}), Vector({0, 0, 1})));
        pc->setToLabTrafo(toLab);
        pc->createParticles(2);
        auto R = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
        auto P = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
        const std::vector<Vector> trackedR{Vector({1e-3, 2e-3, -3e-3}), Vector({-4e-3, 0.0, 5e-3})};
        const std::vector<Vector> trackedP{Vector({0.01, 0.0, 2.0}), Vector({0.0, -0.02, 2.1})};
        for (size_t i = 0; i < 2; ++i) {
            R(i) = trackedR[i];
            P(i) = trackedP[i];
        }
        Kokkos::deep_copy(pc->R.getView(), R);
        Kokkos::deep_copy(pc->P.getView(), P);

        const auto first = stub(0, 0.5, 1e-9), second = stub(0, 0.5, 1e-9);
        const auto stopped = stub(1, 0.5, 1e-9);
        EmittingSamplers samplers({{first, second}, {stopped}});
        bunch_m->resetPcActive();
        bunch_m->setPcAtSStop(1);
        const OpalBeamline beamline;
        for (int step = 1; step <= 2; ++step) {
            samplers.emit(*bunch_m, beamline, step * 1e-12, 1e-12);
            EXPECT_EQ(first->calls, step);
            EXPECT_EQ(second->calls, step);
        }
        EXPECT_EQ(stopped->calls, 0);
        EXPECT_EQ(bunch_m->getParticleContainer(1)->getTotalNum(), 0u);
        ASSERT_EQ(pc->getTotalNum(), 6u);
        EXPECT_TRUE(bunch_m->isPcActive(0));

        R = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
        P = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
        for (size_t i = 0; i < 2; ++i) {
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_EQ(R(i)(d), trackedR[i](d));
                EXPECT_EQ(P(i)(d), trackedP[i](d));
            }
        }
        // The new particles are where the samplers put them, seen from the reference frame.
        const Vector expected = toLab.transformFrom(first->position);
        for (size_t i = 2; i < 6; ++i) {
            for (unsigned d = 0; d < 3; ++d)
                EXPECT_NEAR(R(i)(d), expected(d), 1e-14);
        }
    }

    // With per-particle charge and mass (QM_MODE = ATTRIBUTES) the particles emitted into an
    // empty container get the container's charge and mass.
    TEST_F(EmittingSamplersTest, Attributes_EmittedParticlesKeepCharge) {
        const auto beam = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        PartBunch_t bunch(
                std::vector{opalx::test::bunchBeam(*beam, 8, 2.0, 3.0)}, 1.0,
                opalx::spacecharge::CartesianDomainConfig3D{}, /*qmAttributes=*/true);
        auto pc = bunch.getParticleContainer(0);
        ASSERT_EQ(pc->getQMStorageMode(), ParticleContainer_t::QMStorageMode::Attributes);
        ASSERT_EQ(pc->getTotalNum(), 0u);

        // Each stub appends one particle per emit().
        EmittingSamplers samplers(
                {{std::make_shared<StubSampler>(pc, 0.5, 1e-9),
                  std::make_shared<StubSampler>(pc, 0.5, 1e-9)}});
        bunch.resetPcActive();
        samplers.emit(bunch, OpalBeamline(), 1e-12, 1e-12);
        ASSERT_EQ(pc->getTotalNum(), 2u);

        EXPECT_DOUBLE_EQ(pc->getTotalCharge(), 2.0 * 2.0);
        auto q = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->getQView());
        auto m = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->getMView());
        for (size_t i = 0; i < 2; ++i) {
            EXPECT_DOUBLE_EQ(q(i), 2.0) << "particle " << i;
            EXPECT_DOUBLE_EQ(m(i), 3.0) << "particle " << i;
        }
    }
}  // namespace
