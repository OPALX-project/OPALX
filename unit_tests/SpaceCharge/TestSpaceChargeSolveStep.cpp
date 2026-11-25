//
// Tests of SpaceChargeSolver::solveStep(): the MINBINEMITTED threshold, the step state and solve
// frame it hands to the algorithm, and the reference frame of the containers on return.
//
#include <gtest/gtest.h>

#include "Ippl.h"

#include "../BunchBeam.h"
#include "../KeywordObjects.h"
#include "../TemporaryDirectory.h"
#include "Moments/BeamMoments.h"
#include "PartBunch/PartBunch.h"
#include "SpaceCharge/SpaceChargeFactory.h"
#include "SpaceCharge/SpaceChargeSolver.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Utilities/OpalException.h"

#include <memory>
#include <optional>
#include <vector>

extern Inform* gmsg;

namespace opalx::spacecharge {
    namespace {
        using Vector = Vector_t<double, 3>;

        /// Records the step it is asked to solve; enters the solve frame and, unless told to
        /// stay there, returns to the reference frame as an algorithm must.
        class FrameAlgorithm final : public SpaceChargeAlgorithm {
        public:
            FrameAlgorithm(ParticleContainer_t& primary, SolveFrame frame, bool stay = false)
                : primary_m(primary), frame_m(frame), stay_m(stay) {}

            SpaceChargeSolveResult solve(const SpaceChargeSolveContext& context) override {
                ++calls;
                state = context.stepState();
                activity.assign(context.trackingActive().begin(), context.trackingActive().end());
                if (stay_m) {
                    primary_m.enterFrame(context.stepState().solveFrame);
                } else {
                    const ParticleContainer_t::FrameScope scope(
                            primary_m, context.stepState().solveFrame);
                }
                return {.backendSolves = 1};
            }
            SolveFrame solveFrame() const override { return frame_m; }

            int calls = 0;
            std::optional<SpaceChargeStepState> state;
            std::vector<std::uint8_t> activity;

        private:
            ParticleContainer_t& primary_m;
            SolveFrame frame_m;
            bool stay_m;
        };

        class SpaceChargeSolveStepTest : public ::testing::Test {
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
                directory_m.create("opalx-solve-step");
                const auto beam = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
                bunch_m         = std::make_unique<PartBunch_t>(
                        std::vector(2, opalx::test::bunchBeam(*beam, 8, 1e-15, 0.938)), 1.0,
                        domain_m);
                const std::vector<Vector> positions{
                        Vector({1e-3, 0.0, -1e-3}), Vector({-1e-3, 5e-4, 0.0}),
                        Vector({0.0, -1e-3, 1e-3}), Vector({5e-4, 1e-3, 2e-3})};
                for (const auto& pc : bunch_m->getParticleContainers()) {
                    pc->createParticles(positions.size());
                    auto R = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), pc->R.getView());
                    auto P = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), pc->P.getView());
                    for (size_t i = 0; i < positions.size(); ++i) {
                        R(i) = positions[i];
                        P(i) = Vector({1e-3 * static_cast<double>(i), 0.0, 0.5});
                    }
                    Kokkos::deep_copy(pc->R.getView(), R);
                    Kokkos::deep_copy(pc->P.getView(), P);
                    pc->setRefPartP(Vector({0.0, 0.0, 0.5}));
                }
                bunch_m->resetPcActive();
                bunch_m->setTimeStep(1e-12);
                bunch_m->setT(5e-12);
                bunch_m->setGlobalTrackStep(5);
                sink_m = std::make_unique<DataSink>(
                        DataSinkOptions{}, std::vector<H5PartWrapper*>{}, 2,
                        (directory_m.path() / "sink").string(), "test.in");
            }

            /// The configured CartesianPIC3D solver with Poisson backend @p backend.
            std::unique_ptr<SpaceChargeSolver> makeSolver(PoissonSolverType backend) {
                CartesianPIC3DConfig config;
                config.backend                         = backend;
                config.grid.meshSize                   = domain_m.meshSize;
                config.grid.decomposition              = domain_m.decomposition;
                config.grid.boundingBoxIncreasePercent = domain_m.boundingBoxIncreasePercent;
                return makeSpaceChargeSolver(config, *bunch_m, sink_m.get());
            }

            void expectReferenceFrames() const {
                for (const auto& pc : bunch_m->getParticleContainers())
                    EXPECT_TRUE(pc->frame().isReference()) << pc->frame().name();
            }

            CartesianDomainConfig3D domain_m;
            opalx::test::TemporaryDirectory directory_m;
            std::unique_ptr<PartBunch_t> bunch_m;
            std::unique_ptr<DataSink> sink_m;
            int progressCalls_m                                = 0;
            const std::function<EmissionProgress()> progress_m = [this] {
                ++progressCalls_m;
                return EmissionProgress{.active = true, .fraction = 0.25};
            };
        };

        TEST_F(SpaceChargeSolveStepTest, SolvesOnlyAboveMinBinEmitted) {
            auto algorithm = std::make_unique<FrameAlgorithm>(
                    *bunch_m->getParticleContainer(0), SolveFrame::Reference);
            FrameAlgorithm& observer = *algorithm;
            SpaceChargeSolver solver(std::move(algorithm), 2);

            // Eight particles in all: MINBINEMITTED = 8 still skips the step.
            solver.solveStep(*bunch_m, 8, progress_m);
            EXPECT_EQ(observer.calls, 0);
            EXPECT_EQ(progressCalls_m, 0);
            EXPECT_EQ(solver.backendSolveCount(), 0u);

            bunch_m->setPcAtSStop(1);
            solver.solveStep(*bunch_m, 7, progress_m);
            EXPECT_EQ(observer.calls, 1);
            EXPECT_EQ(progressCalls_m, 1);
            EXPECT_EQ(solver.backendSolveCount(), 1u);
            EXPECT_EQ(observer.activity, (std::vector<std::uint8_t>{1, 0}));
            ASSERT_TRUE(observer.state);
            EXPECT_EQ(observer.state->step, 5u);
            EXPECT_EQ(observer.state->time, 5e-12);
            EXPECT_EQ(observer.state->timeStep, 1e-12);
            EXPECT_TRUE(observer.state->emissionActive);
            EXPECT_EQ(observer.state->emittedFraction, 0.25);
            EXPECT_TRUE(observer.state->solveFrame.isReference());
        }

        TEST_F(SpaceChargeSolveStepTest, CentroidFrameComesFromThePrimaryContainer) {
            auto& primary  = *bunch_m->getParticleContainer(0);
            auto algorithm = std::make_unique<FrameAlgorithm>(primary, SolveFrame::Centroid);
            FrameAlgorithm& observer = *algorithm;
            SpaceChargeSolver solver(std::move(algorithm), 2);
            const auto moments = moments::compute(primary, moments::Moment::MeanPosition);

            solver.solveStep(*bunch_m, 0, progress_m);
            ASSERT_TRUE(observer.state);
            const Vector centroid =
                    observer.state->solveFrame.fromReference().transformTo(moments.meanR());
            for (unsigned d = 0; d < 3; ++d)
                EXPECT_NEAR(centroid(d), 0.0, 1e-18);
            expectReferenceFrames();
        }

        TEST_F(SpaceChargeSolveStepTest, AlgorithmLeavingTheSolveFrameIsAnError) {
            auto& primary = *bunch_m->getParticleContainer(0);
            SpaceChargeSolver solver(
                    std::make_unique<FrameAlgorithm>(primary, SolveFrame::Centroid, true), 2);
            EXPECT_THROW(solver.solveStep(*bunch_m, 0, progress_m), OpalException);
            primary.enterReferenceFrame();
        }

        TEST_F(SpaceChargeSolveStepTest, TypeNoneLeavesParticlesAndFieldsAlone) {
            auto solver = makeSolver(PoissonSolverType::None);
            std::vector<decltype(Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), bunch_m->getParticleContainer(0)->R.getView()))>
                    before;
            for (const auto& pc : bunch_m->getParticleContainers())
                before.push_back(
                        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView()));

            solver->solveStep(*bunch_m, 0, progress_m);
            expectReferenceFrames();
            for (size_t c = 0; c < bunch_m->getNumParticleContainers(); ++c) {
                const auto& pc = bunch_m->getParticleContainer(c);
                auto R = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
                auto E = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->E.getView());
                auto B = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->B.getView());
                for (size_t i = 0; i < pc->getLocalNum(); ++i) {
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_EQ(R(i)(d), before[c](i)(d));
                        EXPECT_EQ(E(i)(d), 0.0);
                        EXPECT_EQ(B(i)(d), 0.0);
                    }
                }
            }
        }

        TEST_F(SpaceChargeSolveStepTest, OpenSolveReturnsEveryContainerToItsReferenceFrame) {
            auto solver = makeSolver(PoissonSolverType::Open);
            ASSERT_EQ(solver->solveFrame(), SolveFrame::Centroid);
            solver->solveStep(*bunch_m, 0, progress_m);
            EXPECT_GE(solver->backendSolveCount(), 1u);
            expectReferenceFrames();
        }
    }  // namespace
}  // namespace opalx::spacecharge
