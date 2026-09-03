#include <gtest/gtest.h>

#include "Algorithms/Quaternion.hpp"
#include "PartBunch/BunchStateHandler.h"
#include "SpaceCharge/ParticleFrameGuard.h"
#include "SpaceCharge/SelfFieldRequestPolicy.h"
#include "SpaceCharge/SelfFieldSystem.h"
#include "Utilities/OpalException.h"

#include <array>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace opalx::spacecharge {
    namespace {

        class RecordingAlgorithm final : public SelfFieldAlgorithm {
        public:
            explicit RecordingAlgorithm(SolverCapabilities capabilities)
                : capabilities_m(capabilities) {}

            [[nodiscard]] SolverCapabilities capabilities() const override {
                return capabilities_m;
            }

            void execute(SolveContext& context, SelfFieldDiagnostics&) override {
                ++executionCount;
                selected.clear();
                for (const ParticleFieldBinding3d& binding : context.particles().bindings()) {
                    selected.push_back(binding.selectedForSolve);
                }
            }

            SolverCapabilities capabilities_m;
            std::size_t executionCount = 0;
            std::vector<bool> selected;
        };

        class SelfFieldSystemTest : public ::testing::Test {
        protected:
            using Container = ::ParticleContainer<double, 3>;

            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
            }

            static void TearDownTestSuite() { ippl::finalize(); }

            void SetUp() override {
                ippl::NDIndex<3> indexDomain;
                for (unsigned dimension = 0; dimension < 3; ++dimension) {
                    indexDomain[dimension] = ippl::Index(4);
                }
                mesh_m = std::make_unique<Mesh_t<3>>(
                        indexDomain, Vector_t<double, 3>(0.25), Vector_t<double, 3>(-0.5));
                layout_m = std::make_unique<FieldLayout_t<3>>(
                        MPI_COMM_WORLD, indexDomain, std::array<bool, 3>{false, false, false},
                        false);
                for (std::size_t index = 0; index < 3; ++index) {
                    auto container = std::make_shared<Container>(*mesh_m, *layout_m);
                    container->setBunchStateHandler(std::make_shared<BunchStateHandler>());
                    containers_m.push_back(std::move(container));
                }
            }

            void TearDown() override {
                containers_m.clear();
                layout_m.reset();
                mesh_m.reset();
            }

            [[nodiscard]] std::vector<ParticleFieldBinding3d> bindings() const {
                std::vector<ParticleFieldBinding3d> result;
                for (const auto& container : containers_m) {
                    result.push_back(bindParticleFields(*container));
                }
                return result;
            }

            [[nodiscard]] SelfFieldConfig pic3dConfig(
                    std::optional<BinningConfig> binning = std::nullopt,
                    CorrectionConfig correction          = CorrectionConfig()) const {
                Pic3DConfigValues values;
                values.backend    = PoissonBackendKind::None;
                values.meshSize   = {4, 4, 4};
                values.binning    = std::move(binning);
                values.correction = std::move(correction);
                ParticleStorageConfig3d storage;
                storage.meshSize = values.meshSize;
                return SelfFieldConfig(Pic3DConfig(std::move(values)), std::move(storage));
            }

            [[nodiscard]] SelfFieldConfig pic2d5Config() const {
                Pic2d5ConfigValues values;
                values.meshSize = {4, 4, 4};
                ParticleStorageConfig3d storage;
                storage.meshSize = values.meshSize;
                return SelfFieldConfig(Pic2d5Config(std::move(values)), std::move(storage));
            }

            std::unique_ptr<Mesh_t<3>> mesh_m;
            std::unique_ptr<FieldLayout_t<3>> layout_m;
            std::vector<std::shared_ptr<Container>> containers_m;
        };

        TEST_F(SelfFieldSystemTest, AcceptsExactBindingsAndRejectsCountOrderAndPrimaryMismatch) {
            auto expected  = bindings();
            auto algorithm = std::make_unique<RecordingAlgorithm>(SolverCapabilities{});
            RecordingAlgorithm* recorder = algorithm.get();
            SelfFieldSystem system(pic3dConfig(), std::move(algorithm), expected);

            auto exact = bindings();
            ParticleSetView exactView(exact, 0);
            SolveContext exactContext(exactView, {});
            EXPECT_NO_THROW(system.solve(exactContext));
            EXPECT_EQ(recorder->executionCount, 1u);

            auto wrongCount = bindings();
            wrongCount.pop_back();
            ParticleSetView wrongCountView(wrongCount, 0);
            SolveContext wrongCountContext(wrongCountView, {});
            EXPECT_THROW(system.solve(wrongCountContext), OpalException);

            auto wrongOrder = bindings();
            std::swap(wrongOrder[0], wrongOrder[1]);
            ParticleSetView wrongOrderView(wrongOrder, 0);
            SolveContext wrongOrderContext(wrongOrderView, {});
            EXPECT_THROW(system.solve(wrongOrderContext), OpalException);

            auto wrongPrimary = bindings();
            ParticleSetView wrongPrimaryView(wrongPrimary, 1);
            SolveContext wrongPrimaryContext(wrongPrimaryView, {});
            EXPECT_THROW(system.solve(wrongPrimaryContext), OpalException);
        }

        TEST_F(SelfFieldSystemTest, RejectsEveryMismatchedNativeIdentity) {
            const auto expected = bindings();
            auto algorithm      = std::make_unique<RecordingAlgorithm>(SolverCapabilities{});
            SelfFieldSystem system(pic3dConfig(), std::move(algorithm), expected);

            for (int mismatch = 0; mismatch < 5; ++mismatch) {
                auto actual = bindings();
                switch (mismatch) {
                    case 0:
                        actual[0].container = containers_m[1].get();
                        break;
                    case 1:
                        actual[0].position = &containers_m[1]->R;
                        break;
                    case 2:
                        actual[0].momentum = &containers_m[1]->P;
                        break;
                    case 3:
                        actual[0].electric = &containers_m[1]->E;
                        break;
                    case 4:
                        actual[0].magnetic = &containers_m[1]->B;
                        break;
                }
                ParticleSetView view(actual, 0);
                SolveContext context(view, {});
                EXPECT_THROW(system.solve(context), OpalException) << "mismatch index " << mismatch;
            }
        }

        TEST_F(SelfFieldSystemTest, AppliesPrimaryOnlySelectionAfterBindingValidation) {
            auto expected = bindings();
            SolverCapabilities capabilities;
            capabilities.particleSelection = ParticleSelectionPolicy::PrimaryOnly;
            auto algorithm                 = std::make_unique<RecordingAlgorithm>(capabilities);
            RecordingAlgorithm* recorder   = algorithm.get();
            SelfFieldSystem system(pic3dConfig(), std::move(algorithm), expected);

            auto actual              = bindings();
            actual[0].trackingActive = true;
            actual[1].trackingActive = true;
            actual[2].trackingActive = false;
            ParticleSetView view(actual, 0);
            SolveContext context(view, {});
            system.solve(context);

            ASSERT_EQ(recorder->selected.size(), 3u);
            EXPECT_TRUE(recorder->selected[0]);
            EXPECT_FALSE(recorder->selected[1]);
            EXPECT_FALSE(recorder->selected[2]);
        }

        TEST_F(SelfFieldSystemTest, AppliesAllTrackingActiveSelection) {
            auto expected = bindings();
            SolverCapabilities capabilities;
            capabilities.particleSelection = ParticleSelectionPolicy::AllTrackingActive;
            auto algorithm                 = std::make_unique<RecordingAlgorithm>(capabilities);
            RecordingAlgorithm* recorder   = algorithm.get();
            SelfFieldSystem system(pic2d5Config(), std::move(algorithm), expected);

            auto actual              = bindings();
            actual[0].trackingActive = true;
            actual[1].trackingActive = true;
            actual[2].trackingActive = false;
            ParticleSetView view(actual, 0);
            SolveContext context(view, {});
            system.solve(context);

            ASSERT_EQ(recorder->selected.size(), 3u);
            EXPECT_TRUE(recorder->selected[0]);
            EXPECT_TRUE(recorder->selected[1]);
            EXPECT_FALSE(recorder->selected[2]);
        }

        TEST_F(SelfFieldSystemTest, RejectsRequestsOutsideTheConfiguredStepPolicy) {
            BinningConfigValues binningValues;
            CorrectionConfig correction(
                    {.kind               = CorrectionKind::ImageCharge,
                     .planeZ             = 0.75,
                     .planeDumpFrequency = 2,
                     .maximumSteps       = 3});
            auto config = pic3dConfig(BinningConfig(std::move(binningValues)), correction);
            const SelfFieldRequestPolicy policy(config);
            auto expected = bindings();
            SolverCapabilities capabilities;
            capabilities.supportsBinning         = true;
            capabilities.supportsImageCharge     = true;
            capabilities.supportsPotentialOutput = true;
            auto algorithm               = std::make_unique<RecordingAlgorithm>(capabilities);
            RecordingAlgorithm* recorder = algorithm.get();
            SelfFieldSystem system(std::move(config), std::move(algorithm), expected);

            auto exact = bindings();
            ParticleSetView exactView(exact, 0);
            StepState exactStep;
            exactStep.step = 2;
            SolveContext exactContext(exactView, exactStep, policy.forStep(2));
            EXPECT_NO_THROW(system.solve(exactContext));
            EXPECT_EQ(recorder->executionCount, 1u);

            auto wrong = bindings();
            ParticleSetView wrongView(wrong, 0);
            RequestedPhysics wrongRequest = policy.forStep(2);
            wrongRequest.correction.planeZ += 1.0;
            SolveContext wrongContext(wrongView, exactStep, wrongRequest);
            EXPECT_THROW(system.solve(wrongContext), OpalException);

            auto expired = bindings();
            ParticleSetView expiredView(expired, 0);
            StepState expiredStep;
            expiredStep.step = 3;
            SolveContext expiredContext(expiredView, expiredStep, policy.forStep(2));
            EXPECT_THROW(system.solve(expiredContext), OpalException);
        }

        TEST_F(SelfFieldSystemTest, FrameGuardRestoresPrimaryPositionAndFieldsDuringUnwinding) {
            auto& particles = *containers_m.front();
            particles.createParticles(1);
            auto position = particles.R.getHostMirror();
            auto electric = particles.E.getHostMirror();
            auto magnetic = particles.B.getHostMirror();
            position(0)   = Vector_t<double, 3>(1.0, 2.0, 3.0);
            electric(0)   = Vector_t<double, 3>(4.0, 5.0, 6.0);
            magnetic(0)   = Vector_t<double, 3>(7.0, 8.0, 9.0);
            Kokkos::deep_copy(particles.R.getView(), position);
            Kokkos::deep_copy(particles.E.getView(), electric);
            Kokkos::deep_copy(particles.B.getView(), magnetic);

            const CoordinateSystemTrafo trackerToSolve(
                    Vector_t<double, 3>(0.25, -0.5, 1.0), Quaternion(1.0, 0.0, 0.0, 0.0));
            const FrameState frames{trackerToSolve, trackerToSolve.inverted()};
            const auto throwDuringSolve = [&] {
                ParticleFrameGuard<double, 3> guard(frames, particles);
                guard.enter();
                guard.markComputedFields();
                throw std::runtime_error("deliberate solve failure");
            };
            EXPECT_THROW(throwDuringSolve(), std::runtime_error);

            const auto restoredPosition =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.R.getView());
            const auto restoredElectric =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.E.getView());
            const auto restoredMagnetic =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.B.getView());
            for (unsigned dimension = 0; dimension < 3; ++dimension) {
                EXPECT_DOUBLE_EQ(restoredPosition(0)[dimension], position(0)[dimension]);
                EXPECT_DOUBLE_EQ(restoredElectric(0)[dimension], electric(0)[dimension]);
                EXPECT_DOUBLE_EQ(restoredMagnetic(0)[dimension], magnetic(0)[dimension]);
            }
        }

    }  // namespace
}  // namespace opalx::spacecharge
