#include <gtest/gtest.h>

#include "Algorithms/Quaternion.hpp"
#include "PartBunch/BunchStateHandler.h"
#include "SpaceCharge/SpaceChargeFrameTransform.h"
#include "SpaceCharge/SpaceChargeRequestSchedule.h"
#include "SpaceCharge/SpaceChargeSolver.h"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <array>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace opalx::spacecharge {
    namespace {

        class RecordingAlgorithm final : public SpaceChargeAlgorithm {
        public:
            explicit RecordingAlgorithm(SpaceChargeCapabilities capabilities)
                : capabilities_m(capabilities) {}

            [[nodiscard]] SpaceChargeCapabilities capabilities() const override {
                return capabilities_m;
            }

            void solve(SpaceChargeSolveContext& context, SpaceChargeDiagnostics&) override {
                ++executionCount;
                selected.clear();
                for (const ParticleFieldBinding3D& binding : context.particles().bindings()) {
                    selected.push_back(binding.selectedForSolve);
                }
            }

            SpaceChargeCapabilities capabilities_m;
            std::size_t executionCount = 0;
            std::vector<bool> selected;
        };

        class ThrowingAlgorithm final : public SpaceChargeAlgorithm {
        public:
            [[nodiscard]] SpaceChargeCapabilities capabilities() const override { return {}; }

            void solve(SpaceChargeSolveContext&, SpaceChargeDiagnostics&) override {
                throw std::runtime_error("deliberate algorithm failure");
            }
        };

        class SpaceChargeSolverTest : public ::testing::Test {
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
                    container->setBunchStateHandler(bunchState_m);
                    containers_m.push_back(std::move(container));
                }
            }

            void TearDown() override {
                containers_m.clear();
                layout_m.reset();
                mesh_m.reset();
            }

            [[nodiscard]] std::vector<ParticleFieldBinding3D> bindings() const {
                std::vector<ParticleFieldBinding3D> result;
                for (const auto& container : containers_m) {
                    result.push_back(makeParticleFieldBinding(*container));
                }
                return result;
            }

            [[nodiscard]] SpaceChargeConfig cartesianPICConfig(
                    std::optional<BinningConfig> binning = std::nullopt,
                    CorrectionConfig correction          = CorrectionConfig()) const {
                CartesianPICConfig::Parameters values;
                values.backend    = PoissonSolverType::None;
                values.meshSize   = {4, 4, 4};
                values.binning    = std::move(binning);
                values.correction = std::move(correction);
                return cartesianPICConfig(std::move(values));
            }

            [[nodiscard]] SpaceChargeConfig cartesianPICConfig(
                    CartesianPICConfig::Parameters values) const {
                CartesianDomainConfig3D storage;
                storage.meshSize                   = values.meshSize;
                storage.decomposition              = values.parallelDimensions;
                storage.boundingBoxIncreasePercent = values.boundingBoxIncreasePercent;
                const bool periodic                = std::all_of(
                        values.boundaryConditions.begin(), values.boundaryConditions.end(),
                        [](FieldBoundaryCondition boundary) {
                            return boundary == FieldBoundaryCondition::Periodic;
                        });
                storage.periodicParticleBoundary = periodic;
                if (values.backend == PoissonSolverType::P3M) {
                    storage.layoutType    = ParticleLayoutType::SpatialOverlap;
                    storage.overlapCutoff = values.p3mCutoff;
                }
                return SpaceChargeConfig(CartesianPICConfig(std::move(values)), std::move(storage));
            }

            [[nodiscard]] SpaceChargeConfig fft2D5Config() const {
                FFT2D5Config::Parameters values;
                values.meshSize = {4, 4, 4};
                CartesianDomainConfig3D storage;
                storage.meshSize = values.meshSize;
                return SpaceChargeConfig(FFT2D5Config(std::move(values)), std::move(storage));
            }

            std::unique_ptr<Mesh_t<3>> mesh_m;
            std::unique_ptr<FieldLayout_t<3>> layout_m;
            std::shared_ptr<BunchStateHandler> bunchState_m = std::make_shared<BunchStateHandler>();
            std::vector<std::shared_ptr<Container>> containers_m;
        };

        TEST_F(SpaceChargeSolverTest, AcceptsExactBindingsAndRejectsCountOrderAndPrimaryMismatch) {
            auto expected  = bindings();
            auto algorithm = std::make_unique<RecordingAlgorithm>(SpaceChargeCapabilities{});
            RecordingAlgorithm* recorder = algorithm.get();
            SpaceChargeSolver system(
                    cartesianPICConfig(), std::move(algorithm), expected, bunchState_m);

            auto exact = bindings();
            ParticleFieldSet exactView(exact, 0);
            SpaceChargeSolveContext exactContext(exactView, {});
            EXPECT_NO_THROW(system.solve(exactContext));
            EXPECT_EQ(recorder->executionCount, 1u);

            auto wrongCount = bindings();
            wrongCount.pop_back();
            ParticleFieldSet wrongCountView(wrongCount, 0);
            SpaceChargeSolveContext wrongCountContext(wrongCountView, {});
            EXPECT_THROW(system.solve(wrongCountContext), OpalException);

            auto wrongOrder = bindings();
            std::swap(wrongOrder[0], wrongOrder[1]);
            ParticleFieldSet wrongOrderView(wrongOrder, 0);
            SpaceChargeSolveContext wrongOrderContext(wrongOrderView, {});
            EXPECT_THROW(system.solve(wrongOrderContext), OpalException);

            auto wrongPrimary = bindings();
            ParticleFieldSet wrongPrimaryView(wrongPrimary, 1);
            SpaceChargeSolveContext wrongPrimaryContext(wrongPrimaryView, {});
            EXPECT_THROW(system.solve(wrongPrimaryContext), OpalException);
        }

        TEST_F(SpaceChargeSolverTest, RejectsEveryMismatchedNativeIdentity) {
            const auto expected = bindings();
            auto algorithm      = std::make_unique<RecordingAlgorithm>(SpaceChargeCapabilities{});
            SpaceChargeSolver system(
                    cartesianPICConfig(), std::move(algorithm), expected, bunchState_m);

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
                ParticleFieldSet view(actual, 0);
                SpaceChargeSolveContext context(view, {});
                EXPECT_THROW(system.solve(context), OpalException) << "mismatch index " << mismatch;
            }
        }

        TEST_F(SpaceChargeSolverTest, AppliesPrimaryOnlySelectionAfterBindingValidation) {
            auto expected = bindings();
            SpaceChargeCapabilities capabilities;
            capabilities.particleSelection = ParticleSelectionMode::PrimaryOnly;
            auto algorithm                 = std::make_unique<RecordingAlgorithm>(capabilities);
            RecordingAlgorithm* recorder   = algorithm.get();
            SpaceChargeSolver system(
                    cartesianPICConfig(), std::move(algorithm), expected, bunchState_m);

            auto actual              = bindings();
            actual[0].trackingActive = true;
            actual[1].trackingActive = true;
            actual[2].trackingActive = false;
            ParticleFieldSet view(actual, 0);
            SpaceChargeSolveContext context(view, {});
            system.solve(context);

            ASSERT_EQ(recorder->selected.size(), 3u);
            EXPECT_TRUE(recorder->selected[0]);
            EXPECT_FALSE(recorder->selected[1]);
            EXPECT_FALSE(recorder->selected[2]);
        }

        TEST_F(SpaceChargeSolverTest, AppliesAllTrackingActiveSelection) {
            auto expected = bindings();
            SpaceChargeCapabilities capabilities;
            capabilities.particleSelection = ParticleSelectionMode::AllTrackingActive;
            auto algorithm                 = std::make_unique<RecordingAlgorithm>(capabilities);
            RecordingAlgorithm* recorder   = algorithm.get();
            SpaceChargeSolver system(fft2D5Config(), std::move(algorithm), expected, bunchState_m);

            auto actual              = bindings();
            actual[0].trackingActive = true;
            actual[1].trackingActive = true;
            actual[2].trackingActive = false;
            ParticleFieldSet view(actual, 0);
            SpaceChargeSolveContext context(view, {});
            system.solve(context);

            ASSERT_EQ(recorder->selected.size(), 3u);
            EXPECT_TRUE(recorder->selected[0]);
            EXPECT_TRUE(recorder->selected[1]);
            EXPECT_FALSE(recorder->selected[2]);
        }

        TEST_F(SpaceChargeSolverTest, RejectsRequestsOutsideTheConfiguredSchedule) {
            BinningConfig::Parameters binningValues;
            CorrectionConfig correction(
                    {.kind               = SpaceChargeCorrectionType::ImageCharge,
                     .planeZ             = 0.75,
                     .planeDumpFrequency = 2,
                     .maximumSteps       = 3});
            auto config = cartesianPICConfig(BinningConfig(std::move(binningValues)), correction);
            const SpaceChargeRequestSchedule schedule(config);
            auto expected = bindings();
            SpaceChargeCapabilities capabilities;
            capabilities.supportsBinning         = true;
            capabilities.supportsImageCharge     = true;
            capabilities.supportsPotentialOutput = true;
            auto algorithm               = std::make_unique<RecordingAlgorithm>(capabilities);
            RecordingAlgorithm* recorder = algorithm.get();
            SpaceChargeSolver system(
                    std::move(config), std::move(algorithm), expected, bunchState_m);

            auto exact = bindings();
            ParticleFieldSet exactView(exact, 0);
            SpaceChargeStepState exactStep;
            exactStep.step = 2;
            SpaceChargeSolveContext exactContext(exactView, exactStep, schedule.requestForStep(2));
            EXPECT_NO_THROW(system.solve(exactContext));
            EXPECT_EQ(recorder->executionCount, 1u);

            auto wrong = bindings();
            ParticleFieldSet wrongView(wrong, 0);
            SpaceChargeRequest wrongRequest = schedule.requestForStep(2);
            wrongRequest.correction.planeZ += 1.0;
            SpaceChargeSolveContext wrongContext(wrongView, exactStep, wrongRequest);
            EXPECT_THROW(system.solve(wrongContext), OpalException);

            auto expired = bindings();
            ParticleFieldSet expiredView(expired, 0);
            SpaceChargeStepState expiredStep;
            expiredStep.step = 3;
            SpaceChargeSolveContext expiredContext(
                    expiredView, expiredStep, schedule.requestForStep(2));
            EXPECT_THROW(system.solve(expiredContext), OpalException);
        }

        TEST_F(SpaceChargeSolverTest, FrameTransformRestoresPrimaryPositionAndFieldsExplicitly) {
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
            const CoordinateFrameTransforms frames{trackerToSolve, trackerToSolve.inverted()};
            SpaceChargeFrameTransform<double, 3> transform(frames, particles);
            transform.enter();
            transform.markComputedFields();
            transform.leave();

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

        TEST_F(SpaceChargeSolverTest, PropagatesAlgorithmFailureWithoutRecoveryContract) {
            auto expected  = bindings();
            auto algorithm = std::make_unique<ThrowingAlgorithm>();
            SpaceChargeSolver system(
                    cartesianPICConfig(), std::move(algorithm), expected, bunchState_m);
            auto actual = bindings();
            ParticleFieldSet view(actual, 0);
            SpaceChargeSolveContext context(view, {});

            EXPECT_THROW(system.solve(context), std::runtime_error);
        }

        TEST_F(SpaceChargeSolverTest, AcceptsFixedDomainForOpenWithoutCorrectionsOrORB) {
            CartesianPICConfig::Parameters values;
            values.backend  = PoissonSolverType::Open;
            values.meshSize = {4, 4, 4};
            auto expected   = bindings();
            SpaceChargeCapabilities capabilities;
            capabilities.supportsFixedCartesianDomain = true;
            auto algorithm               = std::make_unique<RecordingAlgorithm>(capabilities);
            RecordingAlgorithm* recorder = algorithm.get();
            SpaceChargeSolver system(
                    cartesianPICConfig(std::move(values)), std::move(algorithm), expected,
                    bunchState_m);
            bunchState_m->setFixedCartesianDomain({-1.0, -2.0, -3.0}, {1.0, 2.0, 3.0});

            auto actual = bindings();
            ParticleFieldSet view(actual, 0);
            SpaceChargeSolveContext context(view, {});
            EXPECT_NO_THROW(system.solve(context));
            EXPECT_EQ(recorder->executionCount, 1u);
        }

        TEST_F(SpaceChargeSolverTest, RejectsFixedDomainForUnsupportedAlgorithmAndBackends) {
            bunchState_m->setFixedCartesianDomain({-1.0, -2.0, -3.0}, {1.0, 2.0, 3.0});

            {
                auto expected  = bindings();
                auto algorithm = std::make_unique<RecordingAlgorithm>(SpaceChargeCapabilities{});
                SpaceChargeSolver system(
                        fft2D5Config(), std::move(algorithm), expected, bunchState_m);
                auto actual = bindings();
                ParticleFieldSet view(actual, 0);
                SpaceChargeSolveContext context(view, {});
                EXPECT_THROW(system.solve(context), OpalException);
            }

            for (const PoissonSolverType backend :
                 {PoissonSolverType::None, PoissonSolverType::PeriodicFFT,
                  PoissonSolverType::P3M}) {
                CartesianPICConfig::Parameters values;
                values.backend  = backend;
                values.meshSize = {4, 4, 4};
                if (backend == PoissonSolverType::PeriodicFFT) {
                    values.boundaryConditions = {
                            FieldBoundaryCondition::Periodic, FieldBoundaryCondition::Periodic,
                            FieldBoundaryCondition::Periodic};
                } else if (backend == PoissonSolverType::P3M) {
                    values.p3mCutoff = 0.25;
                }
                auto expected = bindings();
                SpaceChargeCapabilities capabilities;
                capabilities.supportsFixedCartesianDomain = true;
                auto algorithm = std::make_unique<RecordingAlgorithm>(capabilities);
                SpaceChargeSolver system(
                        cartesianPICConfig(std::move(values)), std::move(algorithm), expected,
                        bunchState_m);
                auto actual = bindings();
                ParticleFieldSet view(actual, 0);
                SpaceChargeSolveContext context(view, {});
                EXPECT_THROW(system.solve(context), OpalException);
            }
        }

        TEST_F(SpaceChargeSolverTest, RejectsCorrectionsAndORBWhileFixedDomainIsActive) {
            bunchState_m->setFixedCartesianDomain({-1.0, -2.0, -3.0}, {1.0, 2.0, 3.0});
            SpaceChargeCapabilities capabilities;
            capabilities.supportsFixedCartesianDomain = true;
            capabilities.supportsImageCharge          = true;
            capabilities.supportsRedistribution       = true;

            CartesianPICConfig::Parameters correctedValues;
            correctedValues.backend    = PoissonSolverType::Open;
            correctedValues.meshSize   = {4, 4, 4};
            correctedValues.correction = CorrectionConfig(
                    {.kind = SpaceChargeCorrectionType::ImageCharge, .planeZ = 0.5});
            auto correctedConfig = cartesianPICConfig(std::move(correctedValues));
            const SpaceChargeRequest correctedRequest =
                    SpaceChargeRequestSchedule(correctedConfig).requestForStep(0);
            auto correctedBindings  = bindings();
            auto correctedAlgorithm = std::make_unique<RecordingAlgorithm>(capabilities);
            SpaceChargeSolver corrected(
                    std::move(correctedConfig), std::move(correctedAlgorithm), correctedBindings,
                    bunchState_m);
            auto correctedActual = bindings();
            ParticleFieldSet correctedView(correctedActual, 0);
            SpaceChargeSolveContext correctedContext(correctedView, {}, correctedRequest);
            EXPECT_THROW(corrected.solve(correctedContext), OpalException);

            CartesianPICConfig::Parameters orbValues;
            orbValues.backend              = PoissonSolverType::Open;
            orbValues.meshSize             = {4, 4, 4};
            orbValues.repartitionFrequency = 1;
            auto orbBindings               = bindings();
            auto orbAlgorithm              = std::make_unique<RecordingAlgorithm>(capabilities);
            SpaceChargeSolver orb(
                    cartesianPICConfig(std::move(orbValues)), std::move(orbAlgorithm), orbBindings,
                    bunchState_m);
            auto orbActual = bindings();
            ParticleFieldSet orbView(orbActual, 0);
            SpaceChargeSolveContext orbContext(orbView, {});
            EXPECT_THROW(orb.solve(orbContext), OpalException);
        }

    }  // namespace
}  // namespace opalx::spacecharge
