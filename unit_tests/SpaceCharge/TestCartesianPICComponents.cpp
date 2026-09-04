#include <gtest/gtest.h>

#include "PartBunch/CartesianDomain.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"
#include "SpaceCharge/CartesianPIC/CartesianPICFieldStorage.h"
#include "SpaceCharge/Poisson/PoissonSolver.h"
#include "Utilities/OpalException.h"

#include <array>
#include <memory>

namespace opalx::spacecharge {
    namespace {

        class CartesianPICComponentsTest : public ::testing::Test {
        protected:
            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
            }

            static void TearDownTestSuite() { ippl::finalize(); }

            [[nodiscard]] CartesianDomainConfig3D storage(bool periodic = false) const {
                CartesianDomainConfig3D result;
                result.meshSize                 = {8, 8, 8};
                result.decomposition            = {false, false, false};
                result.periodicParticleBoundary = periodic;
                return result;
            }

            void exerciseBackend(
                    PoissonSolverType kind, FieldBoundaryCondition boundary,
                    double cutoff = 0.0) const {
                auto setup = storage(boundary == FieldBoundaryCondition::Periodic);
                CartesianDomain<double, 3> domain(setup);
                CartesianPICFieldStorage<double, 3> workspace(domain);
                workspace.initializeFields("TEST");

                PoissonSolverConfig config;
                config.type               = kind;
                config.p3mCutoff          = cutoff;
                config.boundaryConditions = {boundary, boundary, boundary};
                PoissonSolver adapter(
                        config, {&workspace.chargeDensity(), &workspace.electricField()});
                EXPECT_NO_THROW(adapter.warmup());
                EXPECT_NO_THROW(adapter.solve({}, {.suppressFieldDump = true}));
            }
        };

        TEST_F(CartesianPICComponentsTest, DomainMutatesLayoutAndGeometryInPlace) {
            auto setup     = storage();
            setup.meshSize = {4, 5, 6};
            CartesianDomain<double, 3> domain(setup);
            auto* meshAddress   = &domain.mesh();
            auto* layoutAddress = &domain.layout();
            EXPECT_EQ(domain.layoutExtents(), (std::array<std::size_t, 3>{4, 5, 6}));

            EXPECT_TRUE(domain.rebuildGlobalLayoutInPlace(
                    {4, 5, 12}, std::array<bool, 3>{false, false, true}));
            EXPECT_EQ(&domain.mesh(), meshAddress);
            EXPECT_EQ(&domain.layout(), layoutAddress);
            EXPECT_EQ(domain.layoutExtents(), (std::array<std::size_t, 3>{4, 5, 12}));
            EXPECT_FALSE(domain.rebuildGlobalLayoutInPlace(
                    {4, 5, 12}, std::array<bool, 3>{true, true, true}));

            const Vector_t<double, 3> lower(-1.0, -2.0, -3.0);
            const Vector_t<double, 3> upper(1.0, 2.0, 3.0);
            const Vector_t<double, 3> spacing(0.5, 1.0, 1.5);
            const Vector_t<double, 3> origin(-1.25, -2.5, -3.75);
            domain.setGeometry(lower, upper, spacing, origin);
            for (unsigned dimension = 0; dimension < 3; ++dimension) {
                EXPECT_DOUBLE_EQ(domain.lower()[dimension], lower[dimension]);
                EXPECT_DOUBLE_EQ(domain.upper()[dimension], upper[dimension]);
                EXPECT_DOUBLE_EQ(domain.mesh().getMeshSpacing()[dimension], spacing[dimension]);
                EXPECT_DOUBLE_EQ(domain.mesh().getOrigin()[dimension], origin[dimension]);
            }
        }

        TEST_F(CartesianPICComponentsTest, OverlapDomainUsesCutoffAndConfiguredBoundaries) {
            auto periodicSetup          = storage(true);
            periodicSetup.meshSize      = {4, 4, 4};
            periodicSetup.layoutType    = ParticleLayoutType::SpatialOverlap;
            periodicSetup.overlapCutoff = 0.25;
            CartesianDomain<double, 3> periodicDomain(periodicSetup);
            EXPECT_DOUBLE_EQ(periodicDomain.lower()[0], -0.5);
            EXPECT_DOUBLE_EQ(periodicDomain.upper()[0], 0.5);
            EXPECT_DOUBLE_EQ(periodicDomain.spacing()[0], 0.25);

            using Container = ::ParticleContainer<double, 3>;
            Container periodic(
                    periodicDomain.mesh(), periodicDomain.layout(), false,
                    Container::LayoutType::SpatialOverlap, 0.25, ippl::BC::PERIODIC);
            for (const ippl::BC boundary : periodic.getP3MLayout().getParticleBC()) {
                EXPECT_EQ(boundary, ippl::BC::PERIODIC);
            }

            auto openSetup                     = periodicSetup;
            openSetup.periodicParticleBoundary = false;
            CartesianDomain<double, 3> openDomain(openSetup);
            Container open(
                    openDomain.mesh(), openDomain.layout(), false,
                    Container::LayoutType::SpatialOverlap, 0.25, ippl::BC::NO);
            for (const ippl::BC boundary : open.getP3MLayout().getParticleBC()) {
                EXPECT_EQ(boundary, ippl::BC::NO);
            }
        }

        TEST_F(CartesianPICComponentsTest, FieldStorageRefreshesFieldsAfterDomainLayoutChange) {
            auto setup = storage(true);
            CartesianDomain<double, 3> domain(setup);
            CartesianPICFieldStorage<double, 3> workspace(domain);
            workspace.initializeFields("FFT");
            PoissonSolverConfig backendConfig;
            backendConfig.type               = PoissonSolverType::PeriodicFFT;
            backendConfig.boundaryConditions = {
                    FieldBoundaryCondition::Periodic, FieldBoundaryCondition::Periodic,
                    FieldBoundaryCondition::Periodic};
            PoissonSolver backend(
                    backendConfig, {&workspace.chargeDensity(), &workspace.electricField()});
            EXPECT_EQ(&workspace.mesh(), &domain.mesh());
            EXPECT_EQ(&workspace.layout(), &domain.layout());

            ASSERT_TRUE(domain.rebuildGlobalLayoutInPlace(
                    {8, 8, 16}, std::array<bool, 3>{false, false, true}));
            EXPECT_NO_THROW(workspace.updateFieldLayoutsAfterLayoutChange());
            EXPECT_NO_THROW(backend.rebuildAfterLayoutChange(
                    {&workspace.chargeDensity(), &workspace.electricField()}));
            EXPECT_NO_THROW(backend.solve({}, {.suppressFieldDump = true}));
            EXPECT_EQ(workspace.layoutExtents(), (std::array<std::size_t, 3>{8, 8, 16}));
            EXPECT_EQ(&workspace.electricField().getLayout(), &domain.layout());
            EXPECT_EQ(&workspace.chargeDensity().getLayout(), &domain.layout());

            ASSERT_TRUE(domain.rebuildGlobalLayoutInPlace(
                    {8, 8, 8}, std::array<bool, 3>{true, true, true}));
            EXPECT_NO_THROW(workspace.updateFieldLayoutsAfterLayoutChange());
            EXPECT_NO_THROW(backend.rebuildAfterLayoutChange(
                    {&workspace.chargeDensity(), &workspace.electricField()}));
            EXPECT_NO_THROW(backend.solve({}, {.suppressFieldDump = true}));
            EXPECT_EQ(workspace.layoutExtents(), (std::array<std::size_t, 3>{8, 8, 8}));
        }

        TEST_F(CartesianPICComponentsTest, VariantDispatchesNullPeriodicAndOpenBackends) {
            exerciseBackend(PoissonSolverType::None, FieldBoundaryCondition::Open);
            exerciseBackend(PoissonSolverType::PeriodicFFT, FieldBoundaryCondition::Periodic);
            exerciseBackend(PoissonSolverType::Open, FieldBoundaryCondition::Open);
        }

        TEST_F(CartesianPICComponentsTest, OpenBackendSupportsShiftedGreenSolve) {
            auto setup = storage();
            CartesianDomain<double, 3> domain(setup);
            CartesianPICFieldStorage<double, 3> workspace(domain);
            workspace.initializeFields("OPEN");
            PoissonSolverConfig config;
            config.type = PoissonSolverType::Open;
            PoissonSolver adapter(config, {&workspace.chargeDensity(), &workspace.electricField()});
            adapter.warmup();
            PoissonSolveRequest request;
            request.greenFunctionShift = Vector_t<double, 3>(0.0, 0.0, 0.25);
            EXPECT_NO_THROW(adapter.solve(request, {.suppressFieldDump = true}));
        }

        TEST_F(CartesianPICComponentsTest, VariantDispatchesOpenAndPeriodicP3MBackends) {
            exerciseBackend(PoissonSolverType::P3M, FieldBoundaryCondition::Open, 0.25);
            exerciseBackend(PoissonSolverType::P3M, FieldBoundaryCondition::Periodic, 0.25);
        }

        TEST_F(CartesianPICComponentsTest, CapabilityAndCouplingTablesRejectCG) {
            EXPECT_TRUE(PoissonSolver::capabilitiesFor(PoissonSolverType::None).isNoOp);
            EXPECT_TRUE(
                    PoissonSolver::capabilitiesFor(PoissonSolverType::Open)
                            .supportsShiftedGreenFunction);
            EXPECT_DOUBLE_EQ(
                    PoissonSolver::couplingConstantFor(PoissonSolverType::PeriodicFFT),
                    1.0 / Physics::epsilon_0);
            EXPECT_THROW(
                    static_cast<void>(
                            PoissonSolver::capabilitiesFor(PoissonSolverType::ConjugateGradient)),
                    OpalException);
            EXPECT_THROW(
                    static_cast<void>(PoissonSolver::couplingConstantFor(
                            PoissonSolverType::ConjugateGradient)),
                    OpalException);
        }

    }  // namespace
}  // namespace opalx::spacecharge
