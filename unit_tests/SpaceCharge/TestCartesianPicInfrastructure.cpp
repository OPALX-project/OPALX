#include <gtest/gtest.h>

#include "PartBunch/CartesianDomain.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"
#include "SpaceCharge/Ippl/IpplPoissonAdapter.h"
#include "SpaceCharge/Pic3D/PicWorkspace.h"
#include "Utilities/OpalException.h"

#include <array>
#include <memory>

namespace opalx::spacecharge {
    namespace {

        class CartesianPicInfrastructureTest : public ::testing::Test {
        protected:
            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
            }

            static void TearDownTestSuite() { ippl::finalize(); }

            [[nodiscard]] ParticleStorageConfig3d storage(bool periodic = false) const {
                ParticleStorageConfig3d result;
                result.meshSize                 = {8, 8, 8};
                result.decomposition            = {false, false, false};
                result.periodicParticleBoundary = periodic;
                return result;
            }

            void exerciseBackend(
                    PoissonBackendKind kind, BoundaryConditionKind boundary,
                    double cutoff = 0.0) const {
                auto setup = storage(boundary == BoundaryConditionKind::Periodic);
                CartesianDomain<double, 3> domain(setup);
                PicWorkspace<double, 3> workspace(domain);
                workspace.initializeFields("TEST");

                IpplPoissonBackendConfig config;
                config.kind               = kind;
                config.p3mCutoff          = cutoff;
                config.boundaryConditions = {boundary, boundary, boundary};
                IpplPoissonAdapter adapter(
                        config, {&workspace.chargeDensity(), &workspace.electricField()});
                EXPECT_NO_THROW(adapter.warmup());
                EXPECT_NO_THROW(adapter.solve({}, {.suppressFieldDump = true}));
            }
        };

        TEST_F(CartesianPicInfrastructureTest, DomainMutatesLayoutAndGeometryInPlace) {
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

        TEST_F(CartesianPicInfrastructureTest, OverlapDomainUsesCutoffAndConfiguredBoundaries) {
            auto periodicSetup          = storage(true);
            periodicSetup.meshSize      = {4, 4, 4};
            periodicSetup.layoutKind    = ParticleLayoutKind::SpatialOverlap;
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

        TEST_F(CartesianPicInfrastructureTest, WorkspaceRefreshesFieldsAfterDomainLayoutChange) {
            auto setup = storage(true);
            CartesianDomain<double, 3> domain(setup);
            PicWorkspace<double, 3> workspace(domain);
            workspace.initializeFields("FFT");
            IpplPoissonBackendConfig backendConfig;
            backendConfig.kind               = PoissonBackendKind::FftPeriodic;
            backendConfig.boundaryConditions = {
                    BoundaryConditionKind::Periodic, BoundaryConditionKind::Periodic,
                    BoundaryConditionKind::Periodic};
            IpplPoissonAdapter backend(
                    backendConfig, {&workspace.chargeDensity(), &workspace.electricField()});
            EXPECT_EQ(&workspace.mesh(), &domain.mesh());
            EXPECT_EQ(&workspace.layout(), &domain.layout());

            ASSERT_TRUE(domain.rebuildGlobalLayoutInPlace(
                    {8, 8, 16}, std::array<bool, 3>{false, false, true}));
            EXPECT_NO_THROW(workspace.updateFieldLayoutsAfterLayoutChange());
            EXPECT_NO_THROW(
                    backend.refresh({&workspace.chargeDensity(), &workspace.electricField()}));
            EXPECT_NO_THROW(backend.solve({}, {.suppressFieldDump = true}));
            EXPECT_EQ(workspace.layoutExtents(), (std::array<std::size_t, 3>{8, 8, 16}));
            EXPECT_EQ(&workspace.electricField().getLayout(), &domain.layout());
            EXPECT_EQ(&workspace.chargeDensity().getLayout(), &domain.layout());
        }

        TEST_F(CartesianPicInfrastructureTest, VariantDispatchesNullPeriodicAndOpenBackends) {
            exerciseBackend(PoissonBackendKind::None, BoundaryConditionKind::Open);
            exerciseBackend(PoissonBackendKind::FftPeriodic, BoundaryConditionKind::Periodic);
            exerciseBackend(PoissonBackendKind::Open, BoundaryConditionKind::Open);
        }

        TEST_F(CartesianPicInfrastructureTest, OpenBackendSupportsShiftedGreenSolve) {
            auto setup = storage();
            CartesianDomain<double, 3> domain(setup);
            PicWorkspace<double, 3> workspace(domain);
            workspace.initializeFields("OPEN");
            IpplPoissonBackendConfig config;
            config.kind = PoissonBackendKind::Open;
            IpplPoissonAdapter adapter(
                    config, {&workspace.chargeDensity(), &workspace.electricField()});
            adapter.warmup();
            IpplPoissonSolveRequest request;
            request.greenFunctionShift = Vector_t<double, 3>(0.0, 0.0, 0.25);
            EXPECT_NO_THROW(adapter.solve(request, {.suppressFieldDump = true}));
        }

        TEST_F(CartesianPicInfrastructureTest, VariantDispatchesOpenAndPeriodicP3MBackends) {
            exerciseBackend(PoissonBackendKind::P3M, BoundaryConditionKind::Open, 0.25);
            exerciseBackend(PoissonBackendKind::P3M, BoundaryConditionKind::Periodic, 0.25);
        }

        TEST_F(CartesianPicInfrastructureTest, CapabilityAndCouplingTablesRejectCG) {
            EXPECT_TRUE(IpplPoissonAdapter::capabilitiesFor(PoissonBackendKind::None).isNoOp);
            EXPECT_TRUE(
                    IpplPoissonAdapter::capabilitiesFor(PoissonBackendKind::Open)
                            .supportsShiftedGreenFunction);
            EXPECT_DOUBLE_EQ(
                    IpplPoissonAdapter::couplingConstantFor(PoissonBackendKind::FftPeriodic),
                    1.0 / Physics::epsilon_0);
            EXPECT_THROW(
                    static_cast<void>(IpplPoissonAdapter::capabilitiesFor(
                            PoissonBackendKind::ConjugateGradient)),
                    OpalException);
            EXPECT_THROW(
                    static_cast<void>(IpplPoissonAdapter::couplingConstantFor(
                            PoissonBackendKind::ConjugateGradient)),
                    OpalException);
        }

    }  // namespace
}  // namespace opalx::spacecharge
