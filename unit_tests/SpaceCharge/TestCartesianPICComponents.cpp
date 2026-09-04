#include <gtest/gtest.h>

#include "AbstractObjects/OpalData.h"
#include "Algorithms/Quaternion.hpp"
#include "PartBunch/BunchStateHandler.h"
#include "PartBunch/CartesianDomain.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"
#include "SpaceCharge/CartesianPIC/CartesianDomainUpdater.h"
#include "SpaceCharge/CartesianPIC/CartesianPICAlgorithm.h"
#include "SpaceCharge/CartesianPIC/CartesianPICFieldStorage.h"
#include "SpaceCharge/CartesianPIC/P3MShortRangeInteraction.h"
#include "SpaceCharge/CartesianPIC/ParticleDomainOperations.h"
#include "SpaceCharge/CartesianPIC/ParticleMeshFieldTransfer.h"
#include "SpaceCharge/CartesianPIC/RelativisticFieldComposer.h"
#include "SpaceCharge/Poisson/PoissonSolver.h"
#include "SpaceCharge/SpaceChargeSolveContext.h"
#include "Structure/DataSink.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

namespace opalx::spacecharge {
    namespace {

        class CartesianPICComponentsTest : public ::testing::Test {
        protected:
            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
                OpalData::getInstance()->storeInputFn("cartesian_pic_components.opal");
                gmsg                = new Inform(nullptr, -1);
                Options::enableHDF5 = false;
            }

            static void TearDownTestSuite() {
                delete gmsg;
                gmsg = nullptr;
                std::remove("cartesian_pic_components.stat");
                std::remove("cartesian_pic_components.lbal");
                ippl::finalize();
            }

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

            const auto fillCharge = [&workspace] {
                auto charge = workspace.chargeDensity().getHostMirror();
                for (std::size_t i = 0; i < charge.extent(0); ++i) {
                    for (std::size_t j = 0; j < charge.extent(1); ++j) {
                        for (std::size_t k = 0; k < charge.extent(2); ++k) {
                            charge(i, j, k) = static_cast<double>((i + 1) * (j + 2) + 3 * k);
                        }
                    }
                }
                Kokkos::deep_copy(workspace.chargeDensity().getView(), charge);
            };

            fillCharge();
            adapter.solve({}, {.suppressFieldDump = true});
            const auto standardField = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), workspace.electricField().getView());

            fillCharge();
            PoissonSolveRequest request;
            request.greenFunctionShift = Vector_t<double, 3>(0.0, 0.0, 0.25);
            EXPECT_NO_THROW(adapter.solve(request, {.suppressFieldDump = true}));

            fillCharge();
            EXPECT_NO_THROW(adapter.solve({}, {.suppressFieldDump = true}));
            const auto restoredStandardField = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), workspace.electricField().getView());
            for (std::size_t i = 0; i < standardField.extent(0); ++i) {
                for (std::size_t j = 0; j < standardField.extent(1); ++j) {
                    for (std::size_t k = 0; k < standardField.extent(2); ++k) {
                        for (unsigned dimension = 0; dimension < 3; ++dimension) {
                            const double expected = standardField(i, j, k)[dimension];
                            EXPECT_NEAR(
                                    restoredStandardField(i, j, k)[dimension], expected,
                                    1.0e-12 * std::max(1.0, std::abs(expected)));
                        }
                    }
                }
            }
        }

        TEST_F(CartesianPICComponentsTest, VariantDispatchesOpenAndPeriodicP3MBackends) {
            exerciseBackend(PoissonSolverType::P3M, FieldBoundaryCondition::Open, 0.25);
            exerciseBackend(PoissonSolverType::P3M, FieldBoundaryCondition::Periodic, 0.25);
        }

        TEST_F(CartesianPICComponentsTest, AdaptersExposeMetadataAndRejectReservedCG) {
            auto setup = storage();
            CartesianDomain<double, 3> domain(setup);
            CartesianPICFieldStorage<double, 3> workspace(domain);
            workspace.initializeFields("TEST");
            const PoissonFieldBinding fields{
                    &workspace.chargeDensity(), &workspace.electricField()};

            PoissonSolverConfig config;
            config.type = PoissonSolverType::None;
            PoissonSolver none(config, fields);
            EXPECT_EQ(none.name(), "NONE");
            EXPECT_TRUE(none.capabilities().isNoOp);
            EXPECT_DOUBLE_EQ(
                    none.couplingConstant(), 1.0 / (4.0 * Physics::pi * Physics::epsilon_0));

            config.type = PoissonSolverType::Open;
            PoissonSolver open(config, fields);
            EXPECT_EQ(open.name(), "OPEN");
            EXPECT_TRUE(open.capabilities().supportsShiftedGreenFunction);
            EXPECT_DOUBLE_EQ(open.couplingConstant(), 1.0 / Physics::epsilon_0);

            config.type = PoissonSolverType::ConjugateGradient;
            EXPECT_THROW(PoissonSolver(config, fields), OpalException);

            config.type = static_cast<PoissonSolverType>(255);
            EXPECT_THROW(PoissonSolver(config, fields), OpalException);
        }

        TEST_F(CartesianPICComponentsTest, NonOpenAdapterRejectsShiftedGreenRequest) {
            auto setup = storage();
            CartesianDomain<double, 3> domain(setup);
            CartesianPICFieldStorage<double, 3> workspace(domain);
            workspace.initializeFields("NONE");
            PoissonSolverConfig config;
            config.type = PoissonSolverType::None;
            PoissonSolver adapter(config, {&workspace.chargeDensity(), &workspace.electricField()});
            PoissonSolveRequest request;
            request.greenFunctionShift = Vector_t<double, 3>(0.0, 0.0, 0.25);
            EXPECT_THROW(adapter.solve(request, {.suppressFieldDump = true}), OpalException);
        }

        TEST_F(CartesianPICComponentsTest, FixedDomainOverridesStretchingAndReturnsToFollowing) {
            auto domainConfig          = storage();
            domainConfig.decomposition = {false, false, false};
            CartesianDomain<double, 3> domain(domainConfig);
            CartesianPICFieldStorage<double, 3> workspace(domain);
            workspace.initializeFields("NONE");

            auto bunchState = std::make_shared<BunchStateHandler>();
            using Container = ::ParticleContainer<double, 3>;
            auto particles  = std::make_shared<Container>(domain.mesh(), domain.layout());
            particles->setBunchStateHandler(bunchState);
            auto secondary = std::make_shared<Container>(domain.mesh(), domain.layout());
            secondary->setBunchStateHandler(bunchState);
            secondary->createParticles(1);
            secondary->setM(Physics::m_e);
            auto secondaryPositions = secondary->R.getHostMirror();
            auto secondaryMomenta   = secondary->P.getHostMirror();
            secondaryPositions(0)   = Vector_t<double, 3>(0.0);
            secondaryMomenta(0)     = Vector_t<double, 3>(0.0);
            Kokkos::deep_copy(secondary->R.getView(), secondaryPositions);
            Kokkos::deep_copy(secondary->P.getView(), secondaryMomenta);
            secondary->updateMoments();
            ASSERT_FALSE(secondary->isMomentsDirty());
            particles->createParticles(2);
            auto positions = particles->R.getHostMirror();
            auto momenta   = particles->P.getHostMirror();
            positions(0)   = Vector_t<double, 3>(-0.5, -0.25, -0.125);
            positions(1)   = Vector_t<double, 3>(0.5, 0.25, 0.125);
            momenta(0)     = Vector_t<double, 3>(0.0);
            momenta(1)     = Vector_t<double, 3>(0.0);
            Kokkos::deep_copy(particles->R.getView(), positions);
            Kokkos::deep_copy(particles->P.getView(), momenta);
            particles->setM(Physics::m_e);

            std::vector<ParticleFieldBinding3D> bindings{
                    makeParticleFieldBinding(*particles), makeParticleFieldBinding(*secondary)};
            ParticleFieldSet particleSet(bindings, 0);
            ParticleDomainOperations particleOperations(bindings);

            PoissonSolverConfig poissonConfig;
            poissonConfig.type = PoissonSolverType::None;
            PoissonSolver poisson(
                    poissonConfig, {&workspace.chargeDensity(), &workspace.electricField()});

            CartesianPICConfig::Parameters values;
            values.backend                         = PoissonSolverType::None;
            values.meshSize                        = domainConfig.meshSize;
            values.parallelDimensions              = domainConfig.decomposition;
            values.layoutRebuildParallelDimensions = std::array<bool, 3>{true, true, true};
            values.boundingBoxIncreasePercent      = 50.0;
            values.repartitionFrequency            = 1;
            CartesianDomainUpdater updater(CartesianPICConfig(values), bunchState);

            const std::array<double, 3> fixedLower{-2.0, -3.0, -4.0};
            const std::array<double, 3> fixedUpper{2.0, 3.0, 4.0};
            bunchState->setFixedCartesianDomain(fixedLower, fixedUpper);
            SpaceChargeStepState step;
            step.emissionActive  = true;
            step.emittedFraction = 0.01;
            SpaceChargeSolveContext context(particleSet, step);
            SpaceChargeDiagnostics diagnostics;

            updater.updateForSolve(
                    DomainCoordinateFrame::Beam, context, workspace, particleOperations, poisson,
                    diagnostics);
            EXPECT_EQ(workspace.layoutExtents(), domainConfig.meshSize);
            EXPECT_EQ(domain.decomposition(), domainConfig.decomposition);
            EXPECT_FALSE(particles->isMomentsDirty());
            EXPECT_FALSE(secondary->isMomentsDirty());
            EXPECT_EQ(diagnostics.redistributionCount(), 0u);
            for (unsigned dimension = 0; dimension < 3; ++dimension) {
                EXPECT_DOUBLE_EQ(domain.lower()[dimension], fixedLower[dimension]);
                EXPECT_DOUBLE_EQ(domain.upper()[dimension], fixedUpper[dimension]);
            }

            bunchState->clearFixedCartesianDomain();
            updater.updateForSolve(
                    DomainCoordinateFrame::Beam, context, workspace, particleOperations, poisson,
                    diagnostics);
            EXPECT_NE(domain.lower()[0], fixedLower[0]);
            EXPECT_NE(domain.upper()[0], fixedUpper[0]);
            EXPECT_FALSE(particles->isMomentsDirty());
            EXPECT_FALSE(secondary->isMomentsDirty());
        }

        TEST_F(CartesianPICComponentsTest, CartesianAlgorithmRetainsAndClearsFixedMesh) {
            auto domainConfig          = storage();
            domainConfig.decomposition = {false, false, false};
            CartesianDomain<double, 3> domain(domainConfig);
            auto fieldStorage   = std::make_unique<CartesianPICFieldStorage<double, 3>>(domain);
            auto* fieldObserver = fieldStorage.get();

            auto bunchState = std::make_shared<BunchStateHandler>();
            using Container = ::ParticleContainer<double, 3>;
            auto particles  = std::make_shared<Container>(domain.mesh(), domain.layout());
            particles->setBunchStateHandler(bunchState);
            particles->createParticles(2);
            particles->setQ(1.0e-15);
            particles->setM(Physics::m_e);
            auto positions = particles->R.getHostMirror();
            auto momenta   = particles->P.getHostMirror();
            auto timeSteps = particles->dt.getHostMirror();
            positions(0)   = Vector_t<double, 3>(-0.25, -0.1, -0.05);
            positions(1)   = Vector_t<double, 3>(0.25, 0.1, 0.05);
            momenta(0)     = Vector_t<double, 3>(0.0);
            momenta(1)     = Vector_t<double, 3>(0.0);
            timeSteps(0)   = 1.0e-12;
            timeSteps(1)   = 1.0e-12;
            Kokkos::deep_copy(particles->R.getView(), positions);
            Kokkos::deep_copy(particles->P.getView(), momenta);
            Kokkos::deep_copy(particles->dt.getView(), timeSteps);
            const auto originalPositions = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), particles->R.getView());

            std::vector<ParticleFieldBinding3D> bindings{makeParticleFieldBinding(*particles)};
            ParticleFieldSet particleSet(bindings, 0);
            CartesianPICConfig::Parameters values;
            values.backend                    = PoissonSolverType::Open;
            values.meshSize                   = domainConfig.meshSize;
            values.parallelDimensions         = domainConfig.decomposition;
            values.boundingBoxIncreasePercent = 10.0;
            DataSink dataSink;
            CartesianPICAlgorithm algorithm(
                    CartesianPICConfig(values), bindings, std::move(fieldStorage), &dataSink,
                    bunchState);

            const std::array<double, 3> fixedLower{-1.0, -1.5, -2.0};
            const std::array<double, 3> fixedUpper{1.0, 1.5, 2.0};
            bunchState->setFixedCartesianDomain(fixedLower, fixedUpper);
            SpaceChargeStepState step;
            step.timeStep = 1.0e-12;
            const CoordinateSystemTrafo trackerToSolve(
                    Vector_t<double, 3>(0.1, -0.2, 0.3),
                    Quaternion(std::cos(Physics::pi / 4.0), 0.0, 0.0, std::sin(Physics::pi / 4.0)));
            step.frames = {trackerToSolve, trackerToSolve.inverted()};
            SpaceChargeSolveContext fixedContext(particleSet, step);
            SpaceChargeDiagnostics diagnostics;
            algorithm.solve(fixedContext, diagnostics);
            for (unsigned dimension = 0; dimension < 3; ++dimension) {
                EXPECT_DOUBLE_EQ(fieldObserver->lower()[dimension], fixedLower[dimension]);
                EXPECT_DOUBLE_EQ(fieldObserver->upper()[dimension], fixedUpper[dimension]);
            }
            EXPECT_FALSE(particles->isMomentsDirty());
            const auto restoredPositions = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), particles->R.getView());
            for (std::size_t index = 0; index < particles->getLocalNum(); ++index) {
                for (unsigned dimension = 0; dimension < 3; ++dimension) {
                    EXPECT_NEAR(
                            restoredPositions(index)[dimension],
                            originalPositions(index)[dimension], 1.0e-14);
                }
            }

            bunchState->clearFixedCartesianDomain();
            SpaceChargeSolveContext followingContext(particleSet, step);
            algorithm.solve(followingContext, diagnostics);
            EXPECT_NE(fieldObserver->lower()[0], fixedLower[0]);
            EXPECT_NE(fieldObserver->upper()[0], fixedUpper[0]);
            EXPECT_FALSE(particles->isMomentsDirty());
        }

        TEST_F(CartesianPICComponentsTest, FieldTransferRestoresImageStateAndGathersVectorField) {
            Options::useQMAttributes = false;
            auto setup               = storage();
            CartesianDomain<double, 3> domain(setup);
            CartesianPICFieldStorage<double, 3> workspace(domain);
            workspace.initializeFields("NONE");

            using Container = ::ParticleContainer<double, 3>;
            Container particles(domain.mesh(), domain.layout());
            particles.setBunchStateHandler(std::make_shared<BunchStateHandler>());
            particles.createParticles(2);
            particles.setQ(2.0e-15);
            auto positions = particles.R.getHostMirror();
            auto timeSteps = particles.dt.getHostMirror();
            positions(0)   = Vector_t<double, 3>(-0.2, 0.1, 0.3);
            positions(1)   = Vector_t<double, 3>(0.25, -0.15, 0.45);
            timeSteps(0)   = 1.0e-12;
            timeSteps(1)   = 1.0e-12;
            Kokkos::deep_copy(particles.R.getView(), positions);
            Kokkos::deep_copy(particles.dt.getView(), timeSteps);
            const auto originalCharge =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.getQView());

            ParticleMeshFieldTransfer<double, 3> transfer;
            ParticleMeshFieldTransfer<double, 3>::ChargeNormalization normalization;
            normalization.timeStep       = 1.0e-12;
            normalization.gamma          = 1.0;
            normalization.selectedCharge = 4.0e-15;
            transfer.depositCharge(
                    particles, workspace,
                    ParticleMeshFieldTransfer<double, 3>::DepositKind::PrimaryAndImage,
                    ParticleMeshFieldTransfer<double, 3>::Selection::direct(0, 2), normalization,
                    {.enabled = true, .planeZ = 0.0});

            const auto restoredPositions =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.R.getView());
            const auto restoredTimeSteps = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), particles.dt.getView());
            const auto restoredCharge =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.getQView());
            for (std::size_t index = 0; index < 2; ++index) {
                for (unsigned dimension = 0; dimension < 3; ++dimension) {
                    EXPECT_DOUBLE_EQ(
                            restoredPositions(index)[dimension], positions(index)[dimension]);
                }
                EXPECT_DOUBLE_EQ(restoredTimeSteps(index), timeSteps(index));
            }
            EXPECT_DOUBLE_EQ(restoredCharge(0), originalCharge(0));

            const auto chargeDensity = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), workspace.chargeDensity().getView());
            double absoluteChargeDensity = 0.0;
            for (std::size_t i = 0; i < chargeDensity.extent(0); ++i) {
                for (std::size_t j = 0; j < chargeDensity.extent(1); ++j) {
                    for (std::size_t k = 0; k < chargeDensity.extent(2); ++k) {
                        absoluteChargeDensity += std::abs(chargeDensity(i, j, k));
                    }
                }
            }
            EXPECT_GT(absoluteChargeDensity, 0.0);

            workspace.electricField() = Vector_t<double, 3>(1.0, 2.0, 3.0);
            transfer.gatherVector(
                    particles.E, workspace.electricField(), particles.R,
                    ParticleMeshFieldTransfer<double, 3>::GatherMode::Replace);
            const auto gathered =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.E.getView());
            for (std::size_t index = 0; index < 2; ++index) {
                for (unsigned dimension = 0; dimension < 3; ++dimension) {
                    EXPECT_NEAR(
                            gathered(index)[dimension], static_cast<double>(dimension + 1),
                            1.0e-14);
                }
            }
        }

        TEST_F(CartesianPICComponentsTest, RelativisticComposerHandlesDirectAndMirroredSources) {
            auto setup = storage();
            CartesianDomain<double, 3> domain(setup);
            CartesianPICFieldStorage<double, 3> workspace(domain);
            workspace.initializeFields("OPEN");

            auto source = workspace.electricField().getHostMirror();
            for (std::size_t i = 0; i < source.extent(0); ++i) {
                for (std::size_t j = 0; j < source.extent(1); ++j) {
                    for (std::size_t k = 0; k < source.extent(2); ++k) {
                        source(i, j, k) = Vector_t<double, 3>(
                                static_cast<double>(i + 1), static_cast<double>(j + 2),
                                static_cast<double>(k + 3));
                    }
                }
            }
            Kokkos::deep_copy(workspace.electricField().getView(), source);

            RelativisticFieldComposer<double, 3> composer;
            RelativisticFieldComposer<double, 3>::Policy policy;
            policy.gamma        = 1.0;
            policy.sourceRule   = FieldSourceRule::Direct;
            policy.magneticSign = 1.0;
            composer.clearAccumulation(workspace);
            composer.accumulate(workspace, policy);
            Kokkos::fence();

            const int ghost     = workspace.electricField().getNghost();
            auto directElectric = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), workspace.accumulatedElectricField().getView());
            auto directMagnetic = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), workspace.accumulatedMagneticField().getView());
            for (unsigned dimension = 0; dimension < 3; ++dimension) {
                EXPECT_DOUBLE_EQ(
                        directElectric(ghost, ghost, ghost)[dimension],
                        source(ghost, ghost, ghost)[dimension]);
                EXPECT_DOUBLE_EQ(directMagnetic(ghost, ghost, ghost)[dimension], 0.0);
            }

            composer.clearAccumulation(workspace);
            policy.sourceRule = FieldSourceRule::ShiftedGreenImageZ;
            composer.accumulate(workspace, policy);
            Kokkos::fence();
            const auto mirroredElectric = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), workspace.accumulatedElectricField().getView());
            const std::size_t mirroredK =
                    static_cast<std::size_t>(ghost) + workspace.layoutExtents()[2] - 1;
            EXPECT_DOUBLE_EQ(
                    mirroredElectric(ghost, ghost, ghost)[0], -source(ghost, ghost, mirroredK)[0]);
            EXPECT_DOUBLE_EQ(
                    mirroredElectric(ghost, ghost, ghost)[1], -source(ghost, ghost, mirroredK)[1]);
            EXPECT_DOUBLE_EQ(
                    mirroredElectric(ghost, ghost, ghost)[2], source(ghost, ghost, mirroredK)[2]);
        }

        TEST_F(CartesianPICComponentsTest, P3MShortRangeInteractionProducesFinitePairField) {
            auto setup          = storage();
            setup.layoutType    = ParticleLayoutType::SpatialOverlap;
            setup.overlapCutoff = 0.5;
            CartesianDomain<double, 3> domain(setup);
            using Container = ::ParticleContainer<double, 3>;
            Container particles(
                    domain.mesh(), domain.layout(), false, Container::LayoutType::SpatialOverlap,
                    setup.overlapCutoff, ippl::BC::NO);
            particles.setBunchStateHandler(std::make_shared<BunchStateHandler>());
            particles.createParticles(2);
            particles.setQ(1.0e-15);
            auto positions = particles.R.getHostMirror();
            positions(0)   = Vector_t<double, 3>(-0.1, 0.0, 0.0);
            positions(1)   = Vector_t<double, 3>(0.1, 0.0, 0.0);
            Kokkos::deep_copy(particles.R.getView(), positions);
            particles.E = Vector_t<double, 3>(0.0);
            particles.update();

            P3MShortRangeInteraction interaction(setup.overlapCutoff);
            interaction.apply(particles);
            const auto electric =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles.E.getView());
            double magnitude = 0.0;
            for (std::size_t index = 0; index < particles.getLocalNum(); ++index) {
                for (unsigned dimension = 0; dimension < 3; ++dimension) {
                    EXPECT_TRUE(std::isfinite(electric(index)[dimension]));
                    magnitude += std::abs(electric(index)[dimension]);
                }
            }
            EXPECT_GT(magnitude, 0.0);
        }

    }  // namespace
}  // namespace opalx::spacecharge
