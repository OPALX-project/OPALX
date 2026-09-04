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
#include "SpaceCharge/Poisson/PoissonSolver.h"
#include "SpaceCharge/SpaceChargeSolveContext.h"
#include "Structure/DataSink.h"
#include "Utilities/OpalException.h"
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

        TEST_F(CartesianPICComponentsTest, FixedDomainOverridesStretchingAndReturnsToFollowing) {
            auto domainConfig          = storage();
            domainConfig.decomposition = {true, true, true};
            CartesianDomain<double, 3> domain(domainConfig);
            CartesianPICFieldStorage<double, 3> workspace(domain);
            workspace.initializeFields(PoissonSolverType::None);

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

            std::vector<Container*> particleContainers{particles.get(), secondary.get()};

            PoissonSolverConfig poissonConfig;
            poissonConfig.type = PoissonSolverType::None;
            PoissonSolver poisson(
                    poissonConfig, {&workspace.chargeDensity(), &workspace.electricField()});

            CartesianPICConfig values;
            values.backend                         = PoissonSolverType::None;
            values.grid.meshSize                   = domainConfig.meshSize;
            values.grid.decomposition              = domainConfig.decomposition;
            values.layoutRebuildDecomposition      = std::array<bool, 3>{true, true, true};
            values.grid.boundingBoxIncreasePercent = 50.0;
            values.repartitionFrequency            = 1;
            CartesianDomainUpdater updater(values, particleContainers);

            const std::array<double, 3> fixedLower{-2.0, -3.0, -4.0};
            const std::array<double, 3> fixedUpper{2.0, 3.0, 4.0};
            bunchState->setFixedCartesianDomain(fixedLower, fixedUpper);
            SpaceChargeStepState step;
            step.emissionActive  = true;
            step.emittedFraction = 0.01;
            const std::array<std::uint8_t, 2> activity{1, 0};
            SpaceChargeSolveContext context(activity, step);

            EXPECT_FALSE(updater.updateForSolve(
                    DomainCoordinateFrame::Beam, context, {}, &*bunchState->fixedCartesianDomain(),
                    workspace, poisson));
            EXPECT_EQ(workspace.layoutExtents(), domainConfig.meshSize);
            EXPECT_EQ(domain.decomposition(), domainConfig.decomposition);
            EXPECT_FALSE(particles->isMomentsDirty());
            EXPECT_FALSE(secondary->isMomentsDirty());
            for (unsigned dimension = 0; dimension < 3; ++dimension) {
                EXPECT_DOUBLE_EQ(domain.lower()[dimension], fixedLower[dimension]);
                EXPECT_DOUBLE_EQ(domain.upper()[dimension], fixedUpper[dimension]);
            }

            secondary->markMomentsDirty();
            bunchState->clearFixedCartesianDomain();
            EXPECT_FALSE(updater.updateForSolve(
                    DomainCoordinateFrame::Beam, context, {}, nullptr, workspace, poisson));
            EXPECT_NE(domain.lower()[0], fixedLower[0]);
            EXPECT_NE(domain.upper()[0], fixedUpper[0]);
            EXPECT_FALSE(particles->isMomentsDirty());
            EXPECT_TRUE(secondary->isMomentsDirty());

            EXPECT_FALSE(updater.updateForSolve(
                    DomainCoordinateFrame::Reference, context, {}, nullptr, workspace, poisson));
            EXPECT_FALSE(secondary->isMomentsDirty());
        }

        TEST_F(CartesianPICComponentsTest, CartesianAlgorithmRetainsAndClearsFixedMesh) {
            auto domainConfig          = storage();
            domainConfig.decomposition = {true, true, true};
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
            particles->markMomentsDirty();
            particles->updateMoments();
            const Vector_t<double, 3> originalMean = particles->getMeanR();

            std::vector<Container*> particleContainers{particles.get()};
            CartesianPICConfig values;
            values.backend                         = PoissonSolverType::Open;
            values.grid.meshSize                   = domainConfig.meshSize;
            values.grid.decomposition              = domainConfig.decomposition;
            values.grid.boundingBoxIncreasePercent = 10.0;
            DataSink dataSink;
            CartesianPICAlgorithm algorithm(
                    values, particleContainers, std::move(fieldStorage), &dataSink, bunchState);

            const std::array<double, 3> fixedLower{-1.0, -1.5, -2.0};
            const std::array<double, 3> fixedUpper{1.0, 1.5, 2.0};
            bunchState->setFixedCartesianDomain(fixedLower, fixedUpper);
            SpaceChargeStepState step;
            step.timeStep = 1.0e-12;
            const CoordinateSystemTrafo trackerToSolve(
                    Vector_t<double, 3>(0.1, -0.2, 0.3),
                    Quaternion(std::cos(Physics::pi / 4.0), 0.0, 0.0, std::sin(Physics::pi / 4.0)));
            step.frames = {trackerToSolve, trackerToSolve.inverted()};
            const std::array<std::uint8_t, 1> activity{1};
            SpaceChargeSolveContext fixedContext(activity, step);
            const SpaceChargeSolveResult fixedResult = algorithm.solve(fixedContext);
            EXPECT_EQ(fixedResult.backendSolves, 1u);
            for (unsigned dimension = 0; dimension < 3; ++dimension) {
                EXPECT_DOUBLE_EQ(fieldObserver->lower()[dimension], fixedLower[dimension]);
                EXPECT_DOUBLE_EQ(fieldObserver->upper()[dimension], fixedUpper[dimension]);
            }
            EXPECT_FALSE(particles->isMomentsDirty());
            for (unsigned dimension = 0; dimension < 3; ++dimension) {
                EXPECT_NEAR(particles->getMeanR()[dimension], originalMean[dimension], 1.0e-14);
            }

            bunchState->clearFixedCartesianDomain();
            SpaceChargeSolveContext followingContext(activity, step);
            const SpaceChargeSolveResult followingResult = algorithm.solve(followingContext);
            EXPECT_EQ(followingResult.backendSolves, 1u);
            EXPECT_NE(fieldObserver->lower()[0], fixedLower[0]);
            EXPECT_NE(fieldObserver->upper()[0], fixedUpper[0]);
            EXPECT_FALSE(particles->isMomentsDirty());
        }

    }  // namespace
}  // namespace opalx::spacecharge
