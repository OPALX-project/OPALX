#include <gtest/gtest.h>

#include "Ippl.h"
#include "Utility/Inform.h"

#include "AbsBeamline/BeamBeamDefinitions.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/ElementInteractionManager.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/BeamBeamRep.h"
#include "BeamlineCore/DriftRep.h"
#include "PartBunch/PartBunch.h"
#include "SpaceCharge/CartesianPIC3D/CartesianDomainUpdater.h"
#include "SpaceCharge/Poisson/PoissonSolver.h"
#include "SpaceCharge/SpaceChargeConfigBuilder.h"
#include "Structure/Beam.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/Options.h"

#include <cstdint>
#include <memory>
#include <string>
#include <variant>
#include <vector>

extern Inform* gmsg;

namespace {

    using Bunch_t  = PartBunch<double, 3>;
    using Vector3d = ippl::Vector<double, 3>;
    namespace sc   = opalx::spacecharge;

    /**
     * @brief Test-only FieldSolverCmd helper exposing the minimal attribute setters
     * needed to construct a real PartBunch.
     *
     * The BeamBeam-window unit tests intentionally use the real PartBunch type,
     * because the bugs we fixed were caused by the coupling between:
     * - cached physical bunch bounds,
     * - field-domain bounds, and
     * - the solver-owned fixed Cartesian domain used by BeamBeam.
     *
     * Using a real bunch here gives coverage over those interfaces without pulling
     * in the full ParallelTracker machinery.
     */
    class TestableFieldSolverCmd : public FieldSolverCmd {
    public:
        void setType(const std::string& value) {
            Attributes::setPredefinedString(itsAttr[FIELDSOLVER::TYPE], value);
        }

        void setBoxIncr(double value) {
            Attributes::setReal(itsAttr[FIELDSOLVER::BBOXINCR], value);
        }

        void setParallelDecomposition(bool value) {
            Attributes::setBool(itsAttr[FIELDSOLVER::PARFFTX], value);
            Attributes::setBool(itsAttr[FIELDSOLVER::PARFFTY], value);
            Attributes::setBool(itsAttr[FIELDSOLVER::PARFFTZ], value);
        }
    };

    std::shared_ptr<TestableFieldSolverCmd> fieldSolverForBunch;
    std::shared_ptr<Beam> beamForBunch;

    std::shared_ptr<TestableFieldSolverCmd> makeFieldSolverCmd() {
        auto fieldSolver = std::make_shared<TestableFieldSolverCmd>();
        fieldSolver->setType("OPEN");
        fieldSolver->setNX(8);
        fieldSolver->setNY(8);
        fieldSolver->setNZ(8);
        fieldSolver->setBoxIncr(1.0);
        fieldSolver->setParallelDecomposition(true);
        fieldSolver->execute();
        return fieldSolver;
    }

    std::shared_ptr<Bunch_t> makeBunch() {
        Beam* beam = Beam::find("UNNAMED_BEAM");
        if (beam == nullptr) {
            beam = beamForBunch.get();
        }

        return std::make_shared<Bunch_t>(
                std::vector<double>{-1.0e-15}, std::vector<double>{5.10999e-4},
                std::vector<Beam*>{beam}, std::vector<size_t>{16}, 1.0, "LF2",
                sc::makeCartesianDomainConfig(
                        sc::buildSpaceChargeConfig(*fieldSolverForBunch, {})));
    }

    std::shared_ptr<Bunch_t> makeTwoContainerBunch() {
        Beam* beam = Beam::find("UNNAMED_BEAM");
        if (beam == nullptr) {
            beam = beamForBunch.get();
        }

        return std::make_shared<Bunch_t>(
                std::vector<double>{-1.0e-15, 1.0e-15}, std::vector<double>{5.10999e-4, 5.10999e-4},
                std::vector<Beam*>{beam, beam}, std::vector<size_t>{16, 16}, 1.0, "LF2",
                sc::makeCartesianDomainConfig(
                        sc::buildSpaceChargeConfig(*fieldSolverForBunch, {})));
    }

    void setParticlePositions(
            const std::shared_ptr<Bunch_t>& bunch, const std::vector<Vector3d>& positions,
            size_t containerIndex = 0) {
        auto pc = bunch->getParticleContainer(containerIndex);
        if (pc->getLocalNum() == 0) {
            pc->createParticles(positions.size());
        } else {
            ASSERT_EQ(pc->getLocalNum(), positions.size());
        }

        auto R_host = pc->R.getHostMirror();
        auto P_host = pc->P.getHostMirror();
        for (size_t i = 0; i < positions.size(); ++i) {
            R_host(i) = positions[i];
            P_host(i) = Vector3d(0.0);
        }

        Kokkos::deep_copy(pc->R.getView(), R_host);
        Kokkos::deep_copy(pc->P.getView(), P_host);
        Kokkos::fence();
        pc->markMomentsDirty();
    }

    // Exercise the same domain updater used by the Cartesian solver without a field solve:
    // these tests cover geometry and migration, independently of the Poisson discretization.
    void updateFieldDomain(
            const std::shared_ptr<Bunch_t>& bunch,
            sc::DomainCoordinateFrame frame = sc::DomainCoordinateFrame::Beam) {
        auto config = std::get<sc::CartesianPIC3DConfig>(
                sc::buildSpaceChargeConfig(*fieldSolverForBunch, {}));
        config.backend = sc::PoissonSolverType::None;
        sc::CartesianPIC3DFieldStorage<double, 3> fields(bunch->cartesianDomain());
        fields.initializeFields(config.backend);
        auto poisson = sc::makePoissonSolver(
                sc::makePoissonSolverConfig(config),
                {&fields.chargeDensity(), &fields.electricField()});
        std::vector<Bunch_t::ParticleContainer_t*> particles;
        for (const auto& pc : bunch->getParticleContainers()) {
            particles.push_back(pc.get());
        }
        sc::CartesianDomainUpdater updater(config, particles);
        std::vector<std::uint8_t> activity(particles.size(), 0);
        activity.front() = 1;
        sc::SpaceChargeStepState step;
        step.mpiSize = ippl::Comm->size();
        const sc::SpaceChargeSolveContext context(activity, step);
        const auto& fixed = bunch->getBunchStateHandler()->fixedCartesianDomain();
        const auto* fixedDomain =
                frame == sc::DomainCoordinateFrame::Beam && fixed ? &*fixed : nullptr;
        EXPECT_FALSE(updater.updateForSolve(frame, context, {}, fixedDomain, fields, *poisson));
    }

    void setFixedDomain(
            const std::shared_ptr<Bunch_t>& bunch, const Vector3d& lower, const Vector3d& upper) {
        bunch->getBunchStateHandler()->setFixedCartesianDomain(
                {lower[0], lower[1], lower[2]}, {upper[0], upper[1], upper[2]});
    }

    void expectVectorNear(const Vector3d& actual, const Vector3d& expected, double tolerance) {
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_NEAR(actual[d], expected[d], tolerance) << "Mismatch in component " << d;
        }
    }

    class BeamBeamPartBunchTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
            OpalData::getInstance()->storeInputFn("unit_test.opal");
            // PartBunch teardown writes through the global OPAL Inform stream.
            // The regular `tests/` binary creates this in its custom main(), but
            // the unit_tests tree uses gtest_main, so this fixture must provide
            // the minimal global setup needed by a standalone PartBunch instance.
            gmsg                = new Inform("UnitTests: ", std::cerr);
            Options::enableHDF5 = false;

            fieldSolverForBunch = makeFieldSolverCmd();
            beamForBunch        = std::make_shared<Beam>();
        }

        static void TearDownTestSuite() {
            beamForBunch.reset();
            fieldSolverForBunch.reset();
            delete gmsg;
            gmsg = nullptr;
            ippl::finalize();
        }
    };

    // ----------------------------------------------------------------------------
    // Fixed BeamBeam domain lifecycle
    // ----------------------------------------------------------------------------

    // The shared state owns geometry intent; the Cartesian solver owns its application.
    TEST_F(BeamBeamPartBunchTest, FixedBeamBeamDomainLifecycle) {
        auto bunch = makeBunch();
        auto state = bunch->getBunchStateHandler();

        ASSERT_FALSE(state->fixedCartesianDomain());

        setFixedDomain(bunch, Vector3d(-0.2, -0.3, 1.0), Vector3d(0.2, 0.3, 1.5));
        ASSERT_TRUE(state->fixedCartesianDomain());

        EXPECT_DOUBLE_EQ(state->fixedCartesianDomain()->lower[2], 1.0);
        EXPECT_DOUBLE_EQ(state->fixedCartesianDomain()->upper[2], 1.5);

        state->clearFixedCartesianDomain();
        EXPECT_FALSE(state->fixedCartesianDomain());
    }

    TEST_F(BeamBeamPartBunchTest, OnlyStatefulElementsCreateRuntimeInteractions) {
        auto bunch    = makeBunch();
        auto beamBeam = std::make_shared<BeamBeamRep>("BB");
        auto drift    = std::make_shared<DriftRep>("D");

        ElementInteractionManager manager;
        manager.initialize({beamBeam, drift});

        ASSERT_EQ(manager.size(), 1u);
        EXPECT_FALSE(manager.freezesFieldMesh());
        EXPECT_FALSE(manager.suppressesDefaultSelfField());

        ElementInteractionContext context{*bunch};
        const auto result = manager.execute(ElementInteractionPhase::Diagnostics, context);
        EXPECT_FALSE(result.selfFieldHandled);
    }

    TEST_F(BeamBeamPartBunchTest, WitnessContainerMaskDecodesConfiguredContainers) {
        EXPECT_TRUE(BEAMBEAM::decodeWitnessContainerMask(0.0).empty());

        const auto witnesses = BEAMBEAM::decodeWitnessContainerMask((1ULL << 1U) | (1ULL << 2U));
        ASSERT_EQ(witnesses.size(), 2u);
        EXPECT_EQ(witnesses[0], 1u);
        EXPECT_EQ(witnesses[1], 2u);
    }

    TEST_F(BeamBeamPartBunchTest, RigidSourceOnlyDisablesSourceCollectiveKick) {
        BEAMBEAM::Config config;
        EXPECT_TRUE(BEAMBEAM::sourceCollectiveKickEnabled(config));

        config.rigidSource = true;
        EXPECT_FALSE(BEAMBEAM::sourceCollectiveKickEnabled(config));
    }

    TEST_F(BeamBeamPartBunchTest, InteractionPointIsPlacedElementMidpoint) {
        EXPECT_DOUBLE_EQ(BEAMBEAM::interactionPointAtElementMidpoint(1.0, 1.5), 1.25);
        EXPECT_DOUBLE_EQ(BEAMBEAM::interactionPointAtElementMidpoint(-0.004, 0.012), 0.004);
    }

    TEST_F(BeamBeamPartBunchTest, CollectiveFieldWindowIsTwentyMillimetresAroundIp) {
        constexpr double ip = 0.17;
        EXPECT_DOUBLE_EQ(BEAMBEAM::fieldWindowLength, 20.0e-3);
        EXPECT_DOUBLE_EQ(BEAMBEAM::fieldWindowBegin(ip), 0.16);
        EXPECT_DOUBLE_EQ(BEAMBEAM::fieldWindowEnd(ip), 0.18);
    }

    TEST_F(BeamBeamPartBunchTest, SourceExitsOnlyAfterItsTailPassesFieldWindow) {
        constexpr double ip = 0.17;
        EXPECT_FALSE(BEAMBEAM::sourceFullyExitedFieldWindow(0.179, ip));
        EXPECT_FALSE(BEAMBEAM::sourceFullyExitedFieldWindow(0.180, ip));
        EXPECT_TRUE(BEAMBEAM::sourceFullyExitedFieldWindow(0.181, ip));
    }

    TEST_F(BeamBeamPartBunchTest, WitnessLongitudinalOffsetMapsIpToSourceFrame) {
        const double sourceS  = 30.0e-3;
        const double witnessS = 0.0;
        const double witnessZ = 35.0e-3;

        const double sourceFrameZ =
                witnessZ + BEAMBEAM::longitudinalOffsetToSourceFrame(sourceS, witnessS);

        EXPECT_NEAR(sourceFrameZ, 5.0e-3, 1.0e-15);
    }

    TEST_F(BeamBeamPartBunchTest, CopyTimeUsesConfiguredThreshold) {
        EXPECT_FALSE(BEAMBEAM::copyTimeReached(200.0e-12, std::nullopt));
        const std::optional<double> copyTime = 100.0e-12;
        EXPECT_FALSE(BEAMBEAM::copyTimeReached(99.0e-12, copyTime));
        EXPECT_TRUE(BEAMBEAM::copyTimeReached(100.0e-12, copyTime));
        EXPECT_TRUE(BEAMBEAM::copyTimeReached(101.0e-12, copyTime));
    }

    TEST_F(BeamBeamPartBunchTest, CopiedSourceOverlapUsesMirroredIntervalAtIp) {
        BEAMBEAM::ActualGeometry geometry;
        geometry.interactionPointS = 1.25;

        EXPECT_FALSE(BEAMBEAM::copiedSourceBunchesOverlap(1.00, 1.20, geometry));
        EXPECT_TRUE(BEAMBEAM::copiedSourceBunchesOverlap(1.00, 1.25, geometry));
        EXPECT_TRUE(BEAMBEAM::copiedSourceBunchesOverlap(1.20, 1.30, geometry));
        EXPECT_TRUE(BEAMBEAM::copiedSourceBunchesOverlap(1.25, 1.50, geometry));
        EXPECT_FALSE(BEAMBEAM::copiedSourceBunchesOverlap(1.30, 1.50, geometry));
    }

    // ----------------------------------------------------------------------------
    // Solver-owned fixed window geometry
    // ----------------------------------------------------------------------------

    // The BeamBeam-window mesh switch should keep the transverse field domain
    // unchanged while replacing only the longitudinal mesh extent. This is the
    // explicit Lagrangian-to-Eulerian transition in z that the current model relies on.
    TEST_F(BeamBeamPartBunchTest, FixedBeamBeamDomainOnlyChangesLongitudinalDomain) {
        auto bunch   = makeBunch();
        auto& domain = bunch->cartesianDomain();

        const Vector3d initialRMin = domain.lower();
        const Vector3d initialRMax = domain.upper();
        setFixedDomain(
                bunch, Vector3d(initialRMin[0], initialRMin[1], 1.0),
                Vector3d(initialRMax[0], initialRMax[1], 1.5));
        updateFieldDomain(bunch);

        const Vector3d updatedRMin = domain.lower();
        const Vector3d updatedRMax = domain.upper();
        const Vector3d updatedHr   = domain.spacing();
        const Vector3d meshOrigin  = domain.mesh().getOrigin();

        EXPECT_DOUBLE_EQ(updatedRMin[0], initialRMin[0]);
        EXPECT_DOUBLE_EQ(updatedRMin[1], initialRMin[1]);
        EXPECT_DOUBLE_EQ(updatedRMax[0], initialRMax[0]);
        EXPECT_DOUBLE_EQ(updatedRMax[1], initialRMax[1]);
        // Current master represents RMin/RMax as the first/last cell centers.
        EXPECT_DOUBLE_EQ(updatedHr[0], (initialRMax[0] - initialRMin[0]) / 7.0);
        EXPECT_DOUBLE_EQ(updatedHr[1], (initialRMax[1] - initialRMin[1]) / 7.0);

        EXPECT_DOUBLE_EQ(updatedRMin[2], 1.00);
        EXPECT_DOUBLE_EQ(updatedRMax[2], 1.50);
        EXPECT_DOUBLE_EQ(updatedHr[2], 0.50 / 7.0);

        expectVectorNear(meshOrigin, updatedRMin - 0.5 * updatedHr, 1.0e-14);
    }

    TEST_F(BeamBeamPartBunchTest, FixedBeamBeamDomainUsesExplicitTransverseBounds) {
        auto bunch   = makeBunch();
        auto& domain = bunch->cartesianDomain();

        setFixedDomain(bunch, Vector3d(-2.0e-3, -3.0e-3, 1.0), Vector3d(2.0e-3, 3.0e-3, 1.5));
        updateFieldDomain(bunch);

        const Vector3d updatedRMin = domain.lower();
        const Vector3d updatedRMax = domain.upper();
        const Vector3d updatedHr   = domain.spacing();
        const Vector3d meshOrigin  = domain.mesh().getOrigin();

        EXPECT_DOUBLE_EQ(updatedRMin[0], -2.0e-3);
        EXPECT_DOUBLE_EQ(updatedRMax[0], 2.0e-3);
        EXPECT_DOUBLE_EQ(updatedRMin[1], -3.0e-3);
        EXPECT_DOUBLE_EQ(updatedRMax[1], 3.0e-3);
        EXPECT_DOUBLE_EQ(updatedRMin[2], 1.00);
        EXPECT_DOUBLE_EQ(updatedRMax[2], 1.50);

        EXPECT_DOUBLE_EQ(updatedHr[0], 4.0e-3 / 7.0);
        EXPECT_DOUBLE_EQ(updatedHr[1], 6.0e-3 / 7.0);
        EXPECT_DOUBLE_EQ(updatedHr[2], 0.50 / 7.0);

        expectVectorNear(meshOrigin, updatedRMin - 0.5 * updatedHr, 1.0e-14);
    }

    TEST_F(BeamBeamPartBunchTest, OpenWindowMeshDoesNotWrapLongitudinalParticlePosition) {
        auto bunch = makeBunch();
        auto pc    = bunch->getParticleContainer();

        setParticlePositions(bunch, {Vector3d(0.0, 0.0, 0.30)});
        setFixedDomain(bunch, Vector3d(-1.0, -1.0, -0.25), Vector3d(1.0, 1.0, 0.25));
        updateFieldDomain(bunch);

        ASSERT_EQ(pc->getTotalNum(), static_cast<size_t>(ippl::Comm->size()));
        auto positions = pc->R.getHostMirror();
        Kokkos::deep_copy(positions, pc->R.getView());
        for (size_t i = 0; i < pc->getLocalNum(); ++i) {
            EXPECT_DOUBLE_EQ(positions(i)[2], 0.30);
        }
    }

    TEST_F(BeamBeamPartBunchTest, BeamBeamMinimumTransverseSpanLimitsRestFrameAspect) {
        constexpr double gamma = 480.453039972;
        const double span      = BEAMBEAM::minimumTransverseSpan(20.0e-3, 128, 256, gamma);
        const double dx        = span / 255.0;
        const double dzRest    = gamma * 20.0e-3 / 127.0;

        EXPECT_NEAR(span, 128.625e-6, 5.0e-10);
        EXPECT_NEAR(dzRest / dx, BEAMBEAM::maximumRestFrameCellAspectRatio, 1.0e-9);
        EXPECT_DOUBLE_EQ(BEAMBEAM::minimumTransverseSpan(20.0e-3, 1, 256, gamma), 0.0);
    }

    // ----------------------------------------------------------------------------
    // Return from the fixed field domain to particle-following geometry
    // ----------------------------------------------------------------------------

    // Master no longer snapshots physical bounds alongside fields. Clearing fixed-domain
    // intent lets the solver recompute geometry, while particle statistics remain physical.
    TEST_F(BeamBeamPartBunchTest, ClearFixedDomainRestoresMeshAndPreservesPhysicalBounds) {
        auto bunch   = makeBunch();
        auto& domain = bunch->cartesianDomain();

        const Vector3d physicalRMin(-0.2, -0.3, -0.4);
        const Vector3d physicalRMax(0.5, 0.6, 0.7);
        setParticlePositions(bunch, {physicalRMin, physicalRMax});
        updateFieldDomain(bunch);
        bunch->calcBeamParameters();

        const Vector3d savedLower   = domain.lower();
        const Vector3d savedUpper   = domain.upper();
        const Vector3d savedSpacing = domain.spacing();
        const Vector3d savedOrigin  = domain.origin();
        setFixedDomain(bunch, Vector3d(-1.0), Vector3d(1.0));
        updateFieldDomain(bunch);
        bunch->calcBeamParameters();

        Vector3d fixedPhysicalRMin(0.0), fixedPhysicalRMax(0.0);
        bunch->get_bounds(fixedPhysicalRMin, fixedPhysicalRMax);
        expectVectorNear(fixedPhysicalRMin, physicalRMin, 1.0e-14);
        expectVectorNear(fixedPhysicalRMax, physicalRMax, 1.0e-14);

        bunch->getBunchStateHandler()->clearFixedCartesianDomain();
        updateFieldDomain(bunch);
        bunch->calcBeamParameters();

        Vector3d restoredPhysicalRMin(0.0);
        Vector3d restoredPhysicalRMax(0.0);
        bunch->get_bounds(restoredPhysicalRMin, restoredPhysicalRMax);

        expectVectorNear(domain.lower(), savedLower, 1.0e-14);
        expectVectorNear(domain.upper(), savedUpper, 1.0e-14);
        expectVectorNear(domain.spacing(), savedSpacing, 1.0e-14);
        expectVectorNear(domain.origin(), savedOrigin, 1.0e-14);
        expectVectorNear(restoredPhysicalRMin, physicalRMin, 1.0e-14);
        expectVectorNear(restoredPhysicalRMax, physicalRMax, 1.0e-14);
    }

    // ----------------------------------------------------------------------------
    // Transverse refresh with a fixed longitudinal window
    // ----------------------------------------------------------------------------

    // Moments-only updates cannot mutate the mesh. The interaction explicitly refreshes
    // transverse domain intent, and the solver applies it without changing the fixed z extent.
    TEST_F(BeamBeamPartBunchTest, BeamBeamDomainKeepsZFixedWhileTransverseBoundsExpand) {
        auto bunch   = makeBunch();
        auto& domain = bunch->cartesianDomain();

        setParticlePositions(bunch, {Vector3d(-0.20, -0.10, -0.05), Vector3d(0.10, 0.20, 0.05)});
        setFixedDomain(bunch, Vector3d(-0.25, -0.15, -0.25), Vector3d(0.15, 0.25, 0.25));
        updateFieldDomain(bunch);

        const Vector3d frozenRMin = domain.lower();
        const Vector3d frozenRMax = domain.upper();
        const Vector3d frozenHr   = domain.spacing();

        auto pc        = bunch->getParticleContainer();
        auto positions = pc->R.getHostMirror();
        Kokkos::deep_copy(positions, pc->R.getView());
        for (size_t i = 0; i < pc->getLocalNum(); ++i) {
            positions(i)[0] *= 3.0;
            positions(i)[1] *= 2.0;
        }
        Kokkos::deep_copy(pc->R.getView(), positions);
        pc->markMomentsDirty();
        bunch->updateAllParticleMoments();
        expectVectorNear(domain.lower(), frozenRMin, 1.0e-14);
        expectVectorNear(domain.upper(), frozenRMax, 1.0e-14);

        bunch->getBunchStateHandler()->clearFixedCartesianDomain();
        setFixedDomain(bunch, Vector3d(-0.65, -0.25, -0.25), Vector3d(0.35, 0.45, 0.25));
        updateFieldDomain(bunch);

        const Vector3d updatedRMin = domain.lower();
        const Vector3d updatedRMax = domain.upper();
        const Vector3d updatedHr   = domain.spacing();

        EXPECT_DOUBLE_EQ(updatedRMin[2], frozenRMin[2]);
        EXPECT_DOUBLE_EQ(updatedRMax[2], frozenRMax[2]);
        EXPECT_DOUBLE_EQ(updatedHr[2], frozenHr[2]);

        EXPECT_LT(updatedRMin[0], frozenRMin[0]);
        EXPECT_GT(updatedRMax[0], frozenRMax[0]);
    }

    TEST_F(BeamBeamPartBunchTest, FieldBoundsIncludePassiveWitnessContainer) {
        auto bunch = makeTwoContainerBunch();
        setParticlePositions(
                bunch, {Vector3d(-1.0e-3, -2.0e-3, -3.0e-3), Vector3d(1.0e-3, 2.0e-3, 3.0e-3)});
        setParticlePositions(
                bunch, {Vector3d(-12.0e-3, -8.0e-3, -1.0e-3), Vector3d(15.0e-3, 9.0e-3, 1.0e-3)},
                1);

        updateFieldDomain(bunch, sc::DomainCoordinateFrame::Reference);
        const Vector3d lower = bunch->cartesianDomain().lower();
        const Vector3d upper = bunch->cartesianDomain().upper();

        EXPECT_LT(lower[0], -12.0e-3);
        EXPECT_GT(upper[0], 15.0e-3);
        EXPECT_LT(lower[1], -8.0e-3);
        EXPECT_GT(upper[1], 9.0e-3);
    }

}  // namespace
