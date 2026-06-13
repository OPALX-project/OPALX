#include <gtest/gtest.h>

#include "Ippl.h"
#include "Utility/Inform.h"

#include "AbsBeamline/BeamBeamDefinitions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "PartBunch/PartBunch.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/Options.h"

#include <memory>
#include <string>
#include <vector>

extern Inform* gmsg;

namespace {

    using Bunch_t  = PartBunch<double, 3>;
    using Vector3d = ippl::Vector<double, 3>;

    /**
     * @brief Test-only FieldSolverCmd helper exposing the minimal attribute setters
     * needed to construct a real PartBunch.
     *
     * The BeamBeam-window unit tests intentionally use the real PartBunch type,
     * because the bugs we fixed were caused by the coupling between:
     * - cached physical bunch bounds,
     * - field-domain bounds, and
     * - the temporary BeamBeam-window mesh replacement.
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
    };

    std::shared_ptr<TestableFieldSolverCmd> fieldSolverForBunch;
    std::shared_ptr<DataSink> dataSinkForBunch;
    std::shared_ptr<Beam> beamForBunch;

    std::shared_ptr<TestableFieldSolverCmd> makeFieldSolverCmd() {
        auto fieldSolver = std::make_shared<TestableFieldSolverCmd>();
        fieldSolver->setType("OPEN");
        fieldSolver->setNX(8);
        fieldSolver->setNY(8);
        fieldSolver->setNZ(8);
        fieldSolver->setBoxIncr(1.0);
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
                fieldSolverForBunch.get(), dataSinkForBunch.get());
    }

    void setParticlePositions(
            const std::shared_ptr<Bunch_t>& bunch, const std::vector<Vector3d>& positions) {
        auto pc = bunch->getParticleContainer();
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
        pc->update();
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
            dataSinkForBunch    = std::make_shared<DataSink>();
            beamForBunch        = std::make_shared<Beam>();
        }

        static void TearDownTestSuite() {
            beamForBunch.reset();
            dataSinkForBunch.reset();
            fieldSolverForBunch.reset();
            delete gmsg;
            gmsg = nullptr;
            ippl::finalize();
        }
    };

    // ----------------------------------------------------------------------------
    // BeamBeamWindowConfig lifecycle
    // ----------------------------------------------------------------------------

    // Verifies that the bunch-side BeamBeam-window configuration behaves like a
    // narrow geometry input only: absent by default, present after configuration,
    // and removable again without affecting unrelated state.
    TEST_F(BeamBeamPartBunchTest, BeamBeamWindowConfigLifecycle) {
        auto bunch = makeBunch();

        ASSERT_FALSE(bunch->hasBeamBeamWindowConfig());

        bunch->setBeamBeamWindowConfig(0.50, 1.25, 1.00, 1.50, false);
        ASSERT_TRUE(bunch->hasBeamBeamWindowConfig());

        const auto& config = bunch->getBeamBeamWindowConfig();
        EXPECT_DOUBLE_EQ(config.interactionPointS, 1.25);
        EXPECT_DOUBLE_EQ(config.windowBeginS, 1.00);
        EXPECT_DOUBLE_EQ(config.windowEndS, 1.50);
        EXPECT_FALSE(config.copyModel);

        bunch->clearBeamBeamWindowConfig();
        EXPECT_FALSE(bunch->hasBeamBeamWindowConfig());
    }

    TEST_F(BeamBeamPartBunchTest, WitnessContainerMaskDecodesConfiguredContainers) {
        EXPECT_TRUE(BEAMBEAM::decodeWitnessContainerMask(0.0).empty());

        const auto witnesses = BEAMBEAM::decodeWitnessContainerMask((1ULL << 1U) | (1ULL << 2U));
        ASSERT_EQ(witnesses.size(), 2u);
        EXPECT_EQ(witnesses[0], 1u);
        EXPECT_EQ(witnesses[1], 2u);
    }

    TEST_F(BeamBeamPartBunchTest, WitnessLongitudinalOffsetMapsIpToSourceFrame) {
        const double sourceS  = 30.0e-3;
        const double witnessS = 0.0;
        const double witnessZ = 35.0e-3;

        const double sourceFrameZ =
                witnessZ + BEAMBEAM::longitudinalOffsetToSourceFrame(sourceS, witnessS);

        EXPECT_NEAR(sourceFrameZ, 5.0e-3, 1.0e-15);
    }

    TEST_F(BeamBeamPartBunchTest, SourceRetirementTimeUsesConfiguredThreshold) {
        EXPECT_FALSE(BEAMBEAM::sourceRetireTimeReached(200.0e-12, std::nullopt));
        const std::optional<double> retireTime = 100.0e-12;
        EXPECT_FALSE(BEAMBEAM::sourceRetireTimeReached(99.0e-12, retireTime));
        EXPECT_TRUE(BEAMBEAM::sourceRetireTimeReached(100.0e-12, retireTime));
        EXPECT_TRUE(BEAMBEAM::sourceRetireTimeReached(101.0e-12, retireTime));
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

    TEST_F(BeamBeamPartBunchTest, MarkAllParticlesInvalidRetiresAllParticles) {
        auto bunch = makeBunch();
        setParticlePositions(
                bunch, {Vector3d(-1.0e-3, 0.0, 0.0), Vector3d(0.0, 1.0e-3, 2.0e-3),
                        Vector3d(1.0e-3, 0.0, 4.0e-3)});

        auto pc                  = bunch->getParticleContainer();
        const size_t totalBefore = pc->getTotalNum();
        ASSERT_GT(totalBefore, 0u);

        EXPECT_EQ(pc->markAllParticlesInvalid(), totalBefore);
        EXPECT_EQ(pc->deleteInvalidParticles(), totalBefore);
        EXPECT_EQ(pc->getTotalNum(), 0u);
    }

    // ----------------------------------------------------------------------------
    // enableBeamBeamWindowMesh
    // ----------------------------------------------------------------------------

    // The BeamBeam-window mesh switch should keep the transverse field domain
    // unchanged while replacing only the longitudinal mesh extent. This is the
    // explicit Lagrangian-to-Eulerian transition in z that the current model relies on.
    TEST_F(BeamBeamPartBunchTest, EnableBeamBeamWindowMeshOnlyChangesLongitudinalDomain) {
        auto bunch          = makeBunch();
        auto fieldContainer = bunch->getFieldContainer();

        const Vector3d initialRMin = fieldContainer->getRMin();
        const Vector3d initialRMax = fieldContainer->getRMax();
        const Vector3d initialHr   = fieldContainer->getHr();

        bunch->enableBeamBeamWindowMesh(1.25, 0.50);

        const Vector3d updatedRMin = fieldContainer->getRMin();
        const Vector3d updatedRMax = fieldContainer->getRMax();
        const Vector3d updatedHr   = fieldContainer->getHr();
        const Vector3d meshOrigin  = fieldContainer->getMesh().getOrigin();

        EXPECT_DOUBLE_EQ(updatedRMin[0], initialRMin[0]);
        EXPECT_DOUBLE_EQ(updatedRMin[1], initialRMin[1]);
        EXPECT_DOUBLE_EQ(updatedRMax[0], initialRMax[0]);
        EXPECT_DOUBLE_EQ(updatedRMax[1], initialRMax[1]);
        EXPECT_DOUBLE_EQ(updatedHr[0], initialHr[0]);
        EXPECT_DOUBLE_EQ(updatedHr[1], initialHr[1]);

        EXPECT_DOUBLE_EQ(updatedRMin[2], 1.00);
        EXPECT_DOUBLE_EQ(updatedRMax[2], 1.50);
        EXPECT_DOUBLE_EQ(updatedHr[2], 0.50 / 8.0);

        expectVectorNear(meshOrigin, updatedRMin, 1.0e-14);
    }

    TEST_F(BeamBeamPartBunchTest, EnableBeamBeamWindowMeshUsesExplicitTransverseAperture) {
        auto bunch          = makeBunch();
        auto fieldContainer = bunch->getFieldContainer();

        bunch->enableBeamBeamWindowMesh(1.25, 0.50, 2.0e-3, 3.0e-3);

        const Vector3d updatedRMin = fieldContainer->getRMin();
        const Vector3d updatedRMax = fieldContainer->getRMax();
        const Vector3d updatedHr   = fieldContainer->getHr();
        const Vector3d meshOrigin  = fieldContainer->getMesh().getOrigin();

        EXPECT_DOUBLE_EQ(updatedRMin[0], -2.0e-3);
        EXPECT_DOUBLE_EQ(updatedRMax[0], 2.0e-3);
        EXPECT_DOUBLE_EQ(updatedRMin[1], -3.0e-3);
        EXPECT_DOUBLE_EQ(updatedRMax[1], 3.0e-3);
        EXPECT_DOUBLE_EQ(updatedRMin[2], 1.00);
        EXPECT_DOUBLE_EQ(updatedRMax[2], 1.50);

        EXPECT_DOUBLE_EQ(updatedHr[0], 4.0e-3 / 8.0);
        EXPECT_DOUBLE_EQ(updatedHr[1], 6.0e-3 / 8.0);
        EXPECT_DOUBLE_EQ(updatedHr[2], 0.50 / 8.0);

        expectVectorNear(meshOrigin, updatedRMin, 1.0e-14);
    }

    // ----------------------------------------------------------------------------
    // save/restore field-domain state
    // ----------------------------------------------------------------------------

    // The field-domain state roundtrip must restore both:
    // - the mesh-aligned field-domain quantities, and
    // - the cached physical bunch bounds returned by get_bounds().
    //
    // This test covers the exact contract that the BeamBeam-window mode relies on
    // when it temporarily swaps in an Eulerian longitudinal mesh and later returns to
    // the original bunch-following state.
    TEST_F(BeamBeamPartBunchTest, RestoreFieldDomainStateRestoresMeshAndPhysicalBounds) {
        auto bunch          = makeBunch();
        auto fieldContainer = bunch->getFieldContainer();

        const Vector3d physicalRMin(-0.2, -0.3, -0.4);
        const Vector3d physicalRMax(0.5, 0.6, 0.7);
        bunch->setPhysicalBounds(physicalRMin, physicalRMax);

        const auto savedState = bunch->saveFieldDomainState();
        bunch->enableBeamBeamWindowMesh(1.25, 0.50);
        bunch->setPhysicalBounds(Vector3d(10.0), Vector3d(11.0));

        bunch->restoreFieldDomainState(savedState);

        Vector3d restoredPhysicalRMin(0.0);
        Vector3d restoredPhysicalRMax(0.0);
        bunch->get_bounds(restoredPhysicalRMin, restoredPhysicalRMax);

        expectVectorNear(fieldContainer->getRMin(), savedState.rmin, 1.0e-14);
        expectVectorNear(fieldContainer->getRMax(), savedState.rmax, 1.0e-14);
        expectVectorNear(fieldContainer->getHr(), savedState.hr, 1.0e-14);
        expectVectorNear(fieldContainer->getMesh().getOrigin(), savedState.origin, 1.0e-14);
        expectVectorNear(restoredPhysicalRMin, savedState.partrmin, 1.0e-14);
        expectVectorNear(restoredPhysicalRMax, savedState.partrmax, 1.0e-14);
    }

    // ----------------------------------------------------------------------------
    // bunchUpdate() in BeamBeam-window mode
    // ----------------------------------------------------------------------------

    // Once the BeamBeam-window configuration is active, bunchUpdate() must keep
    // the longitudinal field domain frozen while still allowing the transverse
    // field domain to follow the physical bunch. This is the key invariant behind
    // the frozen Eulerian-in-z mesh used during the BeamBeam solve.
    TEST_F(BeamBeamPartBunchTest, BunchUpdateKeepsLongitudinalFieldDomainFrozen) {
        auto bunch          = makeBunch();
        auto fieldContainer = bunch->getFieldContainer();

        setParticlePositions(bunch, {Vector3d(-0.20, -0.10, -0.05), Vector3d(0.10, 0.20, 0.05)});
        bunch->bunchUpdate();

        // Keep the test particles inside the frozen longitudinal window.  The
        // invariant under test is transverse growth with fixed z bounds; placing
        // particles outside z would exercise undefined MPI layout migration instead.
        bunch->setBeamBeamWindowConfig(0.50, 0.0, -0.25, 0.25, false);
        bunch->enableBeamBeamWindowMesh(0.0, 0.50);

        const Vector3d frozenRMin = fieldContainer->getRMin();
        const Vector3d frozenRMax = fieldContainer->getRMax();
        const Vector3d frozenHr   = fieldContainer->getHr();

        setParticlePositions(bunch, {Vector3d(-0.60, -0.30, -0.20), Vector3d(0.50, 0.40, 0.20)});
        bunch->bunchUpdate();

        const Vector3d updatedRMin = fieldContainer->getRMin();
        const Vector3d updatedRMax = fieldContainer->getRMax();
        const Vector3d updatedHr   = fieldContainer->getHr();

        EXPECT_DOUBLE_EQ(updatedRMin[2], frozenRMin[2]);
        EXPECT_DOUBLE_EQ(updatedRMax[2], frozenRMax[2]);
        EXPECT_DOUBLE_EQ(updatedHr[2], frozenHr[2]);

        EXPECT_LT(updatedRMin[0], frozenRMin[0]);
        EXPECT_GT(updatedRMax[0], frozenRMax[0]);
    }

}  // namespace
