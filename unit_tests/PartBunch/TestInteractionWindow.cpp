#include <gtest/gtest.h>

#include "Ippl.h"
#include "Utility/Inform.h"

#include "Attributes/Attributes.h"
#include "PartBunch/PartBunch.h"
#include "Structure/FieldSolverCmd.h"

#include <memory>
#include <string>

extern Inform* gmsg;

namespace {

using Bunch_t  = PartBunch<double, 3>;
using Vector3d = ippl::Vector<double, 3>;

/**
 * @brief Test-only FieldSolverCmd helper exposing the minimal attribute setters
 * needed to construct a real PartBunch.
 *
 * The interaction-window unit tests intentionally use the real PartBunch type,
 * because the bugs we fixed were caused by the coupling between:
 * - cached physical bunch bounds,
 * - field-domain bounds, and
 * - the temporary interaction-window mesh replacement.
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

std::shared_ptr<FieldSolverCmd> makeFieldSolverCmd() {
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
    auto fieldSolver = makeFieldSolverCmd();
    auto dataSink    = std::shared_ptr<DataSink>();

    return std::make_shared<Bunch_t>(
        -1.0e-15, 5.10999e-4, 16, 1.0, "LF2", fieldSolver, dataSink);
}

void expectVectorNear(
    const Vector3d& actual,
    const Vector3d& expected,
    double tolerance) {
    for (unsigned d = 0; d < 3; ++d) {
        EXPECT_NEAR(actual[d], expected[d], tolerance) << "Mismatch in component " << d;
    }
}

class InteractionWindowPartBunchTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        // PartBunch teardown writes through the global OPAL Inform stream.
        // The regular `tests/` binary creates this in its custom main(), but
        // the unit_tests tree uses gtest_main, so this fixture must provide
        // the minimal global setup needed by a standalone PartBunch instance.
        gmsg = new Inform("UnitTests: ", std::cerr);
    }

    static void TearDownTestSuite() {
        delete gmsg;
        gmsg = nullptr;
        ippl::finalize();
    }
};

// ----------------------------------------------------------------------------
// InteractionWindowConfig lifecycle
// ----------------------------------------------------------------------------

// Verifies that the bunch-side interaction-window configuration behaves like a
// narrow geometry input only: absent by default, present after configuration,
// and removable again without affecting unrelated state.
TEST_F(InteractionWindowPartBunchTest, InteractionWindowConfigLifecycle) {
    auto bunch = makeBunch();

    ASSERT_FALSE(bunch->hasInteractionWindowConfig());

    bunch->setInteractionWindowConfig(1.25, 0.50);
    ASSERT_TRUE(bunch->hasInteractionWindowConfig());

    const auto& config = bunch->getInteractionWindowConfig();
    EXPECT_DOUBLE_EQ(config.interactionPointLocalZ, 1.25);
    EXPECT_DOUBLE_EQ(config.interactionWindowLength, 0.50);

    bunch->clearInteractionWindowConfig();
    EXPECT_FALSE(bunch->hasInteractionWindowConfig());
}

// ----------------------------------------------------------------------------
// enableInteractionWindowMesh
// ----------------------------------------------------------------------------

// The interaction-window mesh switch should keep the transverse field domain
// unchanged while replacing only the longitudinal mesh extent. This is the
// explicit Lagrangian-to-Eulerian transition in z that the current model relies on.
TEST_F(InteractionWindowPartBunchTest, EnableInteractionWindowMeshOnlyChangesLongitudinalDomain) {
    auto bunch = makeBunch();
    auto fieldContainer = bunch->getFieldContainer();

    const Vector3d initialRMin = fieldContainer->getRMin();
    const Vector3d initialRMax = fieldContainer->getRMax();
    const Vector3d initialHr   = fieldContainer->getHr();

    bunch->enableInteractionWindowMesh(1.25, 0.50);

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

// ----------------------------------------------------------------------------
// save/restore field-domain state
// ----------------------------------------------------------------------------

// The field-domain state roundtrip must restore both:
// - the mesh-aligned field-domain quantities, and
// - the cached physical bunch bounds returned by get_bounds().
//
// This test covers the exact contract that the interaction-window mode relies on
// when it temporarily swaps in an Eulerian longitudinal mesh and later returns to
// the original bunch-following state.
TEST_F(InteractionWindowPartBunchTest, RestoreFieldDomainStateRestoresMeshAndPhysicalBounds) {
    auto bunch = makeBunch();
    auto fieldContainer = bunch->getFieldContainer();

    const Vector3d physicalRMin(-0.2, -0.3, -0.4);
    const Vector3d physicalRMax(0.5, 0.6, 0.7);
    bunch->setPhysicalBounds(physicalRMin, physicalRMax);

    const auto savedState = bunch->saveFieldDomainState();
    bunch->enableInteractionWindowMesh(1.25, 0.50);
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

}  // namespace
