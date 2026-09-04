#include "Algorithms/DefaultVisitor.h"
#include "gtest/gtest.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/RBendRep.h"
#include "BeamlineCore/SBendRep.h"
#include "BeamlineGeometry/Geometry.h"
#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"
#include "Elements/OpalRBend.h"
#include "Elements/OpalSBend.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

#include <cmath>

extern Inform* gmsg;

namespace {
    constexpr double tol = 1e-12;

    using Vector3 = ippl::Vector<double, 3>;

    Quaternion rotationAroundY(double angle) {
        return Quaternion(std::cos(0.5 * angle), 0.0, std::sin(0.5 * angle), 0.0);
    }

    void expectVectorNear(const Vector3& actual, const Vector3& expected) {
        EXPECT_NEAR(actual(0), expected(0), tol);
        EXPECT_NEAR(actual(1), expected(1), tol);
        EXPECT_NEAR(actual(2), expected(2), tol);
    }

    // Compare two frames fully: origin plus every rotation-matrix entry.
    void expectTrafoNear(
            const CoordinateSystemTrafo& actual, const CoordinateSystemTrafo& expected) {
        expectVectorNear(actual.getOrigin(), expected.getOrigin());
        const auto& a = actual.getRotationMatrix();
        const auto& e = expected.getRotationMatrix();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                EXPECT_NEAR(a(i, j), e(i, j), 1e-9);
            }
        }
    }

    class DummyBeamline final : public Beamline {
    public:
        DummyBeamline() : Beamline("dummy") {}

        ElementType getType() const override { return ElementType::BEAMLINE; }
        Geometry& getGeometry() override { return geometry_; }
        const Geometry& getGeometry() const override { return geometry_; }
        void accept(BeamlineVisitor& visitor) const override { visitor.visitBeamline(*this); }
        ElementBase* clone() const override { return new DummyBeamline(*this); }
        void iterate(BeamlineVisitor&, bool) const override {}

    private:
        Geometry geometry_{Geometry::makeNull()};
    };
}  // namespace

class OpalBeamlinePlacementTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;

        ippl::initialize(argc, argv);
        gmsg                = new Inform(nullptr, -1);
        Options::enableHDF5 = false;
    }

    static void TearDownTestSuite() {
        delete gmsg;
        gmsg = nullptr;
        ippl::finalize();
    }

    void SetUp() override {
        OpalData::getInstance()->storeInputFn("TestOpalBeamlinePlacement.opal");
        OpalData::getInstance()->setOpenMode(OpalData::OpenMode::WRITE);
    }

    class TestableFieldSolverCmd : public FieldSolverCmd {
    public:
        void setType(const std::string& t) {
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::TYPE], t);
        }

        void setBCX(const std::string& bc) {
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTX], bc);
        }

        void setBCY(const std::string& bc) {
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTY], bc);
        }

        void setBCZ(const std::string& bc) {
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTZ], bc);
        }
    };

    std::shared_ptr<PartBunch_t> makeBunch(const size_t numParticles) {
        dataSink_m       = std::make_shared<DataSink>();
        const auto fsCmd = std::make_shared<TestableFieldSolverCmd>();
        fsCmdBase_m      = fsCmd;
        fsCmd->setType("NONE");
        fsCmd->setNX(8);
        fsCmd->setNY(8);
        fsCmd->setNZ(8);
        fsCmd->setBCX("PERIODIC");
        fsCmd->setBCY("PERIODIC");
        fsCmd->setBCZ("PERIODIC");

        auto beam    = std::make_shared<Beam>();
        Beam* opBeam = Beam::find("UNNAMED_BEAM");
        EXPECT_NE(opBeam, nullptr);

        auto bunch = std::make_shared<PartBunch_t>(
                std::vector{1.0}, std::vector{1.0}, std::vector<Beam*>{opBeam},
                std::vector<size_t>{numParticles}, 1.0, "LF2", fsCmdBase_m.get(), dataSink_m.get());
        bunch->getParticleContainer()->createParticles(numParticles);
        return bunch;
    }

    std::shared_ptr<FieldSolverCmd> fsCmdBase_m;
    std::shared_ptr<DataSink> dataSink_m;
};

TEST_F(OpalBeamlinePlacementTest, BridgeReturnsPlacedElementViewAndPreservesNominalQueries) {
    OpalBeamline beamline;
    auto drift = std::make_shared<DriftRep>("D2");
    drift->getGeometry().setElementLength(2.0);
    drift->setCSTrafoGlobal2Local(
            CoordinateSystemTrafo(Vector3(0.5, -1.0, 4.0), rotationAroundY(M_PI / 8.0)));
    drift->setMisalignment(CoordinateSystemTrafo(Vector3(0.25, 0.0, 0.0), Quaternion()));

    const CoordinateSystemTrafo nominal = drift->getCSTrafoGlobal2Local();
    const Vector3 point(1.0, 2.0, 3.0);

    expectVectorNear(beamline.transformToLocalCS(drift, point), nominal.transformTo(point));
    expectVectorNear(beamline.transformFromLocalCS(drift, point), nominal.transformFrom(point));
    expectVectorNear(beamline.rotateToLocalCS(drift, point), nominal.rotateTo(point));
    expectVectorNear(beamline.rotateFromLocalCS(drift, point), nominal.rotateFrom(point));
    expectVectorNear(beamline.getCSTrafoLab2Local(drift).getOrigin(), nominal.getOrigin());
    // entry edge is identity for a straight element ⇒ entry origin equals nominal origin
    expectVectorNear(beamline.getNominalEntryTransform(drift).getOrigin(), nominal.getOrigin());
    expectVectorNear(
            beamline.getNominalExitTransform(drift).getOrigin(),
            (drift->getGeometry().getEdgeToEnd() * nominal).getOrigin());
    expectVectorNear(beamline.getMisalignment(drift).getOrigin(), Vector3(0.25, 0.0, 0.0));
}

TEST_F(OpalBeamlinePlacementTest, PrepareSectionsComposes6DPoseWithLabFrame) {
    // Mode A (6D pose): the element records its global-to-local pose but stays unpositioned;
    // prepareSections() composes it with the lab frame in one place (no ELEMEDGE set).
    OpalBeamline beamline(Vector3(0.0, 0.0, 5.0), Quaternion());
    DriftRep drift("D3");
    drift.getGeometry().setElementLength(0.4);
    drift.setCSTrafoGlobal2Local(
            CoordinateSystemTrafo(Vector3(1.0, 2.0, 3.0), rotationAroundY(M_PI / 10.0)));

    CoordinateSystemTrafo expected = drift.getCSTrafoGlobal2Local();
    expected *= beamline.getCSTrafoLab2Local();

    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);
    beamline.visit(drift, visitor, *bunch);
    beamline.prepareSections();

    const auto component = *beamline.getElements().begin();
    expectVectorNear(beamline.getCSTrafoLab2Local(component).getOrigin(), expected.getOrigin());
}

TEST_F(OpalBeamlinePlacementTest, PrepareSectionsCompilesElementPositionIntoNominalPlacement) {
    DriftRep drift("D4");
    drift.getGeometry().setElementLength(0.4);
    drift.setElementPosition(1.25);

    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);
    OpalBeamline beamline;
    beamline.visit(drift, visitor, *bunch);

    beamline.prepareSections();

    const auto elements = beamline.getElements();
    ASSERT_EQ(elements.size(), 1u);
    const auto component = *elements.begin();

    expectVectorNear(beamline.getCSTrafoLab2Local(component).getOrigin(), Vector3(0.0, 0.0, 1.25));
    expectVectorNear(
            beamline.getNominalEntryTransform(component).getOrigin(), Vector3(0.0, 0.0, 1.25));
    expectVectorNear(
            beamline.getNominalExitTransform(component).getOrigin(), Vector3(0.0, 0.0, 1.65));
}

TEST_F(OpalBeamlinePlacementTest, BeamlineOwnsPlacedElementAssemblySnapshot) {
    DriftRep drift("D5");
    drift.getGeometry().setElementLength(0.3);
    drift.setElementPosition(0.75);

    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);
    OpalBeamline beamline;
    beamline.visit(drift, visitor, *bunch);
    beamline.prepareSections();

    const auto elements = beamline.getElements();
    ASSERT_EQ(elements.size(), 1u);
    const auto component = *elements.begin();

    const Vector3 assembledOrigin = beamline.getCSTrafoLab2Local(component).getOrigin();
    expectVectorNear(assembledOrigin, Vector3(0.0, 0.0, 0.75));

    component->setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector3(9.0, 8.0, 7.0), Quaternion()));

    expectVectorNear(beamline.getCSTrafoLab2Local(component).getOrigin(), assembledOrigin);
}

TEST_F(OpalBeamlinePlacementTest, SequentialBendThenDriftChainWithoutGap) {
    // Regression: a bend placed by ELEMEDGE must advance the running path length by its arc
    // length, so a drift placed at ELEMEDGE = bendEdge + arc sits flush against the bend exit
    // (its entrance origin equals the bend's exit origin). PlacementResolver previously measured
    // the entrance hard-edge offset from the geometry origin (~half the chord), which drove the
    // bend's path-length contribution to ~0 and inserted a spurious ~chord-length field-free gap
    // after every sequentially placed bend, stretching multi-bend lines (e.g. a chicane).
    const double bendEdge = 0.5, arc = 0.5, angle = 0.3;

    auto bend           = std::make_shared<SBendRep>("BEND_B");
    bend->getGeometry() = Geometry::makeSBend(arc, angle / arc);
    bend->setElementPosition(bendEdge);

    auto drift = std::make_shared<DriftRep>("DRIFT_D");
    drift->getGeometry().setElementLength(1.0);
    drift->setElementPosition(bendEdge + arc);

    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);
    OpalBeamline beamline;
    beamline.visit(*bend, visitor, *bunch);
    beamline.visit(*drift, visitor, *bunch);
    beamline.prepareSections();

    std::shared_ptr<ElementBase> placedBend, placedDrift;
    for (const auto& comp : beamline.getElements()) {
        if (comp->getType() == ElementType::SBEND) {
            placedBend = comp;
        } else if (comp->getType() == ElementType::DRIFT) {
            placedDrift = comp;
        }
    }
    ASSERT_NE(placedBend, nullptr);
    ASSERT_NE(placedDrift, nullptr);

    // The bend is the first element (lab frame at the origin), so its entrance sits at
    // z = bendEdge. A contiguous drift starts exactly one chord (the bend's straight
    // entrance-to-exit extent) downstream. With the bug the bend contributed ~0 path
    // length, so the drift was pushed a further ~chord away.
    const Vector3 bendEntry(0.0, 0.0, bendEdge);
    const Vector3 driftEntry = beamline.getNominalEntryTransform(placedDrift).getOrigin();
    const Vector3 gapVector  = driftEntry - bendEntry;
    EXPECT_NEAR(euclidean_norm(gapVector), bend->getGeometry().getChordLength(), 1e-9);
}

TEST_F(OpalBeamlinePlacementTest, ModeAStoresEntranceFrameVerbatim) {
    // Mode A (X/Y/Z/PHI/PSI/THETA) IS the element's geometrical ENTRANCE frame — the frame the
    // tracker and field kernels use, same as ELEMEDGE placement. It is stored verbatim (here the
    // lab frame is the identity). Feed the entrance frame (0,0,s0) tangent +z and check it is
    // stored unchanged. (The old code treated X/Y/Z as the body centre and shifted it by
    // entranceFromBodyCentre(), so this same input would have been mis-placed by ~half the arc.)
    const double s0 = 1.0, L = 1.0, angle = 0.3, h = angle / L;
    const CoordinateSystemTrafo entrancePose(Vector3(0.0, 0.0, s0), Quaternion());

    auto bend           = std::make_shared<SBendRep>("MA_BEND");
    bend->getGeometry() = Geometry::makeSBend(L, h);
    bend->setCSTrafoGlobal2Local(entrancePose);  // Mode A: no ELEMEDGE

    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);
    OpalBeamline beamline;
    beamline.visit(*bend, visitor, *bunch);
    beamline.prepareSections();

    const auto placed = *beamline.getElements().begin();
    expectTrafoNear(beamline.getCSTrafoLab2Local(placed), entrancePose);
}

// Place a bend by ELEMEDGE (Mode B), read the entrance frame PlacementResolver computes, then
// place an identical bend by Mode A using that frame. Both must store the same frame — Mode A is
// now a verbatim pass-through of the entrance frame Mode B builds. This is the regression guard:
// the old entranceFromBodyCentre() shift made Mode A diverge from Mode B for bends.
TEST_F(OpalBeamlinePlacementTest, SBendModeAMatchesModeB) {
    const double edge = 0.7, arc = 0.5, angle = 0.4, h = angle / arc;

    // Mode B.
    auto bendB           = std::make_shared<SBendRep>("SB_MODEB");
    bendB->getGeometry() = Geometry::makeSBend(arc, h);
    bendB->setElementPosition(edge);
    auto bunchB = makeBunch(0);
    DummyBeamline dbB;
    DefaultVisitor visB(dbB, false, false);
    OpalBeamline beamlineB;
    beamlineB.visit(*bendB, visB, *bunchB);
    beamlineB.prepareSections();
    const CoordinateSystemTrafo frameB =
            beamlineB.getCSTrafoLab2Local(*beamlineB.getElements().begin());

    // Mode A with the pose equal to Mode B's entrance frame (identity lab frame ⇒ verbatim).
    auto bendA           = std::make_shared<SBendRep>("SB_MODEA");
    bendA->getGeometry() = Geometry::makeSBend(arc, h);
    bendA->setCSTrafoGlobal2Local(frameB);
    auto bunchA = makeBunch(0);
    DummyBeamline dbA;
    DefaultVisitor visA(dbA, false, false);
    OpalBeamline beamlineA;
    beamlineA.visit(*bendA, visA, *bunchA);
    beamlineA.prepareSections();
    const CoordinateSystemTrafo frameA =
            beamlineA.getCSTrafoLab2Local(*beamlineA.getElements().begin());

    expectTrafoNear(frameA, frameB);
}

TEST_F(OpalBeamlinePlacementTest, RBendModeAMatchesModeB) {
    // RBend stores the box/chord frame (rotated half the bend angle off the entrance-orbit
    // tangent), so frameB has a non-trivial rotation — a stronger check than the SBend case.
    const double edge = 0.7, boxLength = 0.5, angle = 0.4;

    auto bendB           = std::make_shared<RBendRep>("RB_MODEB");
    bendB->getGeometry() = Geometry::makeRBend(boxLength, angle);
    bendB->setElementPosition(edge);
    auto bunchB = makeBunch(0);
    DummyBeamline dbB;
    DefaultVisitor visB(dbB, false, false);
    OpalBeamline beamlineB;
    beamlineB.visit(*bendB, visB, *bunchB);
    beamlineB.prepareSections();
    const CoordinateSystemTrafo frameB =
            beamlineB.getCSTrafoLab2Local(*beamlineB.getElements().begin());

    auto bendA           = std::make_shared<RBendRep>("RB_MODEA");
    bendA->getGeometry() = Geometry::makeRBend(boxLength, angle);
    bendA->setCSTrafoGlobal2Local(frameB);
    auto bunchA = makeBunch(0);
    DummyBeamline dbA;
    DefaultVisitor visA(dbA, false, false);
    OpalBeamline beamlineA;
    beamlineA.visit(*bendA, visA, *bunchA);
    beamlineA.prepareSections();
    const CoordinateSystemTrafo frameA =
            beamlineA.getCSTrafoLab2Local(*beamlineA.getElements().begin());

    // Sanity: the RBend frame really is rotated (not identity).
    EXPECT_GT(std::abs(frameB.getRotationMatrix()(0, 2)), 1e-3);
    expectTrafoNear(frameA, frameB);
}

TEST_F(OpalBeamlinePlacementTest, SBendRejectsPoleFaceAngleE1) {
    // E1/E2 pole-face rotations are not wired yet; specifying one must error instead of being
    // silently ignored.
    OpalSBend sbend;
    Attributes::setReal(*sbend.findAttribute("L"), 0.5);
    Attributes::setReal(*sbend.findAttribute("ELEMEDGE"), 0.0);
    Attributes::setReal(*sbend.findAttribute("ANGLE"), 0.3);
    Attributes::setReal(*sbend.findAttribute("E1"), 0.05);
    EXPECT_THROW(sbend.update(), OpalException);
}

TEST_F(OpalBeamlinePlacementTest, RBendRejectsPoleFaceAngleE2) {
    OpalRBend rbend;
    Attributes::setReal(*rbend.findAttribute("L"), 0.5);
    Attributes::setReal(*rbend.findAttribute("ELEMEDGE"), 0.0);
    Attributes::setReal(*rbend.findAttribute("ANGLE"), 0.3);
    Attributes::setReal(*rbend.findAttribute("E2"), 0.05);
    EXPECT_THROW(rbend.update(), OpalException);
}

TEST_F(OpalBeamlinePlacementTest, RejectsMixedPlacementConventions) {
    // A beamline must use ONE placement convention. A 6D-posed (Mode A) bend does not deflect
    // the reference orbit that ELEMEDGE (Mode B) elements are walked along, so mixing the two
    // silently misplaces every ELEMEDGE element downstream of the posed bend. resolve() must
    // reject the mix instead of producing a geometrically inconsistent lattice.
    auto poseBend           = std::make_shared<SBendRep>("MIX_BEND_A");
    poseBend->getGeometry() = Geometry::makeSBend(0.5, 0.3 / 0.5);
    poseBend->setCSTrafoGlobal2Local(
            CoordinateSystemTrafo(Vector3(0.0, 0.0, 1.0), Quaternion()));  // Mode A: 6D pose

    auto edgeDrift = std::make_shared<DriftRep>("MIX_DRIFT_B");
    edgeDrift->getGeometry().setElementLength(0.4);
    edgeDrift->setElementPosition(2.0);  // Mode B: ELEMEDGE

    auto bunch = makeBunch(0);
    DummyBeamline db;
    DefaultVisitor visitor(db, false, false);
    OpalBeamline beamline;
    beamline.visit(*poseBend, visitor, *bunch);
    beamline.visit(*edgeDrift, visitor, *bunch);

    EXPECT_THROW(beamline.prepareSections(), OpalException);
}

TEST_F(OpalBeamlinePlacementTest, AcceptsUniform6DPoseAcrossMultipleElements) {
    // The mixed-convention guard must NOT fire when every element uses the same convention:
    // two Mode-A (6D pose) elements are a valid, uniform beamline.
    auto d1 = std::make_shared<DriftRep>("POSE_D1");
    d1->getGeometry().setElementLength(0.4);
    d1->setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector3(0.0, 0.0, 1.0), Quaternion()));

    auto d2 = std::make_shared<DriftRep>("POSE_D2");
    d2->getGeometry().setElementLength(0.4);
    d2->setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector3(0.0, 0.0, 2.0), Quaternion()));

    auto bunch = makeBunch(0);
    DummyBeamline db;
    DefaultVisitor visitor(db, false, false);
    OpalBeamline beamline;
    beamline.visit(*d1, visitor, *bunch);
    beamline.visit(*d2, visitor, *bunch);

    EXPECT_NO_THROW(beamline.prepareSections());
}
