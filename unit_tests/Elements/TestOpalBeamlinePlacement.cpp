#include "AbsBeamline/ScalingFFAMagnet.h"
#include "Algorithms/DefaultVisitor.h"
#include "Elements/PlacementResolver.h"
#include "gtest/gtest.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/MultipoleTRep.h"
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
#include "ValueDefinitions/RealVariable.h"

#include <cmath>
#include <tuple>
#include "../KeywordObjects.h"

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

        auto beam = std::make_shared<Beam>();
        Beam::defineBuiltins();

        Beam::defineBuiltins();
        Beam* opBeam = Beam::find("UNNAMED_BEAM");
        EXPECT_NE(opBeam, nullptr);

        auto bunch = std::make_shared<PartBunch_t>(
                std::vector{1.0}, std::vector{1.0}, std::vector<Beam*>{opBeam},
                std::vector<size_t>{numParticles}, 1.0, "LF2",
                opalx::spacecharge::CartesianDomainConfig3D{.periodicParticleBoundary = true});
        bunch->getParticleContainer()->createParticles(numParticles);
        return bunch;
    }

    std::shared_ptr<FieldSolverCmd> fsCmdBase_m;
    std::shared_ptr<DataSink> dataSink_m;
};

TEST_F(OpalBeamlinePlacementTest, BridgeReturnsPlacedElementViewAndAppliesTheMisalignment) {
    OpalBeamline beamline;
    auto drift = std::make_shared<DriftRep>("D2");
    drift->getGeometry().setElementLength(2.0);
    drift->setCSTrafoGlobal2Local(
            CoordinateSystemTrafo(Vector3(0.5, -1.0, 4.0), rotationAroundY(M_PI / 8.0)));
    drift->setMisalignment(CoordinateSystemTrafo(Vector3(0.25, 0.0, 0.0), Quaternion()));

    const CoordinateSystemTrafo nominal = drift->getCSTrafoGlobal2Local();
    // The local-frame helpers use the frame in which the tracker applies the element.
    const CoordinateSystemTrafo applied = drift->getMisalignment() * nominal;
    const Vector3 point(1.0, 2.0, 3.0);

    expectVectorNear(beamline.transformToLocalCS(drift, point), applied.transformTo(point));
    expectVectorNear(beamline.transformFromLocalCS(drift, point), applied.transformFrom(point));
    expectVectorNear(beamline.rotateToLocalCS(drift, point), applied.rotateTo(point));
    expectVectorNear(beamline.rotateFromLocalCS(drift, point), applied.rotateFrom(point));
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

TEST_F(OpalBeamlinePlacementTest, CurvedMultipoleTIsFramedLikeAnSBend) {
    // Placement frames a bend by its geometry, not by its element type, so a curved
    // MULTIPOLET turns the reference path exactly as an SBEND of the same arc and angle does:
    // its own entrance frame and the frame of the following element agree.
    const double edge = 0.5, arc = 0.5, angle = 0.3;

    const auto place = [&](const std::shared_ptr<ElementBase>& magnet) {
        magnet->setElementPosition(edge);

        auto drift = std::make_shared<DriftRep>("DRIFT_AFTER");
        drift->getGeometry().setElementLength(1.0);
        drift->setElementPosition(edge + arc);

        auto bunch = makeBunch(0);
        DummyBeamline beamlineForVisitor;
        DefaultVisitor visitor(beamlineForVisitor, false, false);
        auto beamline = std::make_shared<OpalBeamline>();
        beamline->visit(*magnet, visitor, *bunch);
        beamline->visit(*drift, visitor, *bunch);
        beamline->prepareSections();

        std::shared_ptr<ElementBase> placedMagnet, placedDrift;
        for (const auto& comp : beamline->getElements()) {
            if (comp->getType() == ElementType::DRIFT) {
                placedDrift = comp;
            } else {
                placedMagnet = comp;
            }
        }
        EXPECT_NE(placedMagnet, nullptr);
        EXPECT_NE(placedDrift, nullptr);
        return std::make_tuple(
                beamline, placedMagnet->getCSTrafoGlobal2Local(),
                placedDrift->getCSTrafoGlobal2Local());
    };

    auto bend                                            = std::make_shared<SBendRep>("BEND_REF");
    bend->getGeometry()                                  = Geometry::makeSBend(arc, angle / arc);
    const auto [bendBeamline, bendFrame, driftAfterBend] = place(bend);

    auto magnet = std::make_shared<MultipoleTRep>("MULTIPOLET_UNDER_TEST");
    magnet->setElementLength(arc);
    magnet->setBendAngle(angle, false);
    magnet->setTransProfile({-1.0});
    magnet->setFringeField(arc / 2, 0.01, 0.01);
    ASSERT_TRUE(magnet->getGeometry().isBend());
    const auto [magnetBeamline, magnetFrame, driftAfterMagnet] = place(magnet);

    expectTrafoNear(magnetFrame, bendFrame);
    expectTrafoNear(driftAfterMagnet, driftAfterBend);
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
    auto sbendObject = opalx::test::makeKeywordObject<OpalSBend>("SBEND", "SBEND");
    OpalSBend& sbend = *sbendObject;
    Attributes::setReal(*sbend.findAttribute("L"), 0.5);
    Attributes::setReal(*sbend.findAttribute("ELEMEDGE"), 0.0);
    Attributes::setReal(*sbend.findAttribute("ANGLE"), 0.3);
    Attributes::setReal(*sbend.findAttribute("E1"), 0.05);
    EXPECT_THROW(opalx::test::buildFrom(sbend), OpalException);
}

TEST_F(OpalBeamlinePlacementTest, RBendDerivesEntranceAngleFromE2) {
    [[maybe_unused]] static auto* p0Owner = &opalx::test::keyword<RealVariable>("REAL_VARIABLE");
    auto rbendObject = opalx::test::makeKeywordObject<OpalRBend>("RBEND", "RBEND");
    OpalRBend& rbend = *rbendObject;
    Attributes::setReal(*rbend.findAttribute("L"), 0.5);
    Attributes::setReal(*rbend.findAttribute("ELEMEDGE"), 0.0);
    Attributes::setReal(*rbend.findAttribute("ANGLE"), 0.3);
    Attributes::setReal(*rbend.findAttribute("E2"), 0.05);
    std::shared_ptr<ElementBase> built;
    ASSERT_NO_THROW(built = opalx::test::buildFrom(rbend));
    auto& bend = dynamic_cast<RBend&>(*built);
    EXPECT_FALSE(bend.hasFieldMap());
    EXPECT_DOUBLE_EQ(bend.getDesignEnergy(), 0.0);
    EXPECT_DOUBLE_EQ(bend.getGeometry().getEntranceAngle(), 0.25);
    EXPECT_DOUBLE_EQ(bend.getGeometry().getExitAngle(), 0.05);
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

// The school F/D cell has overlapping spiral fringe supports. Test placement
// independently of tracking, against the median-plane scaling law in metres.
TEST_F(OpalBeamlinePlacementTest, SchoolScalingFFACellEntrancePosesAndSignedRadii) {
    const double radius = 4.42, pi = std::acos(-1.);
    const double starts[]    = {(13.5 + 0.73) * pi / 180, (2.25 - 0.73) * pi / 180};
    const double lengths[]   = {4.5 * pi / 180 * radius, 2.25 * pi / 180 * radius};
    const double strengths[] = {-0.6532097983351404, 0.502765236428169};
    const std::vector<double> coefficients{-1.3, 3.42, -1.20, 0.23};
    const auto edge = [](double u) {
        return 1 / (1 + std::exp(-1.3 + u * (3.42 + u * (-1.20 + 0.23 * u))));
    };
    for (double sign : {1., -1.})
        for (double spiral : {pi / 4, 55 * pi / 180})
            for (bool pose : {false, true}) {
                ElementList cell;
                for (unsigned i = 0; i < 2; ++i) {
                    const std::string name = "school_fringe_" + std::to_string(i);
                    endfieldmodel::EndFieldModelManager::getEFMManager()->setEndFieldModel(
                            name, std::make_shared<endfieldmodel::AsymmetricEnge>(
                                          coefficients, lengths[i] / 2, 0.12, coefficients,
                                          lengths[i] / 2, 0.12));
                    auto magnet = std::make_shared<ScalingFFAMagnet>(name);
                    magnet->setR0(sign * radius);
                    magnet->setEndFieldName(name);
                    magnet->setDipoleConstant(strengths[i]);
                    magnet->setFieldIndex(7.717022064259446);
                    magnet->setTanDelta(std::tan(spiral));
                    magnet->setMaxOrder(4);
                    magnet->setRMin(radius - 0.8);
                    magnet->setRMax(radius + 1.2);
                    magnet->setVerticalExtent(0.5);
                    magnet->setAzimuthalExtent(2 * pi / 16);
                    magnet->initialise();
                    if (pose) {
                        magnet->setCSTrafoGlobal2Local(CoordinateSystemTrafo(
                                Vector3(sign * radius * (std::cos(starts[i]) - 1), 0,
                                        radius * std::sin(starts[i])),
                                rotationAroundY(sign * starts[i])));
                    } else {
                        magnet->setElementPosition(radius * starts[i]);
                    }
                    cell.push_back(magnet);
                }
                const CoordinateSystemTrafo labFrame(
                        Vector3(0.3, -0.2, 0.1), rotationAroundY(0.17));
                PlacementResolver::resolve(cell, labFrame);
                unsigned i = 0;
                for (const auto& element : cell) {
                    auto& magnet = dynamic_cast<ScalingFFAMagnet&>(*element);
                    EXPECT_NEAR(magnet.getGeometry().getArcLength(), lengths[i], tol);
                    EXPECT_NEAR(magnet.getPhiStart(), lengths[i] / (2 * radius), tol);
                    for (double r : {radius - 0.1, radius, radius + 0.1})
                        for (double offset : {-0.03, 0., 0.03}) {
                            const double phi = starts[i] + lengths[i] / (2 * radius)
                                               + std::tan(spiral) * std::log(r / radius) + offset;
                            Vector3 lab(sign * (r * std::cos(phi) - radius), 0, r * std::sin(phi));
                            const auto local = element->getCSTrafoGlobal2Local().transformTo(
                                    labFrame.transformFrom(lab));
                            Vector3 b;
                            magnet.getFieldValue(local, b);
                            const double s = radius * offset, half = lengths[i] / 2;
                            const double expected =
                                    strengths[i] * std::pow(r / radius, 7.717022064259446)
                                    * (-1 + edge((-s - half) / 0.12) + edge((s - half) / 0.12));
                            EXPECT_NEAR(b(1), expected, 2e-12);
                            EXPECT_NEAR(b(0), 0, tol);
                            EXPECT_NEAR(b(2), 0, tol);
                        }
                    // Maxwell's curl relation at the median plane supplies an
                    // independent first-order off-plane check for this profile.
                    const double centrePhi = starts[i] + lengths[i] / (2 * radius);
                    const Vector3 centre(
                            sign * radius * (std::cos(centrePhi) - 1), 0,
                            radius * std::sin(centrePhi));
                    const auto point = element->getCSTrafoGlobal2Local().transformTo(
                            labFrame.transformFrom(centre));
                    const double h = 1e-6;
                    Vector3 off    = point, offField;
                    off(1)         = h;
                    magnet.getFieldValue(off, offField);
                    for (unsigned d : {0u, 2u}) {
                        Vector3 plus = point, minus = point, bp, bm;
                        plus(d) += h;
                        minus(d) -= h;
                        magnet.getFieldValue(plus, bp);
                        magnet.getFieldValue(minus, bm);
                        EXPECT_NEAR(offField(d) / h, (bp(1) - bm(1)) / (2 * h), 2e-7);
                    }
                    if (i == 1) {
                        magnet.setDipoleConstant(0.);
                        magnet.getFieldValue(off, offField);
                        expectVectorNear(offField, Vector3(0));
                    }
                    ++i;
                }
            }
    endfieldmodel::EndFieldModelManager::clearEFMManager();
}

TEST_F(OpalBeamlinePlacementTest, ScalingFFAExitFramesPlaceFollowingDrift) {
    const std::string name = "exit_frame_fringe";
    endfieldmodel::EndFieldModelManager::getEFMManager()->setEndFieldModel(
            name, std::make_shared<endfieldmodel::Tanh>(0.2, 0.02, 4));
    const CoordinateSystemTrafo labFrame(Vector3(0.3, -0.2, 0.1), rotationAroundY(0.17));
    for (double sign : {1., -1.}) {
        auto magnet = std::make_shared<ScalingFFAMagnet>("ffa");
        magnet->setR0(sign * 4.42);
        magnet->setEndFieldName(name);
        magnet->setupEndField();
        magnet->setElementPosition(0.7);
        auto drift = std::make_shared<DriftRep>("following");
        drift->getGeometry().setElementLength(0.1);
        drift->setElementPosition(0.7 + 0.4 + 0.2);
        ElementList elements{magnet, drift};
        PlacementResolver::resolve(elements, labFrame);
        const double angle = 1.1 / 4.42;
        const Vector3 exit(sign * 4.42 * (std::cos(angle) - 1), 0, 4.42 * std::sin(angle));
        const Vector3 tangent(-sign * std::sin(angle), 0, std::cos(angle));
        expectVectorNear(
                drift->getCSTrafoGlobal2Local().transformTo(
                        labFrame.transformFrom(Vector3(exit + 0.2 * tangent))),
                Vector3(0));
        expectVectorNear(
                drift->getCSTrafoGlobal2Local().rotateTo(labFrame.rotateFrom(tangent)),
                Vector3(0, 0, 1));
    }
    endfieldmodel::EndFieldModelManager::clearEFMManager();
}
