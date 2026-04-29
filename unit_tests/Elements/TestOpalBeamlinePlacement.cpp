#include "Algorithms/DefaultVisitor.h"
#include "gtest/gtest.h"

#include "AbsBeamline/BendBase.h"
#include "AbstractObjects/OpalData.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/RBendRep.h"
#include "BeamlineCore/SBendRep.h"
#include "BeamlineGeometry/NullGeometry.h"
#include "BeamlineGeometry/PlacementPose.h"
#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"
#include "Physics/Physics.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Structure/MeshGenerator.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>

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

    std::filesystem::path outputPath(const std::string& suffix) {
        return std::filesystem::path(OpalData::getInstance()->getAuxiliaryOutputDirectory())
               / ("TestOpalBeamlinePlacement" + suffix);
    }

    class DummyBeamline final : public Beamline {
    public:
        DummyBeamline() : Beamline("dummy") {}

        ElementType getType() const override { return ElementType::BEAMLINE; }
        BGeometryBase& getGeometry() override { return geometry_; }
        const BGeometryBase& getGeometry() const override { return geometry_; }
        void accept(BeamlineVisitor& visitor) const override { visitor.visitBeamline(*this); }
        ElementBase* clone() const override { return new DummyBeamline(*this); }
        void iterate(BeamlineVisitor&, bool) const override {}

    private:
        NullGeometry geometry_;
    };

    class EmptyRBend3D final : public ElementBase {
    public:
        explicit EmptyRBend3D(const std::string& name) : ElementBase(name) {
            setPlacementPose(
                    PlacementPose(CoordinateSystemTrafo(Vector3(0.0, 0.0, 0.0), Quaternion())));
        }

        ElementType getType() const override { return ElementType::RBEND3D; }
        BGeometryBase& getGeometry() override { return geometry_; }
        const BGeometryBase& getGeometry() const override { return geometry_; }
        void accept(BeamlineVisitor&) const override {}
        ElementBase* clone() const override { return new EmptyRBend3D(*this); }

    private:
        NullGeometry geometry_;
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
        std::filesystem::create_directories(OpalData::getInstance()->getAuxiliaryOutputDirectory());
        std::filesystem::remove(outputPath("_ElementPositions.py"));
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
        bunch->getParticleContainer()->create(numParticles);
        return bunch;
    }

    std::shared_ptr<FieldSolverCmd> fsCmdBase_m;
    std::shared_ptr<DataSink> dataSink_m;
};

TEST_F(OpalBeamlinePlacementTest, BridgeReturnsPlacedElementViewAndPreservesNominalQueries) {
    OpalBeamline beamline;
    auto drift = std::make_shared<DriftRep>("D2");
    drift->setElementLength(2.0);
    drift->setCSTrafoGlobal2Local(
            CoordinateSystemTrafo(Vector3(0.5, -1.0, 4.0), rotationAroundY(M_PI / 8.0)));
    drift->setMisalignment(CoordinateSystemTrafo(Vector3(0.25, 0.0, 0.0), Quaternion()));

    const PlacedElement placed = beamline.getPlacedElement(drift);
    const Vector3 point(1.0, 2.0, 3.0);

    expectVectorNear(
            beamline.transformToLocalCS(drift, point),
            placed.getNominalBodyTransform().transformTo(point));
    expectVectorNear(
            beamline.transformFromLocalCS(drift, point),
            placed.getNominalBodyTransform().transformFrom(point));
    expectVectorNear(
            beamline.rotateToLocalCS(drift, point),
            placed.getNominalBodyTransform().rotateTo(point));
    expectVectorNear(
            beamline.rotateFromLocalCS(drift, point),
            placed.getNominalBodyTransform().rotateFrom(point));
    expectVectorNear(
            beamline.getCSTrafoLab2Local(drift).getOrigin(),
            placed.getNominalBodyTransform().getOrigin());
    expectVectorNear(
            beamline.getNominalEntryTransform(drift).getOrigin(),
            placed.getNominalEntryTransform().getOrigin());
    expectVectorNear(
            beamline.getNominalExitTransform(drift).getOrigin(),
            placed.getNominalExitTransform().getOrigin());
    expectVectorNear(
            beamline.getMisalignment(drift).getOrigin(),
            placed.getMisalignment().getNominalToActual().getOrigin());
}

TEST_F(OpalBeamlinePlacementTest, PositionElementRelativeUsesPlacementPoseBridge) {
    OpalBeamline beamline(Vector3(0.0, 0.0, 5.0), Quaternion());
    auto drift = std::make_shared<DriftRep>("D3");
    drift->setPlacementPose(PlacementPose(
            CoordinateSystemTrafo(Vector3(1.0, 2.0, 3.0), rotationAroundY(M_PI / 10.0))));
    drift->fixPosition();

    CoordinateSystemTrafo expected = drift->getPlacementPose().getParentToNominal();
    expected *= beamline.getCSTrafoLab2Local();

    beamline.positionElementRelative(drift);

    expectVectorNear(
            beamline.getPlacedElement(drift).getNominalBodyTransform().getOrigin(),
            expected.getOrigin());
}

TEST_F(OpalBeamlinePlacementTest, PrepareSectionsCompilesElementPositionIntoNominalPlacement) {
    DriftRep drift("D4");
    drift.setElementLength(0.4);
    drift.setElementPosition(1.25);

    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    OpalBeamline beamline;
    beamline.visit(drift, visitor, *bunch);
    beamline.prepareSections();

    const auto elements = beamline.getElements();
    ASSERT_EQ(elements.size(), 1u);
    const auto placed = beamline.getPlacedElement(*elements.begin());

    expectVectorNear(placed.getNominalBodyTransform().getOrigin(), Vector3(0.0, 0.0, 1.25));
    expectVectorNear(placed.getNominalEntryTransform().getOrigin(), Vector3(0.0, 0.0, 1.25));
    expectVectorNear(placed.getNominalExitTransform().getOrigin(), Vector3(0.0, 0.0, 1.65));
}

TEST_F(OpalBeamlinePlacementTest, ExplicitZPositionedDriftDoesNotRequireElementEdge) {
    DriftRep drift("D4Z");
    drift.setElementLength(0.4);
    drift.setPlacementPose(
            PlacementPose(CoordinateSystemTrafo(Vector3(1.0, 2.0, 3.5), Quaternion())));
    drift.fixPosition();

    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    OpalBeamline beamline;
    beamline.visit(drift, visitor, *bunch);
    beamline.prepareSections();

    const auto elements = beamline.getElements();
    ASSERT_EQ(elements.size(), 1u);
    const auto component = *elements.begin();

    double fieldBegin = 0.0;
    double fieldEnd   = 0.0;
    component->getFieldExtend(fieldBegin, fieldEnd);

    EXPECT_NEAR(fieldBegin, 3.5, tol);
    EXPECT_NEAR(fieldEnd, 3.9, tol);
    expectVectorNear(
            beamline.getPlacedElement(component).getNominalBodyTransform().getOrigin(),
            Vector3(1.0, 2.0, 3.5));
}

TEST_F(OpalBeamlinePlacementTest, ExplicitZPositionedAnalyticBendsDoNotRequireElementEdge) {
    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    SBendRep sbend("SBEND_NO_EDGE");
    sbend.getGeometry() = PlanarArcGeometry(0.5, Physics::pi / 8.0);
    sbend.setElementLength(0.5);
    sbend.setBendAngle(Physics::pi / 8.0);
    sbend.setPlacementPose(
            PlacementPose(CoordinateSystemTrafo(Vector3(0.0, 0.0, 2.0), Quaternion())));
    sbend.fixPosition();

    OpalBeamline sbendBeamline;
    EXPECT_NO_THROW(sbendBeamline.visit(sbend, visitor, *bunch));
    {
        const auto elements = sbendBeamline.getElements();
        ASSERT_EQ(elements.size(), 1u);
        expectVectorNear(
                sbendBeamline.getPlacedElement(*elements.begin())
                        .getNominalBodyTransform()
                        .getOrigin(),
                Vector3(0.0, 0.0, 2.0));
    }

    RBendRep rbend("RBEND_NO_EDGE");
    rbend.getGeometry().setElementLength(0.5);
    rbend.getGeometry().setBendAngle(Physics::pi / 8.0);
    rbend.setElementLength(0.5);
    rbend.setBendAngle(Physics::pi / 8.0);
    rbend.setPlacementPose(
            PlacementPose(CoordinateSystemTrafo(Vector3(0.0, 0.0, 2.0), Quaternion())));
    rbend.fixPosition();

    OpalBeamline rbendBeamline;
    EXPECT_NO_THROW(rbendBeamline.visit(rbend, visitor, *bunch));
    {
        const auto elements = rbendBeamline.getElements();
        ASSERT_EQ(elements.size(), 1u);
        expectVectorNear(
                rbendBeamline.getPlacedElement(*elements.begin())
                        .getNominalBodyTransform()
                        .getOrigin(),
                Vector3(0.0, 0.0, 2.0));
    }
}

TEST_F(OpalBeamlinePlacementTest, UnpositionedAnalyticBendsStillRequireElementEdge) {
    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    SBendRep sbend("SBEND_NO_POSITION");
    sbend.getGeometry() = PlanarArcGeometry(0.5, Physics::pi / 8.0);
    sbend.setElementLength(0.5);
    sbend.setBendAngle(Physics::pi / 8.0);

    OpalBeamline sbendBeamline;
    EXPECT_THROW(sbendBeamline.visit(sbend, visitor, *bunch), OpalException);
}

TEST_F(OpalBeamlinePlacementTest, BeamlineOwnsPlacedElementAssemblySnapshot) {
    DriftRep drift("D5");
    drift.setElementLength(0.3);
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

    const Vector3 assembledOrigin =
            beamline.getPlacedElement(component).getNominalBodyTransform().getOrigin();
    expectVectorNear(assembledOrigin, Vector3(0.0, 0.0, 0.75));

    component->setPlacementPose(
            PlacementPose(CoordinateSystemTrafo(Vector3(9.0, 8.0, 7.0), Quaternion())));

    expectVectorNear(
            beamline.getPlacedElement(component).getNominalBodyTransform().getOrigin(),
            assembledOrigin);
}

TEST_F(OpalBeamlinePlacementTest, Save3DInputNormalizesExplicitBodyPoseAndAngles) {
    std::ofstream input("TestOpalBeamlinePlacement.opal");
    ASSERT_TRUE(input.is_open());
    input << "D6: Drift, L = 0.3, ELEMEDGE = 0.75;\n";
    input << "Line1: Line = (D6);\n";
    input.close();

    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    DriftRep drift("D6");
    drift.setElementLength(0.3);
    drift.setPlacementPose(PlacementPose(CoordinateSystemTrafo(
            Vector3(1.0, 2.0, 3.0),
            Quaternion(std::cos(Physics::pi / 8.0), 0.0, 0.0, std::sin(Physics::pi / 8.0)))));
    drift.fixPosition();

    OpalBeamline beamline;
    beamline.visit(drift, visitor, *bunch);
    beamline.save3DInput();

    std::ifstream output(outputPath("_3D.opal"));
    ASSERT_TRUE(output.is_open());
    const std::string rewritten(
            (std::istreambuf_iterator<char>(output)), std::istreambuf_iterator<char>());

    EXPECT_EQ(rewritten.find("ELEMEDGE"), std::string::npos);
    EXPECT_NE(rewritten.find("X = 1"), std::string::npos);
    EXPECT_NE(rewritten.find("Y = 2"), std::string::npos);
    EXPECT_NE(rewritten.find("Z = 3"), std::string::npos);
    EXPECT_NE(rewritten.find("THETA = 0.0 * PI / 180"), std::string::npos);
    EXPECT_NE(rewritten.find("PHI = 0.0 * PI / 180"), std::string::npos);
    EXPECT_NE(rewritten.find("PSI = 315.0 * PI / 180"), std::string::npos);
}

TEST_F(OpalBeamlinePlacementTest, PrintPlacementSummaryReportsBodyPoseOnOneLinePerElement) {
    DriftRep drift("D7");
    drift.setElementLength(0.3);
    drift.setPlacementPose(PlacementPose(CoordinateSystemTrafo(
            Vector3(1.0, 2.0, 3.0),
            Quaternion(std::cos(Physics::pi / 8.0), 0.0, 0.0, std::sin(Physics::pi / 8.0)))));
    drift.fixPosition();

    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    SBendRep bend("B1");
    bend.getGeometry() = PlanarArcGeometry(0.5, Physics::pi / 8.0);
    bend.setElementLength(0.5);
    bend.setBendAngle(Physics::pi / 8.0);
    bend.setFringeHalfGap(0.02);
    bend.setFringeIntegral(0.5);
    bend.setElementPosition(10.0);

    OpalBeamline beamline;
    beamline.visit(drift, visitor, *bunch);
    beamline.visit(bend, visitor, *bunch);
    beamline.prepareSections();

    std::ostringstream summary;
    beamline.printPlacementSummary(summary);
    const std::string text = summary.str();

    EXPECT_NE(text.find("Element"), std::string::npos);
    EXPECT_NE(text.find("X\tY\tZ\tTHETA\tPHI\tPSI"), std::string::npos);
    EXPECT_NE(text.find("------------------------"), std::string::npos);
    EXPECT_NE(text.find("D7"), std::string::npos);
    EXPECT_NE(text.find("1.0\t2.0\t3.0\t0.0\t0.0\t315.0"), std::string::npos);
    EXPECT_NE(text.find("B1"), std::string::npos);
    EXPECT_NE(text.find("Local body/field extents:"), std::string::npos);
    EXPECT_NE(text.find("BODY_BEGIN\tBODY_END\tFIELD_BEGIN\tFIELD_END"), std::string::npos);
    EXPECT_NE(text.find("B1"), std::string::npos);
    EXPECT_NE(text.find("0.0\t0.5\t-0.01\t0.51"), std::string::npos);
    EXPECT_LT(text.find("D7"), text.find("B1"));
}

TEST_F(OpalBeamlinePlacementTest, Save3DLatticeExportsFieldExtentMarkersWhenSupportDiffers) {
    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    SBendRep bend("B2");
    bend.getGeometry() = PlanarArcGeometry(0.5, Physics::pi / 8.0);
    bend.setElementLength(0.5);
    bend.setBendAngle(Physics::pi / 8.0);
    bend.setFringeHalfGap(0.02);
    bend.setFringeIntegral(0.5);
    bend.setElementPosition(0.0);

    OpalBeamline beamline;
    beamline.visit(bend, visitor, *bunch);

    if (ippl::Comm->rank() != 0) {
        return;
    }

    beamline.compute3DLattice();
    beamline.save3DLattice();

    std::ifstream txt(outputPath("_ElementPositions.txt"));
    ASSERT_TRUE(txt.is_open());
    const std::string content(
            (std::istreambuf_iterator<char>(txt)), std::istreambuf_iterator<char>());

    EXPECT_NE(content.find("\"FIELD BEGIN: B2\""), std::string::npos);
    EXPECT_NE(content.find("\"FIELD END: B2\""), std::string::npos);
}

TEST_F(OpalBeamlinePlacementTest, PrepareSectionsPlacesCompatibilitySBendByEntryEdge) {
    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    SBendRep bend("B3");
    bend.getGeometry() = PlanarArcGeometry(1.0, Physics::pi / 4.0);
    bend.setElementLength(1.0);
    bend.setBendAngle(Physics::pi / 4.0);
    bend.setElementPosition(0.0);

    OpalBeamline beamline;
    beamline.visit(bend, visitor, *bunch);
    beamline.prepareSections();

    const auto elements = beamline.getElements();
    ASSERT_EQ(elements.size(), 1u);
    const auto component = *elements.begin();

    expectVectorNear(
            beamline.getNominalEntryTransform(component).getOrigin(), Vector3(0.0, 0.0, 0.0));

    const Vector3 bodyOrigin =
            beamline.getPlacedElement(component).getNominalBodyTransform().getOrigin();
    EXPECT_NEAR(bodyOrigin(0), -0.09691958937035704, 1.0e-12);
    EXPECT_NEAR(bodyOrigin(1), 0.0, 1.0e-12);
    EXPECT_NEAR(bodyOrigin(2), 0.48724767920221634, 1.0e-12);
}

TEST_F(OpalBeamlinePlacementTest, ReportPortContinuityDiagnosticsFlagsGapToFollowingElement) {
    DriftRep drift("D7");
    drift.setElementLength(0.3);
    drift.setPlacementPose(
            PlacementPose(CoordinateSystemTrafo(Vector3(0.0, 0.0, 0.0), Quaternion())));
    drift.fixPosition();

    DriftRep drift2("D8");
    drift2.setElementLength(0.4);
    drift2.setPlacementPose(
            PlacementPose(CoordinateSystemTrafo(Vector3(0.0, 0.0, 0.5), Quaternion())));
    drift2.fixPosition();

    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    OpalBeamline beamline;
    beamline.visit(drift, visitor, *bunch);
    beamline.visit(drift2, visitor, *bunch);
    beamline.prepareSections();

    std::ostringstream diagnostics;
    EXPECT_TRUE(beamline.reportPortContinuityDiagnostics(diagnostics));

    const std::string text = diagnostics.str();
    EXPECT_NE(text.find("Adjacent element port discontinuities:"), std::string::npos);
    EXPECT_NE(text.find("D7 -> D8"), std::string::npos);
    EXPECT_NE(text.find("|gap|=0.2"), std::string::npos);
    EXPECT_NE(text.find("dTHETA=0.0 deg, dPHI=0.0 deg, dPSI=0.0 deg"), std::string::npos);
    EXPECT_NE(text.find("adjust D8 body pose by dX=-0.0, dY=-0.0, dZ=-0.2"), std::string::npos);
}

TEST_F(OpalBeamlinePlacementTest,
       ReportPortContinuityDiagnosticsNormalizesPsiMismatchToSignedAngle) {
    DriftRep drift("D7");
    drift.setElementLength(0.3);
    drift.setPlacementPose(
            PlacementPose(CoordinateSystemTrafo(Vector3(0.0, 0.0, 0.0), Quaternion())));
    drift.fixPosition();

    DriftRep drift2("D8");
    drift2.setElementLength(0.4);
    drift2.setPlacementPose(PlacementPose(CoordinateSystemTrafo(
            Vector3(0.0, 0.0, 0.3),
            Quaternion(std::cos(Physics::pi / 8.0), 0.0, 0.0, std::sin(Physics::pi / 8.0)))));
    drift2.fixPosition();

    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    OpalBeamline beamline;
    beamline.visit(drift, visitor, *bunch);
    beamline.visit(drift2, visitor, *bunch);
    beamline.prepareSections();

    std::ostringstream diagnostics;
    EXPECT_TRUE(beamline.reportPortContinuityDiagnostics(diagnostics));

    const std::string text = diagnostics.str();
    EXPECT_NE(text.find("dTHETA=0.0 deg, dPHI=0.0 deg, dPSI="), std::string::npos);
    EXPECT_EQ(text.find("315.0 deg"), std::string::npos);
    EXPECT_NE(text.find("adjust D8 body pose by"), std::string::npos);
}

TEST_F(OpalBeamlinePlacementTest, BendMeshesAsPinkCurvedTubeInElementPositionsExport) {
    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    SBendRep bend("B1");
    bend.getGeometry() = PlanarArcGeometry(1.0, Physics::pi / 4.0);
    bend.setElementLength(1.0);
    bend.setBendAngle(Physics::pi / 4.0);
    bend.setAperture(ApertureType::ELLIPTICAL, std::vector<double>{0.03, 0.04, 1.0});
    bend.setElementPosition(0.0);

    OpalBeamline beamline;
    beamline.visit(bend, visitor, *bunch);

    if (ippl::Comm->rank() != 0) {
        return;
    }

    beamline.compute3DLattice();
    beamline.save3DLattice();

    std::ifstream py(outputPath("_ElementPositions.py"));
    ASSERT_TRUE(py.is_open());
    const std::string script(
            (std::istreambuf_iterator<char>(py)), std::istreambuf_iterator<char>());
    EXPECT_NE(script.find("color = [1]"), std::string::npos);
    EXPECT_NE(script.find("numVertices = [2424]"), std::string::npos);
    EXPECT_NE(script.find("('Dipole', (1.0, 0.4, 0.7))"), std::string::npos);
    EXPECT_NE(script.find("lookup_table.append([1.0, 0.4, 0.7, 1.0])"), std::string::npos);
}

TEST_F(OpalBeamlinePlacementTest, RuntimeElementQueryRecognizesBendOnExportedCenterline) {
    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    SBendRep bend("B1");
    bend.getGeometry() = PlanarArcGeometry(1.0, Physics::pi / 4.0);
    bend.setElementLength(1.0);
    bend.setBendAngle(Physics::pi / 4.0);
    bend.setElementPosition(0.0);

    OpalBeamline beamline;
    beamline.visit(bend, visitor, *bunch);
    beamline.prepareSections();

    const auto allElements = beamline.getElements();
    ASSERT_EQ(allElements.size(), 1u);
    const auto component    = *allElements.begin();
    const auto* bendElement = dynamic_cast<const BendBase*>(component.get());
    ASSERT_NE(bendElement, nullptr);

    const std::vector<Vector3> designPath = bendElement->getDesignPath(9);
    ASSERT_GE(designPath.size(), 9u);
    const Vector3 labPoint =
            beamline.getPlacedElement(component).getNominalBodyTransform().transformFrom(
                    designPath[4]);

    const auto active = beamline.getElements(labPoint);
    EXPECT_EQ(active.size(), 1u);
    EXPECT_EQ((*active.begin())->getName(), "B1");
}

TEST_F(OpalBeamlinePlacementTest, SBendFieldChartMapsDesignArcToArcLengthCoordinate) {
    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    SBendRep bend("B4");
    bend.getGeometry() = PlanarArcGeometry(1.0, Physics::pi / 4.0);
    bend.setElementLength(1.0);
    bend.setBendAngle(Physics::pi / 4.0);
    bend.setElementPosition(0.0);

    OpalBeamline beamline;
    beamline.visit(bend, visitor, *bunch);
    beamline.prepareSections();

    const auto elements = beamline.getElements();
    ASSERT_EQ(elements.size(), 1u);
    const auto component    = *elements.begin();
    const auto* bendElement = dynamic_cast<const BendBase*>(component.get());
    ASSERT_NE(bendElement, nullptr);

    const std::vector<Vector3> designPath = bendElement->getDesignPath(9);
    ASSERT_GE(designPath.size(), 9u);
    const std::size_t midIndex = designPath.size() / 2u;

    const Vector3 entryLab =
            beamline.getPlacedElement(component).getNominalBodyTransform().transformFrom(
                    designPath.front());
    const Vector3 midLab =
            beamline.getPlacedElement(component).getNominalBodyTransform().transformFrom(
                    designPath[midIndex]);
    const Vector3 exitLab =
            beamline.getPlacedElement(component).getNominalBodyTransform().transformFrom(
                    designPath.back());

    const Vector3 entryLocal = beamline.transformToFieldLocalCS(component, entryLab);
    const Vector3 midLocal   = beamline.transformToFieldLocalCS(component, midLab);
    const Vector3 exitLocal  = beamline.transformToFieldLocalCS(component, exitLab);

    EXPECT_NEAR(entryLocal(2), 0.0, 1.0e-12);
    EXPECT_NEAR(midLocal(2), 0.5, 1.0e-12);
    EXPECT_NEAR(exitLocal(2), 1.0, 1.0e-12);
}

TEST_F(OpalBeamlinePlacementTest, RBendFieldChartMapsChordCoordinatesToBodyLength) {
    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    RBendRep bend("RB4");
    bend.getGeometry().setElementLength(1.0);
    bend.getGeometry().setBendAngle(Physics::pi / 4.0);
    bend.setElementLength(1.0);
    bend.setBendAngle(Physics::pi / 4.0);
    bend.setElementPosition(0.0);

    OpalBeamline beamline;
    beamline.visit(bend, visitor, *bunch);
    beamline.prepareSections();

    const auto elements = beamline.getElements();
    ASSERT_EQ(elements.size(), 1u);
    const auto component = *elements.begin();

    const Vector3 entryLab = beamline.getNominalEntryTransform(component).getOrigin();
    const Vector3 bodyLab =
            beamline.getPlacedElement(component).getNominalBodyTransform().getOrigin();
    const Vector3 exitLab = beamline.getNominalExitTransform(component).getOrigin();

    const Vector3 entryLocal = beamline.transformToFieldLocalCS(component, entryLab);
    const Vector3 bodyLocal  = beamline.transformToFieldLocalCS(component, bodyLab);
    const Vector3 exitLocal  = beamline.transformToFieldLocalCS(component, exitLab);

    const double referenceLength = bend.getGeometry().getArcLength();
    EXPECT_NEAR(entryLocal(2), 0.0, 1.0e-12);
    EXPECT_NEAR(bodyLocal(2), 0.5 * referenceLength, 1.0e-12);
    EXPECT_NEAR(exitLocal(2), referenceLength, 1.0e-12);
}

TEST_F(OpalBeamlinePlacementTest, ElementPositionsScriptStaysValidWithEmptyTriangleMeshBlocks) {
    MeshGenerator generator;

    SBendRep bend("B1");
    bend.getGeometry() = PlanarArcGeometry(1.0, Physics::pi / 4.0);
    bend.setElementLength(1.0);
    bend.setBendAngle(Physics::pi / 4.0);
    bend.setAperture(ApertureType::ELLIPTICAL, std::vector<double>{0.03, 0.04, 1.0});
    bend.setPlacementPose(
            PlacementPose(CoordinateSystemTrafo(Vector3(0.0, 0.0, 0.0), Quaternion())));

    EmptyRBend3D empty("B3D");

    generator.add(bend);
    generator.add(empty);

    if (ippl::Comm->rank() != 0) {
        return;
    }

    generator.write("TestOpalBeamlinePlacementSyntax");

    const std::filesystem::path scriptPath = outputPath("Syntax_ElementPositions.py");
    ASSERT_TRUE(std::filesystem::exists(scriptPath));

    const std::filesystem::path pycachePath = outputPath("_pycache");
    std::filesystem::create_directories(pycachePath);
    const std::string command = "PYTHONPYCACHEPREFIX=\"" + pycachePath.string()
                                + "\" python3 -m py_compile \"" + scriptPath.string() + "\"";
    EXPECT_EQ(std::system(command.c_str()), 0);
}
