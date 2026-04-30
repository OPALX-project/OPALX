#include "Algorithms/DefaultVisitor.h"
#include "gtest/gtest.h"

#include "AbstractObjects/OpalData.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineGeometry/NullGeometry.h"
#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
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
