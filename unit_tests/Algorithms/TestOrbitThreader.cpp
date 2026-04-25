#include "gtest/gtest.h"

#include "AbstractObjects/OpalData.h"
#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/OrbitThreader.h"
#include "Algorithms/PartData.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/MultipoleRep.h"
#include "BeamlineCore/SBendRep.h"
#include "BeamlineGeometry/NullGeometry.h"
#include "BeamlineGeometry/PlacementPose.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

#include <filesystem>
#include <set>

extern Inform* gmsg;

namespace {
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

class OrbitThreaderTest : public ::testing::Test {
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
        OpalData::getInstance()->storeInputFn("TestOrbitThreader.opal");
        OpalData::getInstance()->setOpenMode(OpalData::OpenMode::WRITE);
        std::filesystem::create_directories(OpalData::getInstance()->getAuxiliaryOutputDirectory());
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

    std::shared_ptr<MultipoleRep> makePlacedQuadrupole(
            const std::string& name, const double length, const double entryPosition,
            const double normalComponent) {
        auto quadrupole = std::make_shared<MultipoleRep>(name);
        quadrupole->setElementLength(length);
        quadrupole->setNormalComponent(1, normalComponent);
        quadrupole->setPlacementPose(PlacementPose(
                CoordinateSystemTrafo(Vector_t<double, 3>(0.0, 0.0, entryPosition), Quaternion())));
        quadrupole->fixPosition();
        return quadrupole;
    }

    std::shared_ptr<FieldSolverCmd> fsCmdBase_m;
    std::shared_ptr<DataSink> dataSink_m;
};

namespace {
    class TestableOrbitThreader : public OrbitThreader {
    public:
        using OrbitThreader::OrbitThreader;
        using OrbitThreader::validateVisitedElements;
    };
}  // namespace

TEST_F(OrbitThreaderTest, ExposesEmptyReferencePathModelBeforeExecution) {
    StepSizeConfig stepSizes;
    stepSizes.push_back(1.0e-12, 1.0, 8);
    stepSizes.resetIterator();

    PartData reference(1.0, 9.382720813e8, 1.0e6);
    OpalBeamline beamline;
    OrbitThreader threader(
            reference, Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0, 0.0, 1.0), 0.0, 0.0, 0.0,
            1.0e-12, stepSizes, beamline);

    EXPECT_TRUE(threader.getReferencePathModel().empty());
    EXPECT_TRUE(threader.getActionRangeRegistrationModel().empty());
}

TEST_F(OrbitThreaderTest, ExecutesOverlapAndBuildsTracedAndRegistrationModels) {
    auto bunch = makeBunch(0);

    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    OpalBeamline beamline;
    auto longQuadrupole  = makePlacedQuadrupole("Q_LONG", 0.5, 0.0, 0.5);
    auto shortQuadrupole = makePlacedQuadrupole("Q_SHORT", 0.1, 0.45, 0.8);
    beamline.visit(*longQuadrupole, visitor, *bunch);
    beamline.visit(*shortQuadrupole, visitor, *bunch);
    beamline.prepareSections();

    StepSizeConfig stepSizes;
    stepSizes.push_back(1.0e-11, 0.7, 512);
    stepSizes.resetIterator();

    PartData reference(1.0, 9.382720813e8, 1.0e6);
    OrbitThreader threader(
            reference, Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0, 0.0, 1.0), 0.0, 0.0, 0.0,
            1.0e-11, stepSizes, beamline);

    threader.execute();

    const ReferencePathModel& tracedModel = threader.getReferencePathModel();
    ASSERT_FALSE(tracedModel.empty());

    bool foundOverlap = false;
    for (const auto& segment : tracedModel.getSegments()) {
        if (segment.getActiveElements().size() != 2u) {
            continue;
        }

        std::set<std::string> activeNames;
        for (const auto& element : segment.getActiveElements()) {
            activeNames.insert(element->getName());
        }

        if (activeNames == std::set<std::string>{"Q_LONG", "Q_SHORT"}) {
            foundOverlap = true;
            EXPECT_NEAR(segment.getBegin(), 0.45, 0.02);
            EXPECT_NEAR(segment.getEnd(), 0.5, 0.02);
            break;
        }
    }
    EXPECT_TRUE(foundOverlap);

    const ReferencePathModel& registrationModel = threader.getActionRangeRegistrationModel();
    ASSERT_EQ(registrationModel.size(), 2u);

    std::set<std::string> registeredNames;
    for (const auto& segment : registrationModel.getSegments()) {
        ASSERT_EQ(segment.getActiveElements().size(), 1u);
        EXPECT_TRUE(segment.hasLegacyElementEdge());
        registeredNames.insert((*segment.getActiveElements().begin())->getName());
    }
    EXPECT_EQ(registeredNames, (std::set<std::string>{"Q_LONG", "Q_SHORT"}));
}

TEST_F(OrbitThreaderTest, ThrowsDiagnosticForPlacedElementOutsideTracedReferencePath) {
    auto bunch = makeBunch(0);

    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    OpalBeamline beamline;
    auto drift = std::make_shared<DriftRep>("D0");
    drift->setElementLength(3.292883);
    drift->setPlacementPose(
            PlacementPose(CoordinateSystemTrafo(Vector_t<double, 3>(0.0, 0.0, 0.0), Quaternion())));
    drift->fixPosition();
    beamline.visit(*drift, visitor, *bunch);

    auto bend           = std::make_shared<SBendRep>("B1");
    bend->getGeometry() = PlanarArcGeometry(0.5, Physics::pi / 4.0);
    bend->setElementLength(0.5);
    bend->setBendAngle(Physics::pi / 8.0);
    bend->setPlacementPose(
            PlacementPose(CoordinateSystemTrafo(Vector_t<double, 3>(0.0, 0.0, 3.3), Quaternion())));
    bend->fixPosition();
    beamline.visit(*bend, visitor, *bunch);
    beamline.prepareSections();

    StepSizeConfig stepSizes;
    stepSizes.push_back(1.0e-11, 5.0, 4096);
    stepSizes.resetIterator();

    PartData reference(1.0, 9.382720813e8, 1.0e6);
    TestableOrbitThreader threader(
            reference, Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0, 0.0, 1.0), 0.0, 0.0, 0.0,
            1.0e-11, stepSizes, beamline);

    const FieldList allElements = beamline.getElementByType(ElementType::ANY);
    ASSERT_EQ(allElements.size(), 2u);

    EXPECT_THROW(
            {
                try {
                    threader.validateVisitedElements(allElements, std::set<std::string>{"D0"}, 0.0);
                } catch (const OpalException& ex) {
                    const std::string what = ex.what();
                    EXPECT_NE(what.find("B1"), std::string::npos);
                    EXPECT_NE(
                            what.find("never reached by the traced reference path"),
                            std::string::npos);
                    EXPECT_NE(what.find("nominal entry"), std::string::npos);
                    throw;
                }
            },
            OpalException);
}
