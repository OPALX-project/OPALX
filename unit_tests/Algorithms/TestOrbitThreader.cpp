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
#include "Fields/NullField.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>

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

    struct DesignPathRow {
        double s;
        double rx;
        double ry;
        double rz;
        double px;
        double py;
        double pz;
        double efx;
        double efy;
        double efz;
        double bfx;
        double bfy;
        double bfz;
        double ekin;
        double time;
    };

    DesignPathRow readLastDesignPathRow(const std::filesystem::path& path) {
        std::ifstream stream(path);
        EXPECT_TRUE(stream.is_open());

        DesignPathRow lastRow{};
        std::string line;
        bool foundRow = false;
        while (std::getline(stream, line)) {
            if (line.empty() || line.front() == '#') {
                continue;
            }

            std::istringstream row(line);
            row >> lastRow.s >> lastRow.rx >> lastRow.ry >> lastRow.rz >> lastRow.px >> lastRow.py
                    >> lastRow.pz >> lastRow.efx >> lastRow.efy >> lastRow.efz >> lastRow.bfx
                    >> lastRow.bfy >> lastRow.bfz >> lastRow.ekin >> lastRow.time;
            EXPECT_FALSE(row.fail());
            foundRow = true;
        }

        EXPECT_TRUE(foundRow);
        return lastRow;
    }

    DesignPathRow readFirstDesignPathRow(const std::filesystem::path& path) {
        std::ifstream stream(path);
        EXPECT_TRUE(stream.is_open());

        DesignPathRow firstRow{};
        std::string line;
        while (std::getline(stream, line)) {
            if (line.empty() || line.front() == '#') {
                continue;
            }

            std::istringstream row(line);
            row >> firstRow.s >> firstRow.rx >> firstRow.ry >> firstRow.rz >> firstRow.px
                    >> firstRow.py >> firstRow.pz >> firstRow.efx >> firstRow.efy >> firstRow.efz
                    >> firstRow.bfx >> firstRow.bfy >> firstRow.bfz >> firstRow.ekin
                    >> firstRow.time;
            EXPECT_FALSE(row.fail());
            return firstRow;
        }

        ADD_FAILURE() << "No data rows found in design-path file: " << path;
        return firstRow;
    }

    std::size_t countOccurrences(const std::string& text, const std::string& needle) {
        std::size_t count    = 0;
        std::size_t position = 0;
        while ((position = text.find(needle, position)) != std::string::npos) {
            ++count;
            position += needle.size();
        }
        return count;
    }

    class ScopedInformRedirect {
    public:
        explicit ScopedInformRedirect(std::ostream& stream) : previous_(gmsg) {
            gmsg = new Inform("OPAL-X", stream);
            gmsg->setOutputLevel(5);
        }

        ~ScopedInformRedirect() {
            delete gmsg;
            gmsg = previous_;
        }

    private:
        Inform* previous_;
    };

    /**
     * @brief Mock component with zero body extent but finite field-support extent.
     *
     * This models the case where placement uses the body extent while tracking
     * constraints use the field-support interval.
     */
    class FieldSupportOnlyComponent final : public Component {
    public:
        FieldSupportOnlyComponent(
                const std::string& name, const double fieldBegin, const double fieldEnd)
            : Component(name), fieldBegin_m(fieldBegin), fieldEnd_m(fieldEnd) {}

        void accept(BeamlineVisitor&) const override {}
        ElementBase* clone() const override { return new FieldSupportOnlyComponent(*this); }

        EMField& getField() override { return field_m; }
        const EMField& getField() const override { return field_m; }

        bool apply(const std::shared_ptr<ParticleContainer_t>&) override { return false; }

        bool apply(
                const size_t&, const double&, Vector_t<double, 3>&, Vector_t<double, 3>&) override {
            return false;
        }

        bool apply(
                const Vector_t<double, 3>&, const Vector_t<double, 3>&, const double&,
                Vector_t<double, 3>&, Vector_t<double, 3>&) override {
            return false;
        }

        bool applyToReferenceParticle(
                const Vector_t<double, 3>&, const Vector_t<double, 3>&, const double&,
                Vector_t<double, 3>&, Vector_t<double, 3>&) override {
            return false;
        }

        void initialise(PartBunch_t*, double&, double&) override {}
        void finalise() override {}
        bool bends() const override { return false; }

        void getFieldExtend(double& zBegin, double& zEnd) const override {
            zBegin = fieldBegin_m;
            zEnd   = fieldEnd_m;
        }

        ElementType getType() const override { return ElementType::ANY; }

        BGeometryBase& getGeometry() override { return geometry_m; }
        const BGeometryBase& getGeometry() const override { return geometry_m; }

    private:
        double fieldBegin_m;
        double fieldEnd_m;
        NullGeometry geometry_m;
        NullField field_m;
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
        bunch->getParticleContainer()->createParticles(numParticles);
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

    std::shared_ptr<FieldSupportOnlyComponent> makePlacedFieldSupportOnlyComponent(
            const std::string& name, const double entryPosition, const double fieldLength) {
        auto component = std::make_shared<FieldSupportOnlyComponent>(name, 0.0, fieldLength);
        component->setElementLength(0.0);
        component->setPlacementPose(PlacementPose(
                CoordinateSystemTrafo(Vector_t<double, 3>(0.0, 0.0, entryPosition), Quaternion())));
        component->fixPosition();
        return component;
    }

    std::shared_ptr<SBendRep> makePlacedLinearEdgeSBend(
            const std::string& name, const double entryPosition, const double length,
            const double bendAngle, const double edgeAngle, const double halfGap,
            const double fringeIntegral) {
        auto bend           = std::make_shared<SBendRep>(name);
        bend->getGeometry() = PlanarArcGeometry(length, bendAngle);
        bend->setElementLength(length);
        bend->setBendAngle(bendAngle);
        bend->setEntranceAngle(edgeAngle);
        bend->setExitAngle(edgeAngle);
        bend->setFringeHalfGap(halfGap);
        bend->setFringeIntegral(fringeIntegral);
        bend->setElementPosition(entryPosition);
        bend->setPlacementPose(PlacementPose(
                CoordinateSystemTrafo(Vector_t<double, 3>(0.0, 0.0, entryPosition), Quaternion())));
        bend->fixPosition();
        return bend;
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

TEST_F(OrbitThreaderTest, PrintsOnlyFinalThreadingSummary) {
    auto bunch = makeBunch(0);

    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    OpalBeamline beamline;
    auto firstDrift = std::make_shared<DriftRep>("D0");
    firstDrift->setElementLength(0.2);
    firstDrift->setPlacementPose(
            PlacementPose(CoordinateSystemTrafo(Vector_t<double, 3>(0.0, 0.0, 0.0), Quaternion())));
    firstDrift->fixPosition();
    beamline.visit(*firstDrift, visitor, *bunch);

    auto secondDrift = std::make_shared<DriftRep>("D1");
    secondDrift->setElementLength(0.2);
    secondDrift->setPlacementPose(
            PlacementPose(CoordinateSystemTrafo(Vector_t<double, 3>(0.0, 0.0, 0.4), Quaternion())));
    secondDrift->fixPosition();
    beamline.visit(*secondDrift, visitor, *bunch);
    beamline.prepareSections();

    StepSizeConfig stepSizes;
    stepSizes.push_back(1.0e-11, 0.8, 512);
    stepSizes.resetIterator();

    PartData reference(1.0, 9.382720813e8, 1.0e6);
    std::stringstream output;
    {
        ScopedInformRedirect redirect(output);
        OrbitThreader threader(
                reference, Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0, 0.0, 1.0), 0.0, 0.0,
                0.0, 1.0e-11, stepSizes, beamline);
        ASSERT_NO_THROW(threader.execute());
    }

    if (ippl::Comm->rank() != 0) {
        return;
    }

    const std::string text = output.str();
    EXPECT_EQ(countOccurrences(text, "OrbitThreader maxDistance="), 1u);
    EXPECT_EQ(countOccurrences(text, "OrbitThreader #visited elements"), 1u);
    EXPECT_EQ(countOccurrences(text, "OrbitThreader #elements"), 0u);
    EXPECT_EQ(text.find("1.79769e+308"), std::string::npos);
    EXPECT_NE(text.find("OrbitThreader #visited elements = 2"), std::string::npos);
}

TEST_F(OrbitThreaderTest, UsesFieldSupportExtentForLengthCheck) {
    auto bunch = makeBunch(0);

    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    OpalBeamline beamline;
    auto component = makePlacedFieldSupportOnlyComponent("FSUP", 0.0, 0.4);
    beamline.visit(*component, visitor, *bunch);
    beamline.prepareSections();

    StepSizeConfig stepSizes;
    stepSizes.push_back(1.0e-9, 0.5, 8);
    stepSizes.resetIterator();

    PartData reference(1.0, 9.382720813e8, 1.0e6);
    OrbitThreader threader(
            reference, Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0, 0.0, 1.0), 0.0, 0.0, 0.0,
            1.0e-9, stepSizes, beamline);

    EXPECT_THROW(threader.execute(), OpalException);
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
                    threader.validateVisitedElements(
                            allElements, std::set<std::string>{"D0"}, 0.0, 5.0);
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

TEST_F(OrbitThreaderTest, DoesNotReportUnvisitedElementsBeyondTracedStepLimit) {
    auto bunch = makeBunch(0);

    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    OpalBeamline beamline;
    auto drift = std::make_shared<DriftRep>("D0");
    drift->setElementLength(3.0);
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
    stepSizes.push_back(1.0e-11, 5.0, 1);
    stepSizes.resetIterator();

    PartData reference(1.0, 9.382720813e8, 1.0e6);
    TestableOrbitThreader threader(
            reference, Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0, 0.0, 1.0), 0.0, 0.0, 0.0,
            1.0e-11, stepSizes, beamline);

    const FieldList allElements = beamline.getElementByType(ElementType::ANY);
    ASSERT_NO_THROW(
            threader.validateVisitedElements(allElements, std::set<std::string>{"D0"}, 0.0, 0.01));
}

TEST_F(OrbitThreaderTest, RegistersSBendActionRangeFromFieldSupportExtent) {
    auto bunch = makeBunch(0);

    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    OpalBeamline beamline;
    auto bend = makePlacedLinearEdgeSBend(
            "B1", 0.0, 1.0, Physics::pi / 4.0, Physics::pi / 8.0, 0.015, 0.50);
    beamline.visit(*bend, visitor, *bunch);
    beamline.prepareSections();

    StepSizeConfig stepSizes;
    stepSizes.push_back(5.0e-12, 1.2, 4096);
    stepSizes.resetIterator();

    PartData reference(1.0, 9.382720813e8, 5.90e8);
    OrbitThreader threader(
            reference, Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0, 0.0, 1.0), 0.0, 0.0, 0.0,
            5.0e-12, stepSizes, beamline);

    ASSERT_NO_THROW(threader.execute());

    const ReferencePathModel& registrationModel = threader.getActionRangeRegistrationModel();
    ASSERT_EQ(registrationModel.size(), 1u);

    const auto& segment              = registrationModel.getSegments().front();
    const double expectedEntryFringe = 10.0 * 0.015 / std::abs(std::cos(Physics::pi / 8.0));
    EXPECT_TRUE(segment.hasLegacyElementEdge());
    ASSERT_TRUE(segment.getLegacyElementEdge().has_value());
    EXPECT_NEAR(segment.getLegacyElementEdge().value(), 0.0, 1.0e-12);
    EXPECT_NEAR(segment.getBegin(), -expectedEntryFringe, 1.0e-6);
    EXPECT_NEAR(segment.getEnd(), 1.0 + expectedEntryFringe, 1.0e-6);
}

TEST_F(OrbitThreaderTest, RegistersExplicitPositionedSBendActionRangeWithoutElementEdge) {
    auto bunch = makeBunch(0);

    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    constexpr double length         = 1.0;
    constexpr double bendAngle      = Physics::pi / 4.0;
    constexpr double edgeAngle      = Physics::pi / 8.0;
    constexpr double halfGap        = 0.015;
    constexpr double fringeIntegral = 0.50;

    auto bend           = std::make_shared<SBendRep>("B1");
    bend->getGeometry() = PlanarArcGeometry(length, bendAngle);
    bend->setElementLength(length);
    bend->setBendAngle(bendAngle);
    bend->setEntranceAngle(edgeAngle);
    bend->setExitAngle(edgeAngle);
    bend->setFringeHalfGap(halfGap);
    bend->setFringeIntegral(fringeIntegral);
    bend->setPlacementPose(
            PlacementPose(CoordinateSystemTrafo(Vector_t<double, 3>(0.0, 0.0, 0.0), Quaternion())));
    bend->fixPosition();

    OpalBeamline beamline;
    beamline.visit(*bend, visitor, *bunch);
    beamline.prepareSections();

    StepSizeConfig stepSizes;
    stepSizes.push_back(5.0e-12, 1.2, 4096);
    stepSizes.resetIterator();

    PartData reference(1.0, 9.382720813e8, 5.90e8);
    OrbitThreader threader(
            reference, Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0, 0.0, 1.0), 0.0, 0.0, 0.0,
            5.0e-12, stepSizes, beamline);

    ASSERT_NO_THROW(threader.execute());

    const ReferencePathModel& registrationModel = threader.getActionRangeRegistrationModel();
    ASSERT_EQ(registrationModel.size(), 1u);

    const auto& segment              = registrationModel.getSegments().front();
    const double expectedEntryFringe = 10.0 * halfGap / std::abs(std::cos(edgeAngle));
    ASSERT_EQ(segment.getActiveElements().size(), 1u);
    EXPECT_EQ((*segment.getActiveElements().begin())->getName(), "B1");
    EXPECT_NEAR(segment.getEnd() - segment.getBegin(), length + 2.0 * expectedEntryFringe, 1.0e-6);
}

TEST_F(OrbitThreaderTest, LogsPostStepDesignPathStateAtZStopForSimpleDrift) {
    constexpr double dt     = 1.0e-11;
    const double stepLength = Physics::c * dt / std::sqrt(2.0);
    const double zStop      = 4.0 * stepLength;

    StepSizeConfig stepSizes;
    stepSizes.push_back(dt, zStop, 64);
    stepSizes.resetIterator();

    PartData reference(1.0, 9.382720813e8, 1.0e6);
    OpalBeamline beamline;

    {
        OrbitThreader threader(
                reference, Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0, 0.0, 1.0), 0.0, 0.0,
                0.0, dt, stepSizes, beamline);
        threader.execute();
    }

    if (ippl::Comm->rank() != 0) {
        return;
    }

    const std::filesystem::path designPathFile =
            std::filesystem::path(OpalData::getInstance()->getAuxiliaryOutputDirectory())
            / "TestOrbitThreader_DesignPath.dat";
    const DesignPathRow firstRow = readFirstDesignPathRow(designPathFile);
    const DesignPathRow lastRow  = readLastDesignPathRow(designPathFile);

    EXPECT_NEAR(firstRow.s, 0.0, 1.0e-14);
    EXPECT_NEAR(firstRow.rz, 0.0, 1.0e-14);
    EXPECT_NEAR(firstRow.time, 0.0, 1.0e-14);
    EXPECT_NEAR(lastRow.s, zStop, 1.0e-10);
    EXPECT_NEAR(lastRow.rz, zStop, 1.0e-10);
    EXPECT_NEAR(lastRow.time, 4.0 * dt * Units::s2ns, 1.0e-10);
}

TEST_F(OrbitThreaderTest, LogsInterpolatedTerminalPointWhenSimpleDriftOvershootsZStop) {
    constexpr double dt     = 1.0e-11;
    const double stepLength = Physics::c * dt / std::sqrt(2.0);
    const double zStop      = 3.25 * stepLength;

    StepSizeConfig stepSizes;
    stepSizes.push_back(dt, zStop, 64);
    stepSizes.resetIterator();

    PartData reference(1.0, 9.382720813e8, 1.0e6);
    OpalBeamline beamline;

    {
        OrbitThreader threader(
                reference, Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0, 0.0, 1.0), 0.0, 0.0,
                0.0, dt, stepSizes, beamline);
        threader.execute();
    }

    if (ippl::Comm->rank() != 0) {
        return;
    }

    const std::filesystem::path designPathFile =
            std::filesystem::path(OpalData::getInstance()->getAuxiliaryOutputDirectory())
            / "TestOrbitThreader_DesignPath.dat";
    const DesignPathRow lastRow = readLastDesignPathRow(designPathFile);

    EXPECT_NEAR(lastRow.s, zStop, 1.0e-10);
    EXPECT_NEAR(lastRow.rz, zStop, 1.0e-10);
    EXPECT_NEAR(lastRow.time, 3.25 * dt * Units::s2ns, 1.0e-10);
}
