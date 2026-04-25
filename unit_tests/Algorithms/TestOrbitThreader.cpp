#include "gtest/gtest.h"

#include "AbstractObjects/OpalData.h"
#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/OrbitThreader.h"
#include "Algorithms/PartData.h"
#include "BeamlineCore/MultipoleRep.h"
#include "BeamlineGeometry/NullGeometry.h"
#include "BeamlineGeometry/PlacementPose.h"
#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"
#include "Fields/NullField.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
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

    /**
     * @brief Mock component with zero body extent but finite field-support extent.
     *
     * This models the post-redesign case where placement/geometry uses the body
     * extent while tracking constraints must use the field-support interval.
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

    std::shared_ptr<FieldSupportOnlyComponent> makePlacedFieldSupportOnlyComponent(
            const std::string& name, const double entryPosition, const double fieldLength) {
        auto component = std::make_shared<FieldSupportOnlyComponent>(name, 0.0, fieldLength);
        component->setElementLength(0.0);
        component->setPlacementPose(PlacementPose(
                CoordinateSystemTrafo(Vector_t<double, 3>(0.0, 0.0, entryPosition), Quaternion())));
        component->fixPosition();
        return component;
    }

    std::shared_ptr<FieldSolverCmd> fsCmdBase_m;
    std::shared_ptr<DataSink> dataSink_m;
};

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

TEST_F(OrbitThreaderTest, UsesFieldSupportExtentForLengthCheck) {
    auto bunch = makeBunch(0);

    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    OpalBeamline beamline;
    auto component = makePlacedFieldSupportOnlyComponent("TW_LIKE", 0.0, 0.1);
    beamline.visit(*component, visitor, *bunch);
    beamline.prepareSections();

    StepSizeConfig stepSizes;
    stepSizes.push_back(1.0e-12, 0.2, 64);
    stepSizes.resetIterator();

    PartData reference(1.0, 9.382720813e8, 1.0e6);
    OrbitThreader threader(
            reference, Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0, 0.0, 1.0), 0.0, 0.0, 0.0,
            1.0e-12, stepSizes, beamline);

    EXPECT_NO_THROW(threader.execute());
}
