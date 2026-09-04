#include "gtest/gtest.h"

#include "AbstractObjects/OpalData.h"
#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/OrbitThreader.h"
#include "Algorithms/PartData.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/MultipoleRep.h"
#include "BeamlineGeometry/Geometry.h"
#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

#include <filesystem>
#include <set>

extern Inform* gmsg;

namespace {
    std::set<std::string> names(const IndexMap::value_t& elements) {
        std::set<std::string> result;
        for (const auto& element : elements) {
            result.insert(element->getName());
        }
        return result;
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

    /**
     * @brief Mock component with zero body extent but finite field-support extent.
     *
     * This models the post-redesign case where placement/geometry uses the body
     * extent while tracking constraints must use the field-support interval.
     */
    class FieldSupportOnlyComponent final : public ElementBase {
    public:
        FieldSupportOnlyComponent(
                const std::string& name, const double fieldBegin, const double fieldEnd,
                const double longitudinalElectricField = 0.0)
            : ElementBase(name), fieldBegin_m(fieldBegin), fieldEnd_m(fieldEnd),
              electricField_m(longitudinalElectricField) {}

        void accept(BeamlineVisitor&) const override {}
        ElementBase* clone() const override { return new FieldSupportOnlyComponent(*this); }

        void apply(const std::shared_ptr<ParticleContainer_t>&) override {}

        void apply(
                const Vector_t<double, 3>&, const Vector_t<double, 3>&, const double&,
                Vector_t<double, 3>&, Vector_t<double, 3>&) override {}

        bool applyToReferenceParticle(
                const Vector_t<double, 3>& position, const Vector_t<double, 3>&, const double&,
                Vector_t<double, 3>& electric, Vector_t<double, 3>&) override {
            if (position(2) >= fieldBegin_m && position(2) < fieldEnd_m) {
                electric(2) += electricField_m;
            }
            return false;
        }

        void initialise(PartBunch_t*) override {}
        void finalise() override {}

        void getFieldExtent(double& zBegin, double& zEnd) const override {
            zBegin = fieldBegin_m;
            zEnd   = fieldEnd_m;
        }

        ElementType getType() const override { return ElementType::ANY; }

        Geometry& getGeometry() override { return geometry_m; }
        const Geometry& getGeometry() const override { return geometry_m; }

    private:
        double fieldBegin_m;
        double fieldEnd_m;
        double electricField_m;
        Geometry geometry_m{Geometry::makeNull()};
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
        Options::enableLinearTransferMaps = false;
        OpalData::getInstance()->storeInputFn("TestOrbitThreader.opal");
        OpalData::getInstance()->setOpenMode(OpalData::OpenMode::WRITE);
        std::filesystem::create_directories(OpalData::getInstance()->getAuxiliaryOutputDirectory());
    }

    void TearDown() override { Options::enableLinearTransferMaps = false; }

    class TestableFieldSolverCmd : public FieldSolverCmd {
    public:
        void setType(const std::string& t) {
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::TYPE], t);
        }

        void enableParallelDecomposition() {
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTX], true);
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTY], true);
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTZ], true);
            setFieldSolverCmdType();
            setDomainDecomposition();
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
        // Initialize the layout explicitly; the tests also run with two MPI ranks.
        fsCmd->enableParallelDecomposition();
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
        quadrupole->getGeometry().setElementLength(length);
        quadrupole->setNormalComponent(1, normalComponent);
        quadrupole->setCSTrafoGlobal2Local(
                CoordinateSystemTrafo(Vector_t<double, 3>(0.0, 0.0, entryPosition), Quaternion()));
        quadrupole->fixPosition();
        return quadrupole;
    }

    std::shared_ptr<FieldSupportOnlyComponent> makePlacedFieldSupportOnlyComponent(
            const std::string& name, const double entryPosition, const double fieldLength) {
        auto component = std::make_shared<FieldSupportOnlyComponent>(name, 0.0, fieldLength);
        component->getGeometry().setElementLength(0.0);
        component->setCSTrafoGlobal2Local(
                CoordinateSystemTrafo(Vector_t<double, 3>(0.0, 0.0, entryPosition), Quaternion()));
        component->fixPosition();
        return component;
    }

    std::shared_ptr<FieldSolverCmd> fsCmdBase_m;
    std::shared_ptr<DataSink> dataSink_m;
};

TEST_F(OrbitThreaderTest, MapBuilderDoesNotAttachOrMutateSamples) {
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    auto entrance = LinearTransferMapBuilder::initialFrame(
            beamline, Vector_t<double, 3>(0.0, 0.0, 1.0));
    entrance.momentum = Vector_t<double, 3>(0.0, 0.0, 1.0);
    auto exit = entrance;
    exit.position(2) = exit.pathLength = 0.1;
    exit.time = 0.1 * std::sqrt(2.0) / Physics::c;
    const std::vector<LinearTransferMapBuilder::ReferenceSample> samples{{entrance}, {exit}};
    LinearTransferMapBuilder builder(beamline, reference, 1.0e-11);
    const auto result = builder.build(samples, 0.0);
    ASSERT_EQ(result.segments.size(), 1);
    EXPECT_TRUE(result.segments.front().owners.empty());
    ASSERT_TRUE(result.combined);
    EXPECT_NEAR((*result.combined)(0, 1), 0.1, 1.0e-9);
    EXPECT_NEAR((*result.combined)(4, 5), 0.05, 1.0e-7);
    EXPECT_EQ(samples.front().state.pathLength, 0.0);
    EXPECT_EQ(samples.back().state.pathLength, 0.1);
}

TEST_F(OrbitThreaderTest, RayTrackerZeroStepDoesNotEvaluateFields) {
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    ExternalFieldRayTracker tracker(beamline, reference);
    ExternalFieldRayTracker::State initial{
            Vector_t<double, 3>(1.0, 2.0, 3.0), Vector_t<double, 3>(0.0, 0.0, 1.0), 2.0};
    const auto step = tracker.step(initial, 0.0, [](const auto&, auto&, auto&) {
        ADD_FAILURE() << "A zero step must not evaluate a field";
        return false;
    });
    const Vector_t<double, 3> displacement = step.end.position - initial.position;
    const Vector_t<double, 3> momentumChange = step.end.momentum - initial.momentum;
    EXPECT_EQ(euclidean_norm(displacement), 0.0);
    EXPECT_EQ(euclidean_norm(momentumChange), 0.0);
    EXPECT_EQ(step.end.time, initial.time);
}

TEST_F(OrbitThreaderTest, RayTrackerInitializesAllFieldComponents) {
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    ExternalFieldRayTracker tracker(beamline, reference);
    ExternalFieldRayTracker::State initial;
    initial.momentum(2) = 1.0;
    const auto step = tracker.step(initial, 1.0e-11, [](const auto&, auto& electric, auto& magnetic) {
        for (unsigned component = 0; component < 3; ++component) {
            EXPECT_EQ(electric(component), 0.0);
            EXPECT_EQ(magnetic(component), 0.0);
        }
        return false;
    });
    EXPECT_EQ(step.end.position(0), 0.0);
    EXPECT_EQ(step.end.position(1), 0.0);
    EXPECT_EQ(step.end.momentum(0), 0.0);
    EXPECT_EQ(step.end.momentum(1), 0.0);
    EXPECT_EQ(step.end.momentum(2), 1.0);
}

TEST_F(OrbitThreaderTest, BishopFrameRemainsOrthonormalThroughBendAndReversal) {
    OpalBeamline beamline;
    auto frame = LinearTransferMapBuilder::initialFrame(
            beamline, Vector_t<double, 3>(0.0, 0.0, 1.0));
    for (const auto& momentum : {Vector_t<double, 3>(1.0, 2.0, 3.0),
                                 Vector_t<double, 3>(-1.0, -2.0, -3.0)}) {
        frame = LinearTransferMapBuilder::transportFrame(frame, momentum);
        EXPECT_NEAR(euclidean_norm(frame.xAxis), 1.0, 1.0e-14);
        EXPECT_NEAR(euclidean_norm(frame.yAxis), 1.0, 1.0e-14);
        EXPECT_NEAR(euclidean_norm(frame.sAxis), 1.0, 1.0e-14);
        EXPECT_NEAR(dot(frame.xAxis, frame.yAxis), 0.0, 1.0e-14);
        EXPECT_NEAR(dot(frame.xAxis, frame.sAxis), 0.0, 1.0e-14);
        EXPECT_NEAR(dot(frame.yAxis, frame.sAxis), 0.0, 1.0e-14);
    }
}

TEST_F(OrbitThreaderTest, RayTrackerResolvesThinSupportForEachMomentumAndReverseTime) {
    auto bunch = makeBunch(0);
    DummyBeamline line;
    DefaultVisitor visitor(line, false, false);
    OpalBeamline beamline;
    // Both the entrance and exit are beyond the nominal drift midpoint. The body
    // has zero length: only field support can prevent this region being skipped.
    FieldSupportOnlyComponent field("THIN_E", 0.003, 0.0034, 1.0e8);
    field.setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector_t<double, 3>(0.0), Quaternion()));
    field.fixPosition();
    beamline.visit(field, visitor, *bunch);
    beamline.prepareSections();
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    ExternalFieldRayTracker tracker(beamline, reference);
    for (double momentum : {0.9, 1.0, 1.1}) {
        ExternalFieldRayTracker::State initial;
        initial.momentum(2) = momentum;
        std::vector<ExternalFieldRayTracker::Step> steps;
        const auto final = tracker.advance(initial, 2.0e-11, &steps);
        ASSERT_GT(final.position(2), 0.0034);
        ASSERT_GT(steps.size(), 1);
        // Exact work-energy identity: gamma_out - gamma_in = q E L / (m c^2).
        // 1e-10 in gamma is much smaller than the 4.26e-5 signal; it allows the
        // second-order orbit error and floating-point boundary localization.
        const double expectedGamma = std::sqrt(1.0 + momentum * momentum)
                                     + 1.0e8 * 0.0004 / reference.getM();
        EXPECT_NEAR(std::sqrt(1.0 + dot(final.momentum, final.momentum)), expectedGamma, 1.0e-10);
        EXPECT_NEAR(final.time, 2.0e-11, 1.0e-24);
        const auto recovered = tracker.advance(final, -2.0e-11);
        const Vector_t<double, 3> displacement = recovered.position - initial.position;
        const Vector_t<double, 3> momentumChange = recovered.momentum - initial.momentum;
        EXPECT_LT(euclidean_norm(displacement), 1.0e-10);
        EXPECT_LT(euclidean_norm(momentumChange), 1.0e-9);
    }
}

TEST_F(OrbitThreaderTest, RayTrackerSumsOverlappingSupportsAcrossBoundaries) {
    auto bunch = makeBunch(0);
    DummyBeamline line;
    DefaultVisitor visitor(line, false, false);
    OpalBeamline beamline;
    FieldSupportOnlyComponent first("E1", 0.001, 0.003, 1.0e8);
    FieldSupportOnlyComponent second("E2", 0.002, 0.0035, -2.0e7);
    for (auto* field : {&first, &second}) {
        field->setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector_t<double, 3>(0.0), Quaternion()));
        field->fixPosition();
        beamline.visit(*field, visitor, *bunch);
    }
    beamline.prepareSections();
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    ExternalFieldRayTracker tracker(beamline, reference);
    ExternalFieldRayTracker::State initial;
    initial.momentum(2) = 1.0;
    const auto final = tracker.advance(initial, 2.0e-11);
    ASSERT_GT(final.position(2), 0.0035);
    const double work = 1.0e8 * 0.002 - 2.0e7 * 0.0015;
    EXPECT_NEAR(std::sqrt(1.0 + dot(final.momentum, final.momentum)),
                std::sqrt(2.0) + work / reference.getM(), 1.0e-10);
}

TEST_F(OrbitThreaderTest, ExecutesOverlapAndRecordsBothElements) {
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

    Options::enableLinearTransferMaps = true;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    OrbitThreader threader(
            reference, Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0, 0.0, 1.0), 0.0, 0.0, 0.0,
            1.0e-11, stepSizes, beamline, /*isDesignBeam=*/true);

    threader.execute();

    // The two quadrupoles overlap on s in [0.45, 0.5]; a query centred there returns both.
    IndexMap::value_t elements = threader.query(0.475, 0.01);
    std::set<std::string> activeNames;
    for (const auto& element : elements) {
        activeNames.insert(element->getName());
    }
    EXPECT_EQ(activeNames, (std::set<std::string>{"Q_LONG", "Q_SHORT"}));

    std::shared_ptr<ElementBase> runtimeLong;
    std::shared_ptr<ElementBase> runtimeShort;
    for (const auto& element : beamline.getElements()) {
        if (element->getName() == "Q_LONG") runtimeLong = element;
        if (element->getName() == "Q_SHORT") runtimeShort = element;
    }
    ASSERT_NE(runtimeLong, nullptr);
    ASSERT_NE(runtimeShort, nullptr);
    EXPECT_TRUE(runtimeLong->isOverlapping());
    EXPECT_TRUE(runtimeShort->isOverlapping());
    ASSERT_EQ(runtimeLong->getLinearTransferMaps().size(), 2);
    ASSERT_EQ(runtimeShort->getLinearTransferMaps().size(), 2);

    const auto findOverlap = [](const auto& maps) -> const LinearTransferMap* {
        for (const auto& map : maps) {
            if (map.includesOverlappingFields) return &map;
        }
        return nullptr;
    };
    const auto* longOverlap  = findOverlap(runtimeLong->getLinearTransferMaps());
    const auto* shortOverlap = findOverlap(runtimeShort->getLinearTransferMaps());
    ASSERT_NE(longOverlap, nullptr);
    ASSERT_NE(shortOverlap, nullptr);
    EXPECT_EQ(longOverlap->segment, shortOverlap->segment);
    EXPECT_EQ(
            std::set<std::string>(
                    longOverlap->activeElements.begin(), longOverlap->activeElements.end()),
            (std::set<std::string>{"Q_LONG", "Q_SHORT"}));

    const auto ordered = beamline.getLinearTransferMapsInReferenceOrder();
    EXPECT_EQ(ordered.size(), 3);  // shared overlap segment is returned only once
    ASSERT_TRUE(threader.getCombinedLinearTransferMap().has_value());
    ASSERT_TRUE(threader.getCombinedDeterminantResidual().has_value());
    ASSERT_TRUE(threader.getCombinedSymplecticResidual().has_value());
}

TEST_F(OrbitThreaderTest, BuildsOnePeriodicTurnIndependentOfTrackStepBudget) {
    auto bunch = makeBunch(0);

    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    OpalBeamline beamline;
    DriftRep drift("D_PERIODIC");
    drift.getGeometry().setElementLength(0.8);
    drift.setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector_t<double, 3>(0.0), Quaternion()));
    drift.fixPosition();
    beamline.visit(drift, visitor, *bunch);
    beamline.prepareSections();

    StepSizeConfig stepSizes;
    stepSizes.push_back(1.0e-11, 10.0, 1);
    stepSizes.resetIterator();

    PartData reference(1.0, 9.382720813e8, 1.0e6);
    Options::enableLinearTransferMaps = true;
    OrbitThreader threader(
            reference, Vector_t<double, 3>(0.0, 0.0, 0.2), Vector_t<double, 3>(0.0, 0.0, 1.0), 0.2,
            /*maxDiffZBunch=*/0.1, 0.0, 1.0e-11, stepSizes, beamline, /*isDesignBeam=*/true,
            /*period=*/0.6);

    threader.execute();

    EXPECT_EQ(names(threader.query(0.1, 0.01)), (std::set<std::string>{"D_PERIODIC"}));
    EXPECT_EQ(names(threader.query(0.7, 0.01)), (std::set<std::string>{"D_PERIODIC"}));
    ASSERT_TRUE(threader.getCombinedLinearTransferMap().has_value());
    EXPECT_NEAR((*threader.getCombinedLinearTransferMap())(0, 1), 0.6, 2.0e-5);
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
            1.0e-12, stepSizes, beamline, /*isDesignBeam=*/true);

    EXPECT_NO_THROW(threader.execute());
}

TEST_F(OrbitThreaderTest, CalculatesAndAttachesLinearDriftMap) {
    auto bunch = makeBunch(0);

    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);

    OpalBeamline beamline;
    DriftRep drift("D_MAP");
    constexpr double length = 0.3;
    drift.getGeometry().setElementLength(length);
    drift.setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector_t<double, 3>(0.0), Quaternion()));
    drift.fixPosition();
    beamline.visit(drift, visitor, *bunch);
    beamline.prepareSections();

    StepSizeConfig stepSizes;
    // MAXSTEPS limits production tracking, but an enabled map calculation must still thread the
    // complete requested interval to ZSTOP.
    stepSizes.push_back(1.0e-11, 0.4, 1);
    stepSizes.resetIterator();

    Options::enableLinearTransferMaps = true;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    OrbitThreader threader(
            reference, Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0, 0.0, 1.0), 0.0, 0.0, 0.0,
            1.0e-11, stepSizes, beamline, /*isDesignBeam=*/true);

    ASSERT_NO_THROW(threader.execute());
    const auto runtimeElements = beamline.getElements();
    ASSERT_EQ(runtimeElements.size(), 1);
    const auto& maps = (*runtimeElements.begin())->getLinearTransferMaps();
    ASSERT_EQ(maps.size(), 1);
    const auto& map = maps.front().matrix;
    EXPECT_NEAR(map(0, 0), 1.0, 1.0e-8);
    EXPECT_NEAR(map(0, 1), length, 2.0e-5);
    EXPECT_NEAR(map(2, 2), 1.0, 1.0e-8);
    EXPECT_NEAR(map(2, 3), length, 2.0e-5);
    EXPECT_NEAR(map(4, 5), length / 2.0, 2.0e-5);  // gamma^2 = 1 + |beta*gamma|^2 = 2
    EXPECT_NEAR(map(5, 5), 1.0, 1.0e-8);

    ASSERT_TRUE(threader.getCombinedLinearTransferMap().has_value());
    // The combined LINE map covers the full threaded interval, including the field-free tail to
    // ZSTOP; the element-owned map above covers only the drift itself.
    EXPECT_NEAR((*threader.getCombinedLinearTransferMap())(0, 1), 0.4, 1.0e-3);
    ASSERT_TRUE(threader.getCombinedDeterminantResidual().has_value());
    EXPECT_LT(*threader.getCombinedDeterminantResidual(), 1.0e-8);
    ASSERT_TRUE(threader.getCombinedSymplecticResidual().has_value());
    EXPECT_LT(*threader.getCombinedSymplecticResidual(), 1.0e-8);
}
