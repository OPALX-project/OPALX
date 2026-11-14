#include "gtest/gtest.h"

#include "AbstractObjects/OpalData.h"
#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/PartData.h"
#include "Algorithms/ReferenceTracer.h"
#include "BeamlineCore/MultipoleRep.h"
#include "BeamlineGeometry/Geometry.h"
#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"
#include "Utility/Inform.h"

#include <filesystem>
#include <fstream>
#include <string>

extern Inform* gmsg;

namespace {
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

    /// Zero body length, finite field extent: the length check must use the field extent.
    class FieldSupportOnlyComponent final : public ElementBase {
    public:
        FieldSupportOnlyComponent(const std::string& name, const double fieldLength)
            : ElementBase(name), fieldLength_m(fieldLength) {}

        void accept(BeamlineVisitor&) const override {}
        ElementBase* clone() const override { return new FieldSupportOnlyComponent(*this); }
        void apply(const std::shared_ptr<ParticleContainer_t>&, const StepTime& /*time*/) override {
        }
        void apply(
                const Vector_t<double, 3>&, const Vector_t<double, 3>&, const double&,
                Vector_t<double, 3>&, Vector_t<double, 3>&) override {}
        bool applyToReferenceParticle(
                const Vector_t<double, 3>&, const Vector_t<double, 3>&, const double&,
                Vector_t<double, 3>&, Vector_t<double, 3>&) override {
            return false;
        }
        void initialise(PartBunch_t*) override {}
        void finalise() override {}
        void getFieldExtent(double& zBegin, double& zEnd) const override {
            zBegin = 0.0;
            zEnd   = fieldLength_m;
        }
        ElementType getType() const override { return ElementType::ANY; }
        Geometry& getGeometry() override { return geometry_m; }
        const Geometry& getGeometry() const override { return geometry_m; }

    private:
        double fieldLength_m;
        Geometry geometry_m{Geometry::makeNull()};
    };
}  // namespace

class ReferenceTracerTest : public ::testing::Test {
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
        OpalData::getInstance()->storeInputFn("TestReferenceTracer.opal");
        OpalData::getInstance()->setOpenMode(OpalData::OpenMode::WRITE);
        std::filesystem::create_directories(OpalData::getInstance()->getAuxiliaryOutputDirectory());
        bunch_m = makeBunch();
    }

    class TestableFieldSolverCmd : public FieldSolverCmd {
    public:
        void setType(const std::string& t) {
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::TYPE], t);
        }
        void enableParallelDecomposition() {
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTX], true);
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTY], true);
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTZ], true);
        }
        void setBC(const std::string& bc) {
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTX], bc);
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTY], bc);
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTZ], bc);
        }
    };

    std::shared_ptr<PartBunch_t> makeBunch() {
        dataSink_m       = std::make_shared<DataSink>();
        const auto fsCmd = std::make_shared<TestableFieldSolverCmd>();
        fsCmd_m          = fsCmd;
        fsCmd->setType("NONE");
        fsCmd->enableParallelDecomposition();
        fsCmd->setNX(8);
        fsCmd->setNY(8);
        fsCmd->setNZ(8);
        fsCmd->setBC("PERIODIC");

        auto beam = std::make_shared<Beam>();
        Beam::defineBuiltins();

        Beam::defineBuiltins();
        Beam* opBeam = Beam::find("UNNAMED_BEAM");
        EXPECT_NE(opBeam, nullptr);
        auto bunch = std::make_shared<PartBunch_t>(
                std::vector{1.0}, std::vector{1.0}, std::vector<Beam*>{opBeam},
                std::vector<size_t>{0}, 1.0, "LF2",
                opalx::spacecharge::CartesianDomainConfig3D{.periodicParticleBoundary = true});
        return bunch;
    }

    static std::shared_ptr<MultipoleRep> makeQuadrupole(
            const std::string& name, const double x, const double length) {
        auto quadrupole = std::make_shared<MultipoleRep>(name);
        quadrupole->getGeometry().setElementLength(length);
        quadrupole->setCSTrafoGlobal2Local(
                CoordinateSystemTrafo(Vector_t<double, 3>(x, 0.0, 0.0), Quaternion()));
        quadrupole->fixPosition();
        return quadrupole;
    }

    void add(OpalBeamline& beamline, const ElementBase& element) {
        DummyBeamline line;
        DefaultVisitor visitor(line, false, false);
        beamline.visit(element, visitor, *bunch_m);
    }

    static ReferenceTracer::RayState launch() {
        ReferenceTracer::RayState start;
        start.momentum = Vector_t<double, 3>(0.0, 0.0, 1.0);
        return start;
    }

    static StepSizeConfig schedule(const double dt, const double zStop, const unsigned long steps) {
        StepSizeConfig stepSizes;
        stepSizes.push_back(dt, zStop, steps);
        stepSizes.resetIterator();
        return stepSizes;
    }

    PartData reference_m{1.0, 9.382720813e8, 1.0e6};
    std::shared_ptr<PartBunch_t> bunch_m;
    std::shared_ptr<FieldSolverCmd> fsCmd_m;
    std::shared_ptr<DataSink> dataSink_m;
};

TEST_F(ReferenceTracerTest, LengthCheckUsesTheFieldExtent) {
    OpalBeamline beamline;
    FieldSupportOnlyComponent component("TW_LIKE", 0.1);
    component.setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector_t<double, 3>(0.0), Quaternion()));
    component.fixPosition();
    add(beamline, component);
    beamline.prepareSections();

    // About 0.2 mm per step at beta gamma = 1: ten steps fit into the 0.1 m field extent.
    ReferenceTracer fine(reference_m, beamline, schedule(1.0e-12, 0.2, 1000), 1.0e-12, 0.0);
    EXPECT_NO_THROW(fine.execute(launch(), 0.0));

    // About 2.1 cm per step: ten steps do not fit.
    ReferenceTracer coarse(reference_m, beamline, schedule(1.0e-10, 0.2, 10), 1.0e-10, 0.0);
    EXPECT_THROW(coarse.execute(launch(), 0.0), OpalException);
}

TEST_F(ReferenceTracerTest, ReachCheckStopsWhenAParallelBranchWouldDeleteTheReference) {
    const auto traceWith = [this](const double branchReach) {
        OpalBeamline beamline;
        add(beamline, *makeQuadrupole("Q_HERE", 0.0, 0.5));
        auto branch = makeQuadrupole("Q_BRANCH", 0.3, 0.5);
        branch->setAperture(ApertureType::RECTANGULAR, {0.05, 0.05});
        if (branchReach > 0.0) branch->setReach(branchReach);
        add(beamline, *branch);
        beamline.prepareSections();
        ReferenceTracer tracer(reference_m, beamline, schedule(1.0e-11, 0.7, 512), 1.0e-11, 0.0);
        tracer.execute(launch(), 0.0);
    };

    EXPECT_THROW(traceWith(0.0), OpalException);
    EXPECT_NO_THROW(traceWith(0.2));
}

TEST_F(ReferenceTracerTest, EndsWhenTheReferenceLeavesTheLattice) {
    OpalBeamline beamline;
    add(beamline, *makeQuadrupole("Q", 0.0, 0.5));
    beamline.prepareSections();

    // Without ZSTOP and with a large step budget the pass still ends just after the element.
    ReferenceTracer tracer(reference_m, beamline, schedule(1.0e-11, 1.0e6, 1000000), 1.0e-11, 0.0);
    const auto end = tracer.execute(launch(), 0.0);
    EXPECT_GT(end.pathLength, 0.5);
    EXPECT_LT(end.pathLength, 0.52);
}

TEST_F(ReferenceTracerTest, WritesTheDesignPath) {
    OpalBeamline beamline;
    add(beamline, *makeQuadrupole("Q", 0.0, 0.5));
    beamline.prepareSections();
    ReferenceTracer tracer(reference_m, beamline, schedule(1.0e-11, 1.0e6, 1000), 1.0e-11, 0.0);
    tracer.execute(launch(), 0.0);
    if (ippl::Comm->rank() != 0) return;

    auto opal = OpalData::getInstance();
    std::ifstream file(
            Util::combineFilePath(
                    {opal->getAuxiliaryOutputDirectory(),
                     opal->getInputBasename() + "_DesignPath.dat"}));
    ASSERT_TRUE(file.good());
    std::string line;
    std::getline(file, line);
    EXPECT_EQ(line.front(), '#');
    unsigned rows = 0;
    while (std::getline(file, line))
        ++rows;
    EXPECT_GT(rows, 10u);
}

TEST(ReferenceTracerRay, HitsBoxOnlyAhead) {
    const BoundingBox box = BoundingBox::getBoundingBox(
            {Vector_t<double, 3>(0.0, 0.0, 0.0), Vector_t<double, 3>(1.0, 1.0, 1.0)});
    const Vector_t<double, 3> before(-1.0, 0.5, 0.5);
    EXPECT_TRUE(ReferenceTracer::rayHitsBox(box, before, Vector_t<double, 3>(1.0, 0.0, 0.0)));
    EXPECT_FALSE(ReferenceTracer::rayHitsBox(box, before, Vector_t<double, 3>(-1.0, 0.0, 0.0)));
    EXPECT_FALSE(ReferenceTracer::rayHitsBox(box, before, Vector_t<double, 3>(0.0, 1.0, 0.0)));
    EXPECT_TRUE(
            ReferenceTracer::rayHitsBox(
                    box, Vector_t<double, 3>(0.5, 0.5, 0.5), Vector_t<double, 3>(0.0, 0.0, -1.0)));
    EXPECT_FALSE(
            ReferenceTracer::rayHitsBox(
                    BoundingBox(), Vector_t<double, 3>(0.5), Vector_t<double, 3>(1.0, 0.0, 0.0)));
}
