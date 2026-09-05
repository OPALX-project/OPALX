#include "gtest/gtest.h"

#include "AbstractObjects/OpalData.h"
#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/OrbitThreader.h"
#include "Algorithms/PartData.h"
#include "BasicActions/Option.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/MultipoleRep.h"
#include "BeamlineGeometry/Geometry.h"
#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
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
        unsigned fieldCalls{0};
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
            ++fieldCalls;
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
        resetMapSettings();
        OpalData::getInstance()->storeInputFn("TestOrbitThreader.opal");
        OpalData::getInstance()->setOpenMode(OpalData::OpenMode::WRITE);
        std::filesystem::create_directories(OpalData::getInstance()->getAuxiliaryOutputDirectory());
    }

    void TearDown() override {
        Options::enableLinearTransferMaps = false;
        resetMapSettings();
    }

    static void resetMapSettings() {
        Options::linearTransferMapRichardsonLevels = 0;
        Options::linearTransferMapSteps.fill(1.e-3);
        Options::linearTransferMapIntegrator = "BORIS";
    }

    LinearTransferMap freeDriftMap(const unsigned levels, const double deltaStep,
                                  const double clockOrigin = 0.0, const std::string& method = "BORIS") {
        OpalBeamline beamline;
        PartData reference(1.0, 9.382720813e8, 1.0e6);
        auto entrance = LinearTransferMapBuilder::initialFrame(
                beamline, Vector_t<double, 3>(0.0, 0.0, 1.0));
        entrance.momentum = Vector_t<double, 3>(0.0, 0.0, 1.0);
        entrance.time = clockOrigin;
        auto exit         = entrance;
        exit.position(2) = exit.pathLength = 0.1;
        const double flightTime = 0.1 * std::sqrt(2.0) / Physics::c;
        exit.time = clockOrigin + flightTime;
        exit.timeCorrection = (exit.time - clockOrigin) - flightTime;
        LinearTransferMapBuilder::Settings settings;
        settings.richardsonLevels         = levels;
        settings.finiteDifferenceSteps[5] = deltaStep;
        settings.integrationMethod = ExternalFieldRayTracker::parseIntegrationMethod(method);
        LinearTransferMapBuilder builder(beamline, reference, 1.e-11, settings);
        return builder.build({{entrance}, {exit}}, 0.0).segments.front().map;
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
    // Standalone optics callers get explicit/default settings, not process-wide OPTION state.
    Options::linearTransferMapRichardsonLevels = 2;
    Options::linearTransferMapSteps.fill(0.1);
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
    EXPECT_EQ(result.segments.front().map.richardsonLevels, 0);
    EXPECT_EQ(result.segments.front().map.finiteDifferenceSteps[5], 1.e-3);
    EXPECT_FALSE(result.segments.front().map.richardsonCorrection.has_value());
    EXPECT_TRUE(result.segments.front().owners.empty());
    ASSERT_TRUE(result.combined);
    EXPECT_NEAR((*result.combined)(0, 1), 0.1, 1.0e-9);
    EXPECT_NEAR((*result.combined)(4, 5), 0.05, 1.0e-7);
    EXPECT_EQ(samples.front().state.pathLength, 0.0);
    EXPECT_EQ(samples.back().state.pathLength, 0.1);
}

TEST_F(OrbitThreaderTest, RichardsonHasExpectedDifferentiationOrder) {
    // Exact relativistic drift: R56=L/(1+p0^2)=0.05 for L=0.1 m and p0=1 beta*gamma.
    // Large delta amplitudes keep truncation above integration/roundoff noise in this test.
    for (unsigned levels = 0; levels <= 2; ++levels) {
        const auto coarse          = freeDriftMap(levels, 0.16);
        const auto fine            = freeDriftMap(levels, 0.08);
        const double coarseError   = std::abs(coarse.matrix(4, 5) - 0.05);
        const double fineError     = std::abs(fine.matrix(4, 5) - 0.05);
        const double expectedRatio = std::ldexp(1.0, 2 * int(levels + 1));
        EXPECT_GT(coarseError / fineError, 0.8 * expectedRatio);
        EXPECT_LT(coarseError / fineError, 1.2 * expectedRatio);
        EXPECT_EQ(fine.richardsonLevels, levels);
        EXPECT_EQ(fine.finiteDifferenceSteps[5], 0.08);
        EXPECT_EQ(fine.finestFiniteDifferenceSteps[5], std::ldexp(0.08, -int(levels)));
        EXPECT_EQ(fine.richardsonCorrection.has_value(), levels > 0);
        EXPECT_EQ(fine.integrationMethod, "BORIS");
    }
}

TEST_F(OrbitThreaderTest, RichardsonTableauMatchesIndependentDriftFormula) {
    // No exact map is used by the builder: the analytic flight time is only a test oracle.
    const auto zeta = [](const double delta) {
        const double p = 1.0 + delta;
        return -0.1 / std::sqrt(2.0) * std::sqrt(1.0 + p * p) / p;
    };
    std::vector<double> previous;
    for (unsigned level = 0; level <= LinearTransferMapBuilder::Settings::maximumRichardsonLevels;
         ++level) {
        const double h = std::ldexp(0.16, -int(level));
        std::vector<double> current{(zeta(h) - zeta(-h)) / (2.0 * h)};
        for (unsigned order = 1; order <= level; ++order) {
            const double factor = std::pow(4.0, order);
            current.push_back((factor * current.back() - previous[order - 1]) / (factor - 1.0));
        }
        const auto measured = freeDriftMap(level, 0.16);
        // Numerical tracking + exit-plane root finding add floating-point error to exact drift.
        EXPECT_NEAR(measured.matrix(4, 5), current.back(), 2.e-11);
        if (level > 0) {
            ASSERT_TRUE(measured.richardsonCorrection);
            EXPECT_NEAR(
                    (*measured.richardsonCorrection)[5], std::abs(current.back() - previous.back()),
                    4.e-11);
        }
        previous = std::move(current);
    }
}

TEST_F(OrbitThreaderTest, ValidatesMapSettingsAndIntegrationChoice) {
    LinearTransferMapBuilder::Settings settings;
    EXPECT_NO_THROW(settings.validate());
    for (double invalid :
         {0.0, -1.0, std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::quiet_NaN()}) {
        settings.finiteDifferenceSteps[0] = invalid;
        EXPECT_THROW(settings.validate(), OpalException);
    }
    settings                          = {};
    settings.finiteDifferenceSteps[5] = 1.0;
    EXPECT_THROW(settings.validate(), OpalException);
    settings                  = {};
    settings.richardsonLevels = 5;
    EXPECT_THROW(settings.validate(), OpalException);
    using Method = ExternalFieldRayTracker::IntegrationMethod;
    EXPECT_EQ(ExternalFieldRayTracker::parseIntegrationMethod("LF2"), Method::BORIS);
    EXPECT_EQ(ExternalFieldRayTracker::parseIntegrationMethod("RK4"), Method::RK4);
    EXPECT_EQ(ExternalFieldRayTracker::parseIntegrationMethod("DOP853"), Method::DOP853);
    EXPECT_THROW(ExternalFieldRayTracker::parseIntegrationMethod("UNKNOWN"), OpalException);
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    EXPECT_THROW(
            ExternalFieldRayTracker(beamline, reference, static_cast<Method>(123)), OpalException);
    EXPECT_THROW(LinearTransferMapBuilder(beamline, reference, 0.0), OpalException);
}

TEST_F(OrbitThreaderTest, MapOptionsPersistAcrossStatementsAndRejectInvalidInputs) {
    Option exemplar;
    std::unique_ptr<Option> command(exemplar.clone("MAP_OPTIONS"));
    Attributes::setReal(*command->findAttribute("LINEARTRANSFERMAPRICHARDSON"), 2.0);
    const std::vector<double> steps{1.e-3, 2.e-3, 3.e-3, 4.e-3, 5.e-3, 6.e-3};
    Attributes::setRealArray(*command->findAttribute("LINEARTRANSFERMAPSTEPS"), steps);
    Attributes::setPredefinedString(*command->findAttribute("LINEARTRANSFERMAPINTEGRATOR"), "LF2");
    command->execute();
    EXPECT_EQ(Options::linearTransferMapRichardsonLevels, 2);
    EXPECT_EQ(Options::linearTransferMapSteps[5], 6.e-3);
    EXPECT_EQ(Options::linearTransferMapIntegrator, "BORIS");
    std::unique_ptr<Option> next(exemplar.clone("NEXT_MAP_OPTIONS"));
    next->execute();  // An unrelated later OPTION must not reset these values.
    EXPECT_EQ(Options::linearTransferMapRichardsonLevels, 2);
    EXPECT_EQ(Attributes::getRealArray(*next->findAttribute("LINEARTRANSFERMAPSTEPS")), steps);
    for (double invalid :
         {-1.0, 1.5, 5.0, std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::quiet_NaN()}) {
        Attributes::setReal(*next->findAttribute("LINEARTRANSFERMAPRICHARDSON"), invalid);
        EXPECT_THROW(next->execute(), OpalException);
        EXPECT_EQ(Options::linearTransferMapRichardsonLevels, 2);
    }
    Attributes::setReal(*next->findAttribute("LINEARTRANSFERMAPRICHARDSON"), 1.0);
    for (const auto& invalid :
         {std::vector<double>{}, std::vector<double>{1.e-3}, std::vector<double>(7, 1.e-3),
          std::vector<double>(6, 0.0), std::vector<double>(6, 1.0)}) {
        Attributes::setRealArray(*next->findAttribute("LINEARTRANSFERMAPSTEPS"), invalid);
        EXPECT_THROW(next->execute(), OpalException);
        EXPECT_EQ(Options::linearTransferMapRichardsonLevels, 2);
    }
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

TEST_F(OrbitThreaderTest, RayTrackerRetainsSmallPositionTimeAndPathIncrements) {
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    ExternalFieldRayTracker tracker(beamline, reference);
    ExternalFieldRayTracker::State state;
    state.momentum(2) = 1.0;
    state.position = Vector_t<double, 3>(1.0, 2.0, 3.0);
    state.time = state.pathLength = 1.0;
    const auto initial = state;
    constexpr unsigned steps = 10000;
    constexpr double dt = 1.e-25;  // Each half-drift/time update is below a high-part ulp.
    for (unsigned i = 0; i < steps; ++i)
        state = tracker.step(state, dt, [](const auto&, auto&, auto&) { return false; }).end;
    const double distance = steps * dt * Physics::c / std::sqrt(2.0);
    // Difference includes the retained low part; 1e-25 m is ~1e-12 relative here.
    EXPECT_NEAR((state.position(2) - initial.position(2)) - state.positionCorrection(2), distance, 1.e-25);
    EXPECT_NEAR((state.pathLength - initial.pathLength) - state.pathLengthCorrection, distance, 1.e-25);
    EXPECT_NEAR((state.time - initial.time) - state.timeCorrection, steps * dt, 1.e-32);
    EXPECT_EQ(state.position(0), initial.position(0));
    EXPECT_EQ(state.position(1), initial.position(1));
}

TEST_F(OrbitThreaderTest, RayTrackerMagneticPathUsesSpeedNotChordAndReverses) {
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    ExternalFieldRayTracker tracker(beamline, reference);
    ExternalFieldRayTracker::State initial;
    initial.momentum(2) = 1.0;
    constexpr double dt = 1.e-9;
    const auto field = [](const auto&, auto&, auto& magnetic) {
        magnetic(1) = 10.0;
        return false;
    };
    const auto step = tracker.step(initial, dt, field);
    const double length = Physics::c * dt / std::sqrt(2.0);
    EXPECT_NEAR(step.midpoint.pathLength, 0.5 * length, 1.e-15);
    EXPECT_NEAR(step.end.pathLength, length, 1.e-15);
    const Vector_t<double, 3> chord = step.end.position - initial.position;
    EXPECT_GT(length - euclidean_norm(chord), 1.e-3);
    EXPECT_NEAR(euclidean_norm(step.end.momentum), 1.0, 1.e-15);
    const auto recovered = tracker.step(step.end, -dt, field).end;
    EXPECT_NEAR(euclidean_norm(recovered.position), 0.0, 1.e-15);
    EXPECT_NEAR(recovered.pathLength, 0.0, 1.e-15);
    EXPECT_NEAR(recovered.time, 0.0, 1.e-24);
}

TEST_F(OrbitThreaderTest, PathStopIsTrackedUnderAccelerationAndReverseTime) {
    auto bunch = makeBunch(0);
    DummyBeamline line;
    DefaultVisitor visitor(line, false, false);
    OpalBeamline beamline;
    FieldSupportOnlyComponent field("PATH_E", -1.0, 1.0, 1.e10);
    field.setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector_t<double, 3>(0.0), Quaternion()));
    field.fixPosition();
    beamline.visit(field, visitor, *bunch);
    beamline.prepareSections();
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    ExternalFieldRayTracker tracker(beamline, reference);
    ExternalFieldRayTracker::State initial;
    initial.momentum(2) = 0.1;
    constexpr double dt = 1.e-9;
    const auto end = tracker.advance(initial, dt);
    const double target = 0.37 * end.pathLength;
    const auto clipped = tracker.advanceToPathLength(initial, dt, target);
    EXPECT_NEAR(clipped.pathLength, target, 2.e-12);  // c times the existing time-root tolerance.
    EXPECT_NEAR(clipped.position(2), target, 2.e-12);
    EXPECT_GT(std::abs(clipped.time - 0.37 * dt), 1.e-11);  // A time fraction is not a path root.
    const auto reverse = tracker.advanceToPathLength(end, -dt, target);
    EXPECT_NEAR(reverse.pathLength, target, 2.e-12);
    EXPECT_THROW(tracker.advanceToPathLength(initial, dt, 2.0 * end.pathLength), OpalException);
    EXPECT_THROW(tracker.advanceToPathLength(initial, dt, -1.0), OpalException);
}

TEST_F(OrbitThreaderTest, CompensatedClockPreservesFullDriftMapAtLargeTimeOrigin) {
    const auto zeroOrigin = freeDriftMap(1, 0.03);
    const auto lateOrigin = freeDriftMap(1, 0.03, 1000.0);
    for (unsigned row = 0; row < 6; ++row)
        for (unsigned column = 0; column < 6; ++column)
            EXPECT_NEAR(zeroOrigin.matrix(row, column), lateOrigin.matrix(row, column), 2.e-11);
}

TEST_F(OrbitThreaderTest, AcceleratedPathQuadratureConvergesAtSecondOrder) {
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    ExternalFieldRayTracker tracker(beamline, reference);
    constexpr double electric = 1.e10, duration = 1.e-9, initialMomentum = 0.1;
    const double rate = reference.getQ() * electric * Physics::c / reference.getM();
    const double finalMomentum = initialMomentum + rate * duration;
    const double exactLength = Physics::c / rate *
            (std::sqrt(1.0 + finalMomentum * finalMomentum) - std::sqrt(1.0 + initialMomentum * initialMomentum));
    std::vector<double> errors;
    for (unsigned steps : {100u, 200u}) {
        ExternalFieldRayTracker::State ray;
        ray.momentum(2) = initialMomentum;
        for (unsigned i = 0; i < steps; ++i)
            ray = tracker.step(ray, duration / steps, [](const auto&, auto& field, auto&) {
                field(2) = electric;
                return false;
            }).end;
        EXPECT_NEAR(ray.momentum(2), finalMomentum, 1.e-12);
        errors.push_back(std::abs(ray.pathLength - exactLength));
    }
    EXPECT_NEAR(errors[0] / errors[1], 4.0, 0.01);
}

TEST_F(OrbitThreaderTest, MapBuilderClipsPrerollWithCompensatedReferenceClock) {
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    auto before = LinearTransferMapBuilder::initialFrame(
            beamline, Vector_t<double, 3>(0.0, 0.0, 1.0));
    before.momentum(2) = 1.0;
    before.position(2) = before.pathLength = -0.05;
    before.time = 1000.0;
    auto after = before;
    after.position(2) = after.pathLength = 0.05;
    const double flightTime = 0.1 * std::sqrt(2.0) / Physics::c;
    after.time = before.time + flightTime;
    after.timeCorrection = (after.time - before.time) - flightTime;
    LinearTransferMapBuilder builder(beamline, reference, 1.e-11);
    const auto result = builder.build({{before}, {after}}, 0.0);
    ASSERT_EQ(result.segments.size(), 1);
    const auto& map = result.segments.front().map;
    EXPECT_EQ(map.entrance.pathLength, 0.0);
    EXPECT_NEAR(map.entrance.position(2), 0.0, 2.e-13);
    const double elapsed = (map.entrance.time - before.time) - map.entrance.timeCorrection;
    EXPECT_NEAR(elapsed, 0.5 * flightTime, 1.e-21);
    EXPECT_NEAR(map.matrix(0, 1), 0.05, 1.e-10);
    EXPECT_NEAR(map.matrix(4, 5), 0.025, 1.e-7);  // Existing centered delta truncation.
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
    for (const std::string method : {"BORIS", "RK4", "DOP853"}) {
    SCOPED_TRACE(method);
    ExternalFieldRayTracker tracker(beamline, reference, ExternalFieldRayTracker::parseIntegrationMethod(method));
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
    for (const std::string method : {"BORIS", "RK4", "DOP853"}) {
    SCOPED_TRACE(method);
    ExternalFieldRayTracker tracker(beamline, reference, ExternalFieldRayTracker::parseIntegrationMethod(method));
    ExternalFieldRayTracker::State initial;
    initial.momentum(2) = 1.0;
    const auto final = tracker.advance(initial, 2.0e-11);
    ASSERT_GT(final.position(2), 0.0035);
    const double work = 1.0e8 * 0.002 - 2.0e7 * 0.0015;
    EXPECT_NEAR(std::sqrt(1.0 + dot(final.momentum, final.momentum)),
                std::sqrt(2.0) + work / reference.getM(), 1.0e-10);
    }
}

TEST_F(OrbitThreaderTest, RungeKuttaMagneticHelixHasExpectedGlobalOrder) {
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    constexpr double magnetic = 10.0, transverse = 0.2, longitudinal = 0.1;
    const double gamma = std::sqrt(1.0 + transverse * transverse + longitudinal * longitudinal);
    const double omega = reference.getQ() * Physics::c * Physics::c * magnetic / (reference.getM() * gamma);
    const double duration = 4.0 / omega;
    for (const std::string method : {"RK4", "DOP853"}) {
        SCOPED_TRACE(method);
        ExternalFieldRayTracker tracker(beamline, reference, ExternalFieldRayTracker::parseIntegrationMethod(method));
        std::vector<double> errors;
        for (unsigned n : {4u, 8u, 16u}) {
            ExternalFieldRayTracker::State ray;
            ray.momentum = Vector_t<double, 3>(transverse, 0.0, longitudinal);
            for (unsigned i = 0; i < n; ++i)
                ray = tracker.step(ray, duration / n, [](const auto&, auto&, auto& b) {
                    b(2) = magnetic;
                    return false;
                }).end;
            const Vector_t<double, 3> exactP(transverse * std::cos(4.0), -transverse * std::sin(4.0), longitudinal);
            const Vector_t<double, 3> exactR(
                    Physics::c * transverse / (gamma * omega) * std::sin(4.0),
                    Physics::c * transverse / (gamma * omega) * (std::cos(4.0) - 1.0),
                    Physics::c * longitudinal / gamma * duration);
            const Vector_t<double, 3> dp = ray.momentum - exactP;
            const Vector_t<double, 3> dr = ray.position - exactR;
            const double exactPath = Physics::c * std::hypot(transverse, longitudinal) / gamma * duration;
            errors.push_back(std::max({euclidean_norm(dp), euclidean_norm(dr), std::abs(ray.pathLength - exactPath)}));
            std::cout << std::scientific << std::setprecision(12)
                      << "HELIX " << method << " steps=" << n << " dt_s=" << duration / n
                      << " position_error_m=" << euclidean_norm(dr)
                      << " momentum_error=" << euclidean_norm(dp)
                      << " path_error_m=" << std::abs(ray.pathLength - exactPath)
                      << " gamma_error=" << std::abs(std::sqrt(1.0 + dot(ray.momentum, ray.momentum)) - gamma)
                      << std::defaultfloat << '\n';
        }
        const double order = std::log2(errors[1] / errors[2]);
        EXPECT_GT(order, method == "RK4" ? 3.7 : 7.4);
        EXPECT_LT(order, method == "RK4" ? 4.4 : 8.6);
    }
}

TEST_F(OrbitThreaderTest, RungeKuttaTimeDependentElectricForceUsesStageClocksAndChargeUnits) {
    OpalBeamline beamline;
    constexpr double duration = 1.e-9, amplitude = 3.e8;
    for (double charge : {-1.0, 1.0}) {
        PartData reference(charge, 9.382720813e8, 1.0e6);
        for (const std::string method : {"RK4", "DOP853"}) {
            SCOPED_TRACE(method);
            ExternalFieldRayTracker tracker(beamline, reference, ExternalFieldRayTracker::parseIntegrationMethod(method));
            ExternalFieldRayTracker::State ray;
            ray.momentum(2) = 1.0;
            const auto fields = [](const auto& state, auto& e, auto&) {
                e(2) = amplitude * std::cos(state.time / duration);
                return false;
            };
            for (unsigned i = 0; i < 20; ++i) ray = tracker.step(ray, duration / 20, fields).end;
            const double exact = 1.0 + charge * Physics::c * amplitude / reference.getM() * duration * std::sin(1.0);
            EXPECT_NEAR(ray.momentum(2), exact, method == "RK4" ? 1.e-9 : 3.e-15);
            for (unsigned i = 0; i < 20; ++i) ray = tracker.step(ray, -duration / 20, fields).end;
            EXPECT_NEAR(ray.momentum(2), 1.0, 3.e-14);
            EXPECT_NEAR(ray.position(2), 0.0, 1.e-10);
            EXPECT_NEAR(ray.pathLength, 0.0, 1.e-10);
        }
    }
}

TEST_F(OrbitThreaderTest, RungeKuttaZeroStepsMaterialHitsAndFieldInitialization) {
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    for (const std::string method : {"RK4", "DOP853"}) {
        SCOPED_TRACE(method);
        ExternalFieldRayTracker tracker(beamline, reference, ExternalFieldRayTracker::parseIntegrationMethod(method));
        ExternalFieldRayTracker::State ray;
        ray.momentum(2) = 1.0;
        unsigned calls = 0;
        const auto field = [&](const auto&, auto& e, auto& b) {
            ++calls;
            for (unsigned d = 0; d < 3; ++d) { EXPECT_EQ(e(d), 0.0); EXPECT_EQ(b(d), 0.0); }
            e(2) = 1.0;
            return false;
        };
        tracker.step(ray, 0.0, field);
        EXPECT_EQ(calls, 0u);
        tracker.step(ray, 1.e-11, field);
        EXPECT_EQ(calls, method == "RK4" ? 9u : 25u);
        EXPECT_TRUE(tracker.step(ray, 1.e-11, [](const auto&, auto&, auto&) { return true; }).hitMaterial);
        EXPECT_THROW(tracker.step(ray, std::numeric_limits<double>::quiet_NaN(), field), OpalException);
    }
}

TEST_F(OrbitThreaderTest, RungeKuttaMapsPreserveDriftAndCompensatedClock) {
    for (const std::string method : {"RK4", "DOP853"}) {
        SCOPED_TRACE(method);
        Option exemplar;
        std::unique_ptr<Option> option(exemplar.clone("RK_MAP_OPTION"));
        Attributes::setPredefinedString(*option->findAttribute("LINEARTRANSFERMAPINTEGRATOR"), method);
        option->execute();
        EXPECT_EQ(Options::linearTransferMapIntegrator, method);
        const auto origin = freeDriftMap(1, .03, 0.0, method);
        const auto late = freeDriftMap(1, .03, 1000.0, method);
        EXPECT_EQ(origin.integrationMethod, method);
        EXPECT_NEAR(origin.matrix(0, 1), .1, 1.e-10);
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                EXPECT_NEAR(origin.matrix(i, j), late.matrix(i, j), 2.e-11);
    }
}

TEST_F(OrbitThreaderTest, MapIntegratorOptionDoesNotAffectDisabledMapsOrSecondaryThreading) {
    auto bunch = makeBunch(0);
    DummyBeamline line;
    DefaultVisitor visitor(line, false, false);
    OpalBeamline beamline;
    FieldSupportOnlyComponent field("ISOLATION_E", 0.0, .1, 1.e8);
    field.setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector_t<double, 3>(0.0), Quaternion()));
    field.fixPosition();
    beamline.visit(field, visitor, *bunch);
    beamline.prepareSections();
    auto* runtime = dynamic_cast<FieldSupportOnlyComponent*>(beamline.getElements().begin()->get());
    ASSERT_NE(runtime, nullptr);
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    for (bool design : {false, true}) {
        // Secondary threading remains Boris even when map computation is enabled.
        Options::enableLinearTransferMaps = !design;
        unsigned borisCalls = 0;
        for (const std::string method : {"BORIS", "RK4", "DOP853"}) {
            SCOPED_TRACE(method);
            Options::linearTransferMapIntegrator = method;
            runtime->fieldCalls = 0;
            StepSizeConfig steps;
            steps.push_back(1.e-11, .1, 20);
            steps.resetIterator();
            OrbitThreader threader(reference, Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0, 0.0, 1.0),
                    0.0, 0.0, 0.0, 1.e-11, steps, beamline, design);
            threader.execute();
            EXPECT_FALSE(threader.getCombinedLinearTransferMap());
            EXPECT_GT(runtime->fieldCalls, 0u);
            if (method == "BORIS") borisCalls = runtime->fieldCalls;
            else EXPECT_EQ(runtime->fieldCalls, borisCalls);
        }
    }
}

TEST_F(OrbitThreaderTest, ExecutesOverlapAndRecordsBothElements) {
    Options::linearTransferMapRichardsonLevels = 1;
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
    EXPECT_EQ(longOverlap->richardsonLevels, 1);
    EXPECT_EQ(shortOverlap->richardsonLevels, 1);
    EXPECT_EQ(longOverlap->finestFiniteDifferenceSteps[0], 5.e-4);
    EXPECT_TRUE(longOverlap->richardsonCorrection.has_value());
    for (unsigned row = 0; row < 6; ++row)
        for (unsigned column = 0; column < 6; ++column)
            EXPECT_EQ(longOverlap->matrix(row, column), shortOverlap->matrix(row, column));
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
    Options::linearTransferMapRichardsonLevels = 1;
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
    const auto maps = beamline.getLinearTransferMapsInReferenceOrder();
    ASSERT_FALSE(maps.empty());
    EXPECT_EQ(maps.front().map->richardsonLevels, 1);
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
