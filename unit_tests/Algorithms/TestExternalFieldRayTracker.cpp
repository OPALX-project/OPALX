#include "../BunchBeam.h"
#include "../TemporaryDirectory.h"
#include "gtest/gtest.h"

#include "AbsBeamline/Drift.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/SBend.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/ExternalFieldRayTracker.h"
#include "Algorithms/PartData.h"
#include "BeamlineGeometry/Geometry.h"
#include "Elements/OpalBeamline.h"
#include "Structure/Beam.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

#include <array>
#include <bit>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>

extern Inform* gmsg;

namespace {
    // Freeze the pre-cache traversal as a regression oracle. Public step() still
    // supplies the numerical method; this independently recomputes every support
    // set and checks that reuse changes neither the trial tree nor any ray bits.
    ExternalFieldRayTracker::State advanceWithoutMembershipReuse(
            ExternalFieldRayTracker& tracker, OpalBeamline& beamline,
            const ExternalFieldRayTracker::State& initial, double dt,
            std::vector<ExternalFieldRayTracker::Step>& accepted) {
        using Ray  = ExternalFieldRayTracker::State;
        double cap = std::numeric_limits<double>::max();
        for (const auto& element : beamline.getElements()) {
            if (element->getType() == ElementType::MARKER
                || element->getType() == ElementType::MONITOR)
                continue;
            double begin = 0.0, end = 0.0;
            element->getFieldExtent(begin, end);
            if (std::abs(end - begin) > 0.0)
                cap = std::min(cap, std::abs(end - begin) / (4 * Physics::c));
            const double length = element->getGeometry().getArcLength();
            if (length > 0.0) cap = std::min(cap, length / (4 * Physics::c));
        }
        const double tolerance = 1.0e-12 * std::abs(dt);
        std::function<Ray(const Ray&, double, unsigned)> advance;
        advance = [&](const Ray& start, double h, unsigned depth) -> Ray {
            const double floor = 64.0 * std::numeric_limits<double>::epsilon()
                                 * std::max(1.0, euclidean_norm(start.position)) / Physics::c;
            const bool resolved = std::abs(h) <= std::max(tolerance, floor);
            const auto split    = [&]() {
                if (depth >= 64) throw std::runtime_error("Uncached oracle did not converge");
                const Ray middle = advance(start, 0.5 * h, depth + 1);
                return advance(middle, 0.5 * h, depth + 1);
            };
            if (!resolved && std::abs(h) > cap) return split();
            const auto initialSet = beamline.getElements(start.position);
            bool crossed          = false;
            const auto trial =
                    tracker.step(start, h, [&](const Ray& ray, auto& electric, auto& magnetic) {
                        const auto stageSet = beamline.getElements(ray.position);
                        crossed             = crossed || stageSet != initialSet;
                        for (const auto& element : stageSet) {
                            if (element->getType() == ElementType::MARKER
                                || element->getType() == ElementType::MONITOR)
                                continue;
                            const auto localR = beamline.transformToLocalCS(element, ray.position);
                            const auto localP = beamline.rotateToLocalCS(element, ray.momentum);
                            Vector_t<double, 3> localE(0.0), localB(0.0);
                            if (element->applyToReferenceParticle(
                                        localR, localP, ray.time, localE, localB))
                                return true;
                            electric += beamline.rotateFromLocalCS(element, localE);
                            magnetic += beamline.rotateFromLocalCS(element, localB);
                        }
                        return false;
                    });
            if (!resolved && (crossed || initialSet != beamline.getElements(trial.end.position)))
                return split();
            if (trial.hitMaterial) throw std::runtime_error("Uncached oracle hit material");
            accepted.push_back(trial);
            return trial.end;
        };
        return advance(initial, dt, 0);
    }

    void expectSameRayBits(
            const ExternalFieldRayTracker::State& a, const ExternalFieldRayTracker::State& b) {
        const auto same = [](double x, double y) {
            EXPECT_EQ(std::bit_cast<std::uint64_t>(x), std::bit_cast<std::uint64_t>(y));
        };
        for (unsigned d = 0; d < 3; ++d) {
            same(a.position(d), b.position(d));
            same(a.positionCorrection(d), b.positionCorrection(d));
            same(a.momentum(d), b.momentum(d));
        }
        same(a.time, b.time);
        same(a.timeCorrection, b.timeCorrection);
        same(a.pathLength, b.pathLength);
        same(a.pathLengthCorrection, b.pathLengthCorrection);
    }

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
                const double longitudinalElectricField = 0.0, const double magneticField = 0.0)
            : ElementBase(name),
              fieldBegin_m(fieldBegin),
              fieldEnd_m(fieldEnd),
              electricField_m(longitudinalElectricField),
              magneticField_m(magneticField) {
            getGeometry() = Geometry::makeNull();
        }

        void apply(const std::shared_ptr<ParticleContainer_t>&, const StepTime& /*time*/) override {
        }

        void apply(
                const Vector_t<double, 3>&, const Vector_t<double, 3>&, const double&,
                Vector_t<double, 3>&, Vector_t<double, 3>&) override {}

        bool applyToReferenceParticle(
                const Vector_t<double, 3>& position, const Vector_t<double, 3>&, const double&,
                Vector_t<double, 3>& electric, Vector_t<double, 3>& magnetic) override {
            ++fieldCalls;
            if (position(2) >= fieldBegin_m && position(2) < fieldEnd_m) {
                electric(2) += electricField_m;
                magnetic(1) += magneticField_m;
            }
            return false;
        }

        void initialise(PartBunch_t*) override {}

        void getFieldExtent(double& zBegin, double& zEnd) const override {
            zBegin = fieldBegin_m;
            zEnd   = fieldEnd_m;
        }

        ElementType getType() const override { return ElementType::ANY; }

    private:
        double fieldBegin_m;
        double fieldEnd_m;
        double electricField_m;
        double magneticField_m;
    };
}  // namespace

class ExternalFieldRayTrackerTest : public ::testing::Test {
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
        savedInput_m = OpalData::getInstance()->getInputFn();
        directory_m.create("opalx-external-field-ray-tracker", ippl::Comm->getCommunicator());
        OpalData::getInstance()->storeInputFn(
                (directory_m.path() / "TestExternalFieldRayTracker.opal").string());
        OpalData::getInstance()->setOpenMode(OpalData::OpenMode::WRITE);
    }

    void TearDown() override {
        OpalData::getInstance()->storeInputFn(savedInput_m);
        ippl::Comm->barrier();
        directory_m.reset();
    }

    opalx::test::TemporaryDirectory directory_m;
    std::string savedInput_m;

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

        auto beam = std::make_shared<Beam>();
        Beam::defineBuiltins();
        Beam* opBeam = Beam::find("UNNAMED_BEAM");
        EXPECT_NE(opBeam, nullptr);

        auto bunch = std::make_shared<PartBunch_t>(
                std::vector{opalx::test::bunchBeam(*opBeam, numParticles)}, 1.0,
                opalx::spacecharge::CartesianDomainConfig3D{.periodicParticleBoundary = true});
        bunch->getParticleContainer()->createParticles(numParticles);
        return bunch;
    }

    std::shared_ptr<FieldSolverCmd> fsCmdBase_m;
};

TEST_F(ExternalFieldRayTrackerTest, RayTrackerZeroStepDoesNotEvaluateFields) {
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    ExternalFieldRayTracker tracker(beamline, reference);
    ExternalFieldRayTracker::State initial{
            Vector_t<double, 3>(1.0, 2.0, 3.0), Vector_t<double, 3>(0.0, 0.0, 1.0), 2.0};
    const auto step = tracker.step(initial, 0.0, [](const auto&, auto&, auto&) {
        ADD_FAILURE() << "A zero step must not evaluate a field";
        return false;
    });
    const Vector_t<double, 3> displacement   = step.end.position - initial.position;
    const Vector_t<double, 3> momentumChange = step.end.momentum - initial.momentum;
    EXPECT_EQ(euclidean_norm(displacement), 0.0);
    EXPECT_EQ(euclidean_norm(momentumChange), 0.0);
    EXPECT_EQ(step.end.time, initial.time);
}

TEST_F(ExternalFieldRayTrackerTest, RayTrackerInitializesAllFieldComponents) {
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    ExternalFieldRayTracker tracker(beamline, reference);
    ExternalFieldRayTracker::State initial;
    initial.momentum(2) = 1.0;
    const auto step =
            tracker.step(initial, 1.0e-11, [](const auto&, auto& electric, auto& magnetic) {
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

TEST_F(ExternalFieldRayTrackerTest, RayTrackerRetainsSmallPositionTimeAndPathIncrements) {
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    ExternalFieldRayTracker tracker(beamline, reference);
    ExternalFieldRayTracker::State state;
    state.momentum(2) = 1.0;
    state.position    = Vector_t<double, 3>(1.0, 2.0, 3.0);
    state.time = state.pathLength = 1.0;
    const auto initial            = state;
    constexpr unsigned steps      = 10000;
    constexpr double dt = 1.e-25;  // Each half-drift/time update is below a high-part ulp.
    for (unsigned i = 0; i < steps; ++i)
        state = tracker.step(state, dt,
                             [](const auto&, auto&, auto&) {
                                 return false;
                             })
                        .end;
    const double distance = steps * dt * Physics::c / std::sqrt(2.0);
    // Difference includes the retained low part; 1e-25 m is ~1e-12 relative here.
    EXPECT_NEAR(
            (state.position(2) - initial.position(2)) - state.positionCorrection(2), distance,
            1.e-25);
    EXPECT_NEAR(
            (state.pathLength - initial.pathLength) - state.pathLengthCorrection, distance, 1.e-25);
    EXPECT_NEAR((state.time - initial.time) - state.timeCorrection, steps * dt, 1.e-32);
    EXPECT_EQ(state.position(0), initial.position(0));
    EXPECT_EQ(state.position(1), initial.position(1));
}

TEST_F(ExternalFieldRayTrackerTest, RayTrackerMagneticPathUsesSpeedNotChordAndReverses) {
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    ExternalFieldRayTracker tracker(beamline, reference);
    ExternalFieldRayTracker::State initial;
    initial.momentum(2) = 1.0;
    constexpr double dt = 1.e-9;
    const auto field    = [](const auto&, auto&, auto& magnetic) {
        magnetic(1) = 10.0;
        return false;
    };
    const auto step     = tracker.step(initial, dt, field);
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

TEST_F(ExternalFieldRayTrackerTest, AcceleratedPathQuadratureConvergesAtSecondOrder) {
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    ExternalFieldRayTracker tracker(beamline, reference);
    constexpr double electric = 1.e10, duration = 1.e-9, initialMomentum = 0.1;
    const double rate          = reference.getQ() * electric * Physics::c / reference.getM();
    const double finalMomentum = initialMomentum + rate * duration;
    const double exactLength   = Physics::c / rate
                               * (std::sqrt(1.0 + finalMomentum * finalMomentum)
                                  - std::sqrt(1.0 + initialMomentum * initialMomentum));
    std::vector<double> errors;
    for (unsigned steps : {100u, 200u}) {
        ExternalFieldRayTracker::State ray;
        ray.momentum(2) = initialMomentum;
        for (unsigned i = 0; i < steps; ++i)
            ray = tracker.step(ray, duration / steps,
                               [](const auto&, auto& field, auto&) {
                                   field(2) = electric;
                                   return false;
                               })
                          .end;
        EXPECT_NEAR(ray.momentum(2), finalMomentum, 1.e-12);
        errors.push_back(std::abs(ray.pathLength - exactLength));
    }
    EXPECT_NEAR(errors[0] / errors[1], 4.0, 0.01);
}

TEST_F(ExternalFieldRayTrackerTest, RayTrackerResolvesThinSupportForEachMomentumAndReverseTime) {
    auto bunch = makeBunch(0);
    OpalBeamline beamline;
    // Both the entrance and exit are beyond the nominal drift midpoint. The body
    // has zero length: only field support can prevent this region being skipped.
    auto field = std::make_shared<FieldSupportOnlyComponent>("THIN_E", 0.003, 0.0034, 1.0e8);
    field->setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector_t<double, 3>(0.0), Quaternion()));
    field->fixPosition();
    beamline.addElement(field, *bunch);
    beamline.prepareSections();
    beamline.buildElementMapping(defaultElementMappingGrid, defaultElementReach);
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    for (const std::string method : {"BORIS", "RK4", "DOP853"}) {
        SCOPED_TRACE(method);
        ExternalFieldRayTracker tracker(
                beamline, reference, ExternalFieldRayTracker::parseIntegrationMethod(method));
        for (double momentum : {0.9, 1.0, 1.1}) {
            ExternalFieldRayTracker::State initial;
            initial.momentum(2) = momentum;
            std::vector<ExternalFieldRayTracker::Step> steps;
            const auto final = tracker.advance(initial, 2.0e-11, &steps);
            {
                std::vector<ExternalFieldRayTracker::Step> measuredSteps;
                const auto measured = tracker.advance(initial, 2.0e-11, &measuredSteps);
                EXPECT_EQ(measuredSteps.size(), steps.size());
                for (unsigned i = 0; i < 3; ++i) {
                    EXPECT_EQ(measured.position(i), final.position(i));
                    EXPECT_EQ(measured.momentum(i), final.momentum(i));
                    EXPECT_EQ(measured.positionCorrection(i), final.positionCorrection(i));
                }
                EXPECT_EQ(measured.time, final.time);
                EXPECT_EQ(measured.pathLength, final.pathLength);
                EXPECT_EQ(measured.timeCorrection, final.timeCorrection);
                EXPECT_EQ(measured.pathLengthCorrection, final.pathLengthCorrection);
            }
            ASSERT_GT(final.position(2), 0.0034);
            ASSERT_GT(steps.size(), 1);
            // Exact work-energy identity: gamma_out - gamma_in = q E L / (m c^2).
            // 1e-10 in gamma is much smaller than the 4.26e-5 signal; it allows the
            // second-order orbit error and floating-point boundary localization.
            const double expectedGamma =
                    std::sqrt(1.0 + momentum * momentum) + 1.0e8 * 0.0004 / reference.getM();
            EXPECT_NEAR(
                    std::sqrt(1.0 + dot(final.momentum, final.momentum)), expectedGamma, 1.0e-10);
            EXPECT_NEAR(final.time, 2.0e-11, 1.0e-24);
            const auto recovered                     = tracker.advance(final, -2.0e-11);
            const Vector_t<double, 3> displacement   = recovered.position - initial.position;
            const Vector_t<double, 3> momentumChange = recovered.momentum - initial.momentum;
            EXPECT_LT(euclidean_norm(displacement), 1.0e-10);
            EXPECT_LT(euclidean_norm(momentumChange), 1.0e-9);
        }
    }
}

TEST_F(ExternalFieldRayTrackerTest, MembershipReusePreservesUncachedSubstepsExactly) {
    auto bunch = makeBunch(0);
    OpalBeamline beamline;
    auto electric = std::make_shared<FieldSupportOnlyComponent>("E", 0.003, 0.0034, 1.0e8);
    auto magnetic = std::make_shared<FieldSupportOnlyComponent>("B", 0.0032, 0.004, 0.0, 0.02);
    for (const auto& field : {electric, magnetic}) {
        field->setCSTrafoGlobal2Local(
                CoordinateSystemTrafo(Vector_t<double, 3>(0.0), Quaternion()));
        field->fixPosition();
        beamline.addElement(field, *bunch);
    }
    beamline.prepareSections();
    beamline.buildElementMapping(defaultElementMappingGrid, defaultElementReach);
    PartData particle(1.0, 9.382720813e8, 1.0e6);
    for (const std::string method : {"BORIS", "RK4", "DOP853"}) {
        SCOPED_TRACE(method);
        ExternalFieldRayTracker tracker(
                beamline, particle, ExternalFieldRayTracker::parseIntegrationMethod(method));
        const std::array<std::pair<double, double>, 3> cases{{
                {0.0, 4.e-11},
                {0.006, -4.e-11},
                {0.003, 4.e-12},
        }};
        for (const auto& [z, dt] : cases) {
            SCOPED_TRACE(z);
            for (const double momentum : {0.9, 1.1}) {
                ExternalFieldRayTracker::State initial, expected, actual;
                initial.position(2) = z;
                initial.momentum(2) = momentum;
                // Nonzero compensated low parts are passed through all recursive states.
                initial.positionCorrection(2) = 1.e-20;
                initial.timeCorrection        = 1.e-26;
                std::vector<ExternalFieldRayTracker::Step> oldSteps, newSteps;
                expected = advanceWithoutMembershipReuse(tracker, beamline, initial, dt, oldSteps);
                actual   = tracker.advance(initial, dt, &newSteps);
                expectSameRayBits(expected, actual);
                ASSERT_EQ(oldSteps.size(), newSteps.size());
                for (std::size_t i = 0; i < oldSteps.size(); ++i) {
                    expectSameRayBits(oldSteps[i].midpoint, newSteps[i].midpoint);
                    expectSameRayBits(oldSteps[i].end, newSteps[i].end);
                    EXPECT_EQ(oldSteps[i].duration, newSteps[i].duration);
                    EXPECT_EQ(oldSteps[i].hitMaterial, newSteps[i].hitMaterial);
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_EQ(oldSteps[i].electric(d), newSteps[i].electric(d));
                        EXPECT_EQ(oldSteps[i].magnetic(d), newSteps[i].magnetic(d));
                    }
                }
            }
        }
    }
}

TEST_F(ExternalFieldRayTrackerTest, RayTrackerSumsOverlappingSupportsAcrossBoundaries) {
    auto bunch = makeBunch(0);
    OpalBeamline beamline;
    auto first  = std::make_shared<FieldSupportOnlyComponent>("E1", 0.001, 0.003, 1.0e8);
    auto second = std::make_shared<FieldSupportOnlyComponent>("E2", 0.002, 0.0035, -2.0e7);
    for (const auto& field : {first, second}) {
        field->setCSTrafoGlobal2Local(
                CoordinateSystemTrafo(Vector_t<double, 3>(0.0), Quaternion()));
        field->fixPosition();
        beamline.addElement(field, *bunch);
    }
    beamline.prepareSections();
    beamline.buildElementMapping(defaultElementMappingGrid, defaultElementReach);
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    for (const std::string method : {"BORIS", "RK4", "DOP853"}) {
        SCOPED_TRACE(method);
        ExternalFieldRayTracker tracker(
                beamline, reference, ExternalFieldRayTracker::parseIntegrationMethod(method));
        ExternalFieldRayTracker::State initial;
        initial.momentum(2) = 1.0;
        const auto final    = tracker.advance(initial, 2.0e-11);
        ASSERT_GT(final.position(2), 0.0035);
        const double work = 1.0e8 * 0.002 - 2.0e7 * 0.0015;
        EXPECT_NEAR(
                std::sqrt(1.0 + dot(final.momentum, final.momentum)),
                std::sqrt(2.0) + work / reference.getM(), 1.0e-10);
    }
}

TEST_F(ExternalFieldRayTrackerTest, RungeKuttaMagneticHelixHasExpectedGlobalOrder) {
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    constexpr double magnetic = 10.0, transverse = 0.2, longitudinal = 0.1;
    const double gamma = std::sqrt(1.0 + transverse * transverse + longitudinal * longitudinal);
    const double omega =
            reference.getQ() * Physics::c * Physics::c * magnetic / (reference.getM() * gamma);
    const double duration = 4.0 / omega;
    for (const std::string method : {"RK4", "DOP853"}) {
        SCOPED_TRACE(method);
        ExternalFieldRayTracker tracker(
                beamline, reference, ExternalFieldRayTracker::parseIntegrationMethod(method));
        std::vector<double> errors;
        for (unsigned n : {4u, 8u, 16u}) {
            ExternalFieldRayTracker::State ray;
            ray.momentum = Vector_t<double, 3>(transverse, 0.0, longitudinal);
            for (unsigned i = 0; i < n; ++i)
                ray = tracker.step(ray, duration / n,
                                   [](const auto&, auto&, auto& b) {
                                       b(2) = magnetic;
                                       return false;
                                   })
                              .end;
            const Vector_t<double, 3> exactP(
                    transverse * std::cos(4.0), -transverse * std::sin(4.0), longitudinal);
            const Vector_t<double, 3> exactR(
                    Physics::c * transverse / (gamma * omega) * std::sin(4.0),
                    Physics::c * transverse / (gamma * omega) * (std::cos(4.0) - 1.0),
                    Physics::c * longitudinal / gamma * duration);
            const Vector_t<double, 3> dp = ray.momentum - exactP;
            const Vector_t<double, 3> dr = ray.position - exactR;
            const double exactPath =
                    Physics::c * std::hypot(transverse, longitudinal) / gamma * duration;
            errors.push_back(
                    std::max(
                            {euclidean_norm(dp), euclidean_norm(dr),
                             std::abs(ray.pathLength - exactPath)}));
            std::cout << std::scientific << std::setprecision(12) << "HELIX " << method
                      << " steps=" << n << " dt_s=" << duration / n
                      << " position_error_m=" << euclidean_norm(dr)
                      << " momentum_error=" << euclidean_norm(dp)
                      << " path_error_m=" << std::abs(ray.pathLength - exactPath) << " gamma_error="
                      << std::abs(std::sqrt(1.0 + dot(ray.momentum, ray.momentum)) - gamma)
                      << std::defaultfloat << '\n';
        }
        const double order = std::log2(errors[1] / errors[2]);
        EXPECT_GT(order, method == "RK4" ? 3.7 : 7.4);
        EXPECT_LT(order, method == "RK4" ? 4.4 : 8.6);
    }
}

TEST_F(ExternalFieldRayTrackerTest,
       RungeKuttaTimeDependentElectricForceUsesStageClocksAndChargeUnits) {
    OpalBeamline beamline;
    constexpr double duration = 1.e-9, amplitude = 3.e8;
    for (double charge : {-1.0, 1.0}) {
        PartData reference(charge, 9.382720813e8, 1.0e6);
        for (const std::string method : {"RK4", "DOP853"}) {
            SCOPED_TRACE(method);
            ExternalFieldRayTracker tracker(
                    beamline, reference, ExternalFieldRayTracker::parseIntegrationMethod(method));
            ExternalFieldRayTracker::State ray;
            ray.momentum(2)   = 1.0;
            const auto fields = [](const auto& state, auto& e, auto&) {
                e(2) = amplitude * std::cos(state.time / duration);
                return false;
            };
            for (unsigned i = 0; i < 20; ++i)
                ray = tracker.step(ray, duration / 20, fields).end;
            const double exact =
                    1.0
                    + charge * Physics::c * amplitude / reference.getM() * duration * std::sin(1.0);
            EXPECT_NEAR(ray.momentum(2), exact, method == "RK4" ? 1.e-9 : 3.e-15);
            for (unsigned i = 0; i < 20; ++i)
                ray = tracker.step(ray, -duration / 20, fields).end;
            EXPECT_NEAR(ray.momentum(2), 1.0, 3.e-14);
            EXPECT_NEAR(ray.position(2), 0.0, 1.e-10);
            EXPECT_NEAR(ray.pathLength, 0.0, 1.e-10);
        }
    }
}

TEST_F(ExternalFieldRayTrackerTest, RungeKuttaZeroStepsMaterialHitsAndFieldInitialization) {
    OpalBeamline beamline;
    PartData reference(1.0, 9.382720813e8, 1.0e6);
    for (const std::string method : {"RK4", "DOP853"}) {
        SCOPED_TRACE(method);
        ExternalFieldRayTracker tracker(
                beamline, reference, ExternalFieldRayTracker::parseIntegrationMethod(method));
        ExternalFieldRayTracker::State ray;
        ray.momentum(2)  = 1.0;
        unsigned calls   = 0;
        const auto field = [&](const auto&, auto& e, auto& b) {
            ++calls;
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_EQ(e(d), 0.0);
                EXPECT_EQ(b(d), 0.0);
            }
            e(2) = 1.0;
            return false;
        };
        tracker.step(ray, 0.0, field);
        EXPECT_EQ(calls, 0u);
        tracker.step(ray, 1.e-11, field);
        EXPECT_EQ(calls, method == "RK4" ? 9u : 17u);
        EXPECT_TRUE(tracker.step(ray, 1.e-11,
                                 [](const auto&, auto&, auto&) {
                                     return true;
                                 })
                            .hitMaterial);
        // Every sampled stage, including dense-output stages and midpoint,
        // must propagate a material hit and stop further field evaluations.
        const unsigned samples = method == "RK4" ? 9u : 17u;
        for (unsigned hit = 1; hit <= samples; ++hit) {
            calls             = 0;
            const auto result = tracker.step(ray, 1.e-11, [&](const auto&, auto&, auto&) {
                return ++calls == hit;
            });
            EXPECT_TRUE(result.hitMaterial);
            EXPECT_EQ(calls, hit);
        }
        EXPECT_THROW(
                tracker.step(ray, std::numeric_limits<double>::quiet_NaN(), field), OpalException);
    }
}

TEST_F(ExternalFieldRayTrackerTest, Dop853DenseMidpointAccuracyAndWork) {
    OpalBeamline beamline;
    PartData reference(1., 9.382720813e8, 1.e6);
    ExternalFieldRayTracker tracker(
            beamline, reference, ExternalFieldRayTracker::IntegrationMethod::DOP853);
    constexpr double b = 10., transverse = .2, longitudinal = .1;
    const double gamma = std::sqrt(1. + transverse * transverse + longitudinal * longitudinal);
    const double omega =
            reference.getQ() * Physics::c * Physics::c * b / (reference.getM() * gamma);
    std::vector<double> errors;
    for (double phase : {.8, .4}) {
        unsigned evaluations = 0;
        ExternalFieldRayTracker::State initial;
        initial.momentum = {transverse, 0., longitudinal};
        const auto step =
                tracker.step(initial, phase / omega, [&](const auto&, auto&, auto& magnetic) {
                    ++evaluations;
                    magnetic(2) = b;
                    return false;
                });
        EXPECT_EQ(evaluations, 17u);
        EXPECT_FALSE(step.hitMaterial);
        const double theta = .5 * phase, t = theta / omega;
        const Vector_t<double, 3> exactP(
                transverse * std::cos(theta), -transverse * std::sin(theta), longitudinal);
        const Vector_t<double, 3> exactR(
                Physics::c * transverse / (gamma * omega) * std::sin(theta),
                Physics::c * transverse / (gamma * omega) * (std::cos(theta) - 1),
                Physics::c * longitudinal / gamma * t);
        const double path  = Physics::c * std::hypot(transverse, longitudinal) / gamma * t;
        const double error = std::max(
                {euclidean_norm(Vector_t<double, 3>(step.midpoint.momentum - exactP)),
                 euclidean_norm(Vector_t<double, 3>(step.midpoint.position - exactR)),
                 std::abs(step.midpoint.pathLength - path)});
        errors.push_back(error);
        std::cout << "DENSE_MIDPOINT phase=" << phase << " error=" << std::setprecision(17) << error
                  << "\n";
        EXPECT_NEAR(step.midpoint.time, t, 1e-23);
        EXPECT_LT(error, 1e-6);
    }
    const double order = std::log2(errors[0] / errors[1]);
    std::cout << "DENSE_MIDPOINT local order=" << order << "\n";
    EXPECT_GT(order, 7.3);
    EXPECT_LT(order, 8.7);
}
