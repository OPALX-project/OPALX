/** @file TestSpaceChargeTransport.cpp @brief C++ fixtures for actual NS/Strang RMS transport. */
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>

#include "../TemporaryDirectory.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/PartData.h"
#include "Algorithms/SpaceChargeTransport.h"
#include "BeamlineGeometry/Geometry.h"
#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"
#include "Ippl.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"
#include "Structure/Beam.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

namespace {
    using namespace opalx::maps;
    using Matrix              = SpaceChargeMatrix;
    using Chart               = SpaceChargeChart;
    using Profile             = opalx::spacecharge::analytic3d::Profile;
    constexpr double momentum = 2, mass = 9.382720813e8, dt = 2e-11;
    const double gamma = std::sqrt(1 + momentum * momentum), beta = momentum / gamma;
    constexpr Profile profiles[]{Profile::UniformEllipsoid, Profile::GaussianRms};

    Matrix identity() {
        Matrix result{};
        for (unsigned i = 0; i < 6; ++i)
            result[i][i] = 1;
        return result;
    }
    Matrix multiply(const Matrix& a, const Matrix& b) {
        Matrix result{};
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                for (unsigned k = 0; k < 6; ++k)
                    result[i][j] += a[i][k] * b[k][j];
        return result;
    }
    Matrix transpose(const Matrix& a) {
        Matrix result{};
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                result[i][j] = a[j][i];
        return result;
    }
    Matrix derivative(const tps::Map<double>& map) {
        const auto j = map.jacobian(std::array<double, 6>{});
        Matrix result{};
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned k = 0; k < 6; ++k)
                result[i][k] = j(i, k);
        return result;
    }
    double difference(const Matrix& a, const Matrix& b, const Matrix& scale) {
        double sum = 0;
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                sum = std::hypot(sum, (a[i][j] - b[i][j]) / std::sqrt(scale[i][i] * scale[j][j]));
        return sum;
    }
    Matrix sigma() {
        const double widths[]{.002, .0002, .003, .0003, .004, .0004};
        Matrix l{};
        for (unsigned i = 0; i < 6; ++i)
            l[i][i] = widths[i];
        l[1][0] = .0001;
        l[2][0] = .0007;
        l[3][2] = -.0001;
        l[5][4] = .0001;
        return multiply(l, transpose(l));
    }
    SpaceChargeTransportSettings settings(std::size_t ns) {
        SpaceChargeTransportSettings s;
        s.ns                         = ns;
        s.external.richardsonLevels  = 1;
        s.external.integrationMethod = ExternalFieldRayTracker::IntegrationMethod::DOP853;
        return s;
    }
    LinearTransferMapReference station(double z, double path) {
        LinearTransferMapReference r;
        r.position   = {0, 0, z};
        r.momentum   = {0, 0, momentum};
        r.pathLength = path;
        r.time       = path / (beta * Physics::c);
        return r;
    }
    LinearTransferMapBuilder::Result sourcePass(
            OpalBeamline& beamline, PartData& species, double length = .25, double startZ = 0,
            unsigned samples = 1) {
        std::vector<LinearTransferMapBuilder::ReferenceSample> points;
        for (unsigned i = 0; i <= samples; ++i) {
            const double path = double(i) / samples * length;
            points.push_back({station(startZ + path, path)});
        }
        const auto s = settings(1);
        LinearTransferMapBuilder builder(beamline, species, dt, s.external);
        return builder.build(std::move(points), 0);
    }
    Matrix drift(double length) {
        auto result  = identity();
        result[0][1] = result[2][3] = length / momentum;
        result[4][5]                = length / (momentum * gamma * gamma);
        return result;
    }

    // Independently integrate dSigma/ds=A(Sigma)*Sigma+Sigma*A(Sigma)^T by fine RK4.
    // Shares the separately validated field tensor, but no production maps, kicks,
    // covariance transport or Strang schedule. All dense arithmetic here is scalar.
    Matrix momentRhs(
            const Matrix& covariance, Profile profile, double charge, double path,
            double externalStrength) {
        const double d[]{1, 1, gamma};
        opalx::spacecharge::analytic3d::Matrix rest{};
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j)
                rest[i][j] = d[i] * covariance[2 * i][2 * j] * d[j];
        const auto tensor =
                opalx::spacecharge::analytic3d::evaluateRestFrameTensor(profile, rest, charge);
        Matrix a{};
        a[0][1] = a[2][3] = 1 / momentum;
        a[4][5]           = 1 / (momentum * gamma * gamma);
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j)
                a[2 * i + 1][2 * j] = d[i] * tensor.gradient[i][j] * d[j] / (gamma * beta * mass);
        const double focusing = externalStrength * (1 + 8 * (.02 + path)) / (beta * mass);
        a[1][0] -= focusing;
        a[3][2] += focusing;
        const auto left = multiply(a, covariance), right = multiply(covariance, transpose(a));
        Matrix result{};
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                result[i][j] = left[i][j] + right[i][j];
        return result;
    }
    Matrix shifted(const Matrix& a, const Matrix& k, double h) {
        Matrix result{};
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                result[i][j] = a[i][j] + h * k[i][j];
        return result;
    }
    Matrix referenceMoments(
            Matrix covariance, Profile profile, double charge, unsigned count,
            double externalStrength = 0) {
        const double h = .25 / count;
        for (unsigned step = 0; step < count; ++step) {
            const double path = step * h;
            const auto k1     = momentRhs(covariance, profile, charge, path, externalStrength);
            const auto k2     = momentRhs(
                    shifted(covariance, k1, h / 2), profile, charge, path + h / 2,
                    externalStrength);
            const auto k3 = momentRhs(
                    shifted(covariance, k2, h / 2), profile, charge, path + h / 2,
                    externalStrength);
            const auto k4 = momentRhs(
                    shifted(covariance, k3, h), profile, charge, path + h, externalStrength);
            for (unsigned i = 0; i < 6; ++i)
                for (unsigned j = 0; j < 6; ++j)
                    covariance[i][j] += h / 6 * (k1[i][j] + 2 * k2[i][j] + 2 * k3[i][j] + k4[i][j]);
        }
        return covariance;
    }

    class DummyLine : public Beamline {
    public:
        DummyLine() : Beamline("SC_TEST") {}
        ElementType getType() const override { return ElementType::BEAMLINE; }
        Geometry& getGeometry() override { return geometry; }
        const Geometry& getGeometry() const override { return geometry; }
        void accept(BeamlineVisitor& v) const override { v.visitBeamline(*this); }
        ElementBase* clone() const override { return new DummyLine(*this); }
        void iterate(BeamlineVisitor&, bool) const override {}
        Geometry geometry{Geometry::makeNull()};
    };
    // Vacuum electrostatic quadrupole potential phi=k(z)*(x^2-y^2)/2, k linear in z.
    // Includes Ez=-k'(z)*(x^2-y^2)/2, so the deliberately varying field is conservative.
    class VaryingQuadrupole : public ElementBase {
    public:
        VaryingQuadrupole(std::string name, double length, double gradient = 1e7)
            : ElementBase(name),
              geometry(Geometry::makeStraight(length)),
              gradient(gradient),
              calls(std::make_shared<unsigned>(0)) {}
        void accept(BeamlineVisitor&) const override {}
        ElementBase* clone() const override { return new VaryingQuadrupole(*this); }
        void initialise(PartBunch_t*) override {}
        void finalise() override {}
        void apply(const std::shared_ptr<ParticleContainer_t>&) override {}
        void apply(
                const Vector_t<double, 3>&, const Vector_t<double, 3>&, const double&,
                Vector_t<double, 3>&, Vector_t<double, 3>&) override {}
        bool applyToReferenceParticle(
                const Vector_t<double, 3>& r, const Vector_t<double, 3>&, const double&,
                Vector_t<double, 3>& e, Vector_t<double, 3>&) override {
            ++*calls;
            const double k = gradient * (1 + 8 * r(2));
            e(0)           = -k * r(0);
            e(1)           = k * r(1);
            e(2)           = -4 * gradient * (r(0) * r(0) - r(1) * r(1));
            return false;
        }
        void getFieldExtent(double& start, double& end) const override {
            start = 0;
            end   = extent;
        }
        ElementType getType() const override { return type; }
        Geometry& getGeometry() override { return geometry; }
        const Geometry& getGeometry() const override { return geometry; }
        Geometry geometry;
        double gradient, extent = .5;
        ElementType type = ElementType::ANY;
        std::shared_ptr<unsigned> calls;
    };
    class StrangTransportTest : public ::testing::Test {
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
            oldOrder = tps::Tps<double>::getGlobalTruncOrder();
            oldInput = OpalData::getInstance()->getInputFn();
            directory.create("opalx-strang", ippl::Comm->getCommunicator());
            OpalData::getInstance()->storeInputFn((directory.path() / "strang.opal").string());
        }
        void TearDown() override {
            tps::Tps<double>::setGlobalTruncOrder(oldOrder);
            OpalData::getInstance()->storeInputFn(oldInput);
            ippl::Comm->barrier();
            directory.reset();
        }
        void add(OpalBeamline& lattice, VaryingQuadrupole& field) {
            Beam beam;
            PartBunch_t bunch(
                    {1}, {1}, {Beam::find("UNNAMED_BEAM")}, {0}, 1, "LF2",
                    opalx::spacecharge::CartesianDomainConfig3D{.periodicParticleBoundary = true});
            DummyLine line;
            DefaultVisitor visitor(line, false, false);
            field.setCSTrafoGlobal2Local(
                    CoordinateSystemTrafo(Vector_t<double, 3>(0), Quaternion()));
            field.fixPosition();
            lattice.visit(field, visitor, bunch);
            lattice.prepareSections();
        }
        int oldOrder;
        std::string oldInput;
        opalx::test::TemporaryDirectory directory;
    };
}  // namespace

TEST_F(StrangTransportTest, NsControlsActualHalfMapsAndPhysicalPathCoverage) {
    OpalBeamline lattice;
    PartData species(1, mass, momentum * mass);
    const auto source = sourcePass(lattice, species);
    ASSERT_TRUE(source.pass);
    // advanceToPathLength stops at max(1e-12*DT,64*epsilon*1m/c).
    // Allow twice that spatial floor for station arithmetic and twice again for
    // a difference of two independently localized endpoints. No tracker tolerance changes.
    const double pathTolerance =
            2 * std::max(1e-12 * dt * Physics::c, 64 * std::numeric_limits<double>::epsilon());
    for (std::size_t ns : {1u, 2u, 4u, 8u, 64u}) {
        const auto plan = buildSpaceChargeTransportPlan(
                lattice, species, *source.pass, dt, settings(ns), "coverage");
        ASSERT_EQ(plan.slices().size(), ns);
        EXPECT_EQ(plan.ns(), ns);
        EXPECT_NEAR(plan.pathLength(), .25, 1e-15);
        double length = 0;
        for (std::size_t i = 0; i < ns; ++i) {
            const auto& slice = plan.slices()[i];
            EXPECT_EQ(slice.sourceSegment, 0u);
            EXPECT_EQ(slice.subdivision, i);
            EXPECT_NEAR(slice.pathMidpoint, (i + .5) * .25 / ns, pathTolerance);
            EXPECT_NEAR(slice.integrationLength, .25 / ns, 2 * pathTolerance);
            EXPECT_NEAR(slice.leftExternal.exit.position(2), slice.pathMidpoint, 1e-13);
            EXPECT_NEAR(slice.rightExternal.exit.position(2), slice.pathEnd, 1e-13);
            EXPECT_EQ(slice.left.context().output(), slice.right.context().input());
            EXPECT_LT(difference(derivative(slice.left), drift(.125 / ns), identity()), 2e-10);
            EXPECT_LT(difference(derivative(slice.right), drift(.125 / ns), identity()), 2e-10);
            length += slice.integrationLength;
        }
        EXPECT_NEAR(length, plan.pathLength(), 1e-15);
    }
}

TEST_F(StrangTransportTest, VaryingExternalFieldHasDistinctActualHalves) {
    OpalBeamline lattice;
    PartData species(1, mass, momentum * mass);
    VaryingQuadrupole field("VARYING", .5);
    add(lattice, field);
    const auto source = sourcePass(lattice, species, .25, .02);
    ASSERT_TRUE(source.pass);
    *field.calls    = 0;
    const auto plan = buildSpaceChargeTransportPlan(
            lattice, species, *source.pass, dt, settings(1), "varying");
    ASSERT_EQ(plan.slices().size(), 1u);
    const auto& slice = plan.slices()[0];
    EXPECT_GT(*field.calls, 100u);
    EXPECT_GT(std::abs(slice.leftExternal.matrix(1, 0) - slice.rightExternal.matrix(1, 0)), 1e-5);
    EXPECT_NEAR(slice.leftExternal.exit.position(2), .145, 1e-13);
    const auto result = transportSpaceCharge(plan, sigma(), Profile::GaussianRms, 0);
    // With zero charge the composed actual halves reproduce independently built full transport.
    const auto& whole    = source.pass->polynomial();
    const auto canonical = tps::compose(
            plan.exit().map(Chart::SectionSlopes, Chart::CommonTimeCanonical, "out"),
            tps::compose(
                    whole,
                    plan.entrance().map(Chart::CommonTimeCanonical, Chart::SectionSlopes, "in")));
    EXPECT_LT(difference(derivative(result.polynomial), derivative(canonical), identity()), 2e-8);
}

TEST_F(StrangTransportTest, MidpointDensityIsTransportedAndRecomputedAcrossSlicesAndRuns) {
    OpalBeamline lattice;
    PartData species(1, mass, momentum * mass);
    const auto source = sourcePass(lattice, species);
    const auto plan   = buildSpaceChargeTransportPlan(
            lattice, species, *source.pass, dt, settings(4), "density");
    const auto input       = sigma();
    const auto result      = transportSpaceCharge(plan, input, Profile::GaussianRms, 1e-9);
    const auto firstHalf   = drift(.25 / 8);
    const auto expectedMid = multiply(multiply(firstHalf, input), transpose(firstHalf));
    EXPECT_LT(difference(result.steps[0].midpointCovariance, expectedMid, input), 2e-10);
    EXPECT_GT(std::abs(result.steps[0].midpointCovariance[0][0] - input[0][0]), 1e-9);
    EXPECT_GT(
            std::abs(
                    result.steps.back().kick.restSpatialCovariance[0][0]
                    - result.steps.front().kick.restSpatialCovariance[0][0]),
            1e-8);
    auto wider = input;
    for (auto& row : wider)
        for (auto& entry : row)
            entry *= 4;
    const auto repeated = transportSpaceCharge(plan, wider, Profile::GaussianRms, 1e-9);
    EXPECT_NEAR(
            repeated.steps[0].kick.forceGradient[0][0],
            result.steps[0].kick.forceGradient[0][0] / 8,
            2e-10 * std::abs(result.steps[0].kick.forceGradient[0][0]));
    const auto j = derivative(result.polynomial);
    EXPECT_LT(
            difference(result.covariance, multiply(multiply(j, input), transpose(j)), input),
            2e-11);
    // A kick changes momenta/correlations, but leaves its full spatial block unchanged.
    const auto afterKick = transportSpaceChargeCovariance(
            result.steps[0].kick.polynomial, result.steps[0].midpointCovariance);
    for (unsigned i : {0u, 2u, 4u})
        for (unsigned k : {0u, 2u, 4u})
            EXPECT_NEAR(afterKick[i][k], result.steps[0].midpointCovariance[i][k], 1e-18);
}

TEST_F(StrangTransportTest, BothProfilesConvergeAtSecondOrderAgainstIndependentMomentOde) {
    OpalBeamline lattice;
    PartData species(1, mass, momentum * mass);
    const auto source = sourcePass(lattice, species);
    const auto input  = sigma();
    for (auto profile : profiles) {
        const auto reference = referenceMoments(input, profile, 1e-9, 1024);
        const auto refined   = referenceMoments(input, profile, 1e-9, 2048);
        EXPECT_LT(difference(reference, refined, input), 1e-10);
        std::vector<double> errors;
        for (std::size_t ns : {1u, 2u, 4u, 8u}) {
            const auto plan = buildSpaceChargeTransportPlan(
                    lattice, species, *source.pass, dt, settings(ns), "convergence");
            const auto result  = transportSpaceCharge(plan, input, profile, 1e-9);
            const double error = difference(result.covariance, refined, input);
            errors.push_back(error);
            std::cout << "STRANG profile=" << int(profile) << " NS=" << ns
                      << " scaled_error=" << std::setprecision(17) << error << '\n';
        }
        for (unsigned i = 1; i < errors.size(); ++i) {
            const double ratio = errors[i - 1] / errors[i];
            EXPECT_GT(ratio, 3.7);
            EXPECT_LT(ratio, 4.3);
        }
    }
}

TEST_F(StrangTransportTest, VaryingExternalHalvesAlsoGiveSecondOrderSelfConsistentTransport) {
    OpalBeamline lattice;
    PartData species(1, mass, momentum * mass);
    VaryingQuadrupole field("VARYING_CONVERGENCE", .5);
    add(lattice, field);
    const auto source = sourcePass(lattice, species, .25, .02);
    const auto input  = sigma();
    for (auto profile : profiles) {
        const auto reference = referenceMoments(input, profile, 1e-9, 1024, 1e7);
        const auto refined   = referenceMoments(input, profile, 1e-9, 2048, 1e7);
        EXPECT_LT(difference(reference, refined, input), 1e-10);
        double previous = 0;
        for (std::size_t ns : {1u, 2u, 4u, 8u}) {
            const auto plan = buildSpaceChargeTransportPlan(
                    lattice, species, *source.pass, dt, settings(ns), "varying-convergence");
            const auto result  = transportSpaceCharge(plan, input, profile, 1e-9);
            const double error = difference(result.covariance, refined, input);
            std::cout << "STRANG_VARYING profile=" << int(profile) << " NS=" << ns
                      << " scaled_error=" << std::setprecision(17) << error << '\n';
            if (previous > 0) {
                EXPECT_GT(previous / error, 3.7);
                EXPECT_LT(previous / error, 4.3);
            }
            previous = error;
        }
    }
}

TEST_F(StrangTransportTest, ZeroCurrentPreservesCanonicalFormAndAmbientOrder) {
    OpalBeamline lattice;
    PartData species(1, mass, momentum * mass);
    const auto source = sourcePass(lattice, species);
    const auto plan =
            buildSpaceChargeTransportPlan(lattice, species, *source.pass, dt, settings(2), "zero");
    Matrix form{};
    for (unsigned i = 0; i < 6; i += 2) {
        form[i][i + 1] = 1;
        form[i + 1][i] = -1;
    }
    for (int order : {0, 1, 4, tps::Tps<double>::EXACT}) {
        tps::Tps<double>::setGlobalTruncOrder(order);
        const auto result = transportSpaceCharge(plan, sigma(), Profile::UniformEllipsoid, 0);
        const auto j      = derivative(result.polynomial);
        EXPECT_LT(difference(j, drift(.25), identity()), 2e-10);
        EXPECT_LT(difference(multiply(multiply(transpose(j), form), j), form, identity()), 2e-10);
        EXPECT_EQ(tps::Tps<double>::getGlobalTruncOrder(), order);
        EXPECT_EQ(
                result.polynomial.context().input(),
                plan.entrance().chart(Chart::CommonTimeCanonical));
        EXPECT_EQ(
                result.polynomial.context().output(),
                plan.exit().chart(Chart::CommonTimeCanonical));
        EXPECT_EQ(result.polynomial.context().segments().size(), 7 * plan.slices().size());
        for (const auto& step : result.steps)
            EXPECT_EQ(step.kick.restTensor.integrandEvaluations, 0u);
    }
}

TEST_F(StrangTransportTest, IncludesBodyAndOwnerlessDriftIntervalsExactlyOnce) {
    OpalBeamline lattice;
    PartData species(1, mass, momentum * mass);
    VaryingQuadrupole field("BODY", .1, 0);
    field.extent = .1;
    add(lattice, field);
    const auto source = sourcePass(lattice, species, .25, -.05, 20);
    ASSERT_TRUE(source.pass);
    ASSERT_EQ(source.segments.size(), 3u);
    const auto plan = buildSpaceChargeTransportPlan(
            lattice, species, *source.pass, dt, settings(2), "drifts");
    EXPECT_EQ(plan.slices().size(), 6u);
    EXPECT_NEAR(plan.pathLength(), .25, 1e-13);
    EXPECT_TRUE(source.segments.front().owners.empty());
    EXPECT_TRUE(source.segments.back().owners.empty());
    EXPECT_EQ(source.segments[1].owners.size(), 1u);
    for (unsigned sourceIndex = 0; sourceIndex < 3; ++sourceIndex)
        for (unsigned n = 0; n < 2; ++n)
            EXPECT_EQ(plan.slices()[2 * sourceIndex + n].sourceSegment, sourceIndex);
}

TEST_F(StrangTransportTest, RejectsRecordedAndNewlySampledOverlapsWithoutPublishing) {
    OpalBeamline lattice;
    PartData species(1, mass, momentum * mass);
    VaryingQuadrupole first("FIRST", .5), second("SECOND", .5);
    add(lattice, first);
    add(lattice, second);
    const auto overlapping = sourcePass(lattice, species, .25, .02);
    ASSERT_TRUE(overlapping.pass);
    EXPECT_THROW(
            (void)buildSpaceChargeTransportPlan(
                    lattice, species, *overlapping.pass, dt, settings(1), "overlap"),
            OpalException);
    OpalBeamline changed;
    add(changed, first);
    const auto original = sourcePass(changed, species, .25, .02);
    ASSERT_FALSE(original.pass->segments()[0]->metadata().includesOverlappingFields);
    // Existing metadata cannot certify a changed lattice. Strict sampling must still reject it.
    add(changed, second);
    EXPECT_THROW(
            (void)buildSpaceChargeTransportPlan(
                    changed, species, *original.pass, dt, settings(1), "new-overlap"),
            OpalException);
    EXPECT_TRUE(original.pass->available());
    OpalBeamline rfLattice;
    const auto beforeRf = sourcePass(rfLattice, species);
    VaryingQuadrupole rf("UNSUPPORTED_RF", .5);
    rf.type = ElementType::TRAVELINGWAVE;
    add(rfLattice, rf);
    EXPECT_THROW(
            (void)buildSpaceChargeTransportPlan(
                    rfLattice, species, *beforeRf.pass, dt, settings(1), "rf"),
            OpalException);
}

TEST_F(StrangTransportTest, RejectsInvalidSettingsInputsAndNonlinearCovarianceTransport) {
    OpalBeamline lattice;
    PartData species(1, mass, momentum * mass);
    const auto source = sourcePass(lattice, species);
    for (unsigned mode = 0; mode < 5; ++mode) {
        auto invalid = settings(1);
        if (mode == 0) invalid.ns = 0;
        if (mode == 1) invalid.maximumSlices = 0;
        if (mode == 2) {
            invalid.ns            = 2;
            invalid.maximumSlices = 1;
        }
        if (mode == 3)
            invalid.external.integrationMethod = ExternalFieldRayTracker::IntegrationMethod::RK4;
        if (mode == 4) invalid.ns = std::numeric_limits<std::size_t>::max();
        EXPECT_THROW(
                (void)buildSpaceChargeTransportPlan(
                        lattice, species, *source.pass, dt, invalid, "invalid"),
                OpalException);
    }
    EXPECT_THROW(
            (void)buildSpaceChargeTransportPlan(
                    lattice, species, *source.pass, 0, settings(1), "invalid"),
            OpalException);
    EXPECT_THROW(
            (void)buildSpaceChargeTransportPlan(
                    lattice, species, *source.pass, dt, settings(1), ""),
            OpalException);
    const auto plan =
            buildSpaceChargeTransportPlan(lattice, species, *source.pass, dt, settings(1), "valid");
    auto invalid  = sigma();
    invalid[0][0] = 0;
    EXPECT_THROW(
            (void)transportSpaceCharge(plan, invalid, Profile::GaussianRms, 1e-9), OpalException);
    EXPECT_THROW(
            (void)transportSpaceCharge(
                    plan, sigma(), Profile::GaussianRms, std::numeric_limits<double>::infinity()),
            OpalException);
    EXPECT_THROW(
            (void)transportSpaceCharge(plan, sigma(), static_cast<Profile>(99), 0), OpalException);
    auto polynomial  = plan.slices()[0].left.polynomial();
    polynomial[0][0] = 1;
    tps::Map<double> shifted(plan.slices()[0].left.context(), polynomial);
    EXPECT_THROW((void)transportSpaceChargeCovariance(shifted, sigma()), OpalException);
    const auto constant = plan.slices()[0].left.truncated(0);
    EXPECT_THROW((void)transportSpaceChargeCovariance(constant, sigma()), OpalException);
    tps::VpsMap<double> nonlinear(6, 6);
    for (unsigned i = 0; i < 6; ++i) {
        nonlinear[i]        = tps::Tps<double>(2, 6);
        nonlinear[i][i + 1] = 1;
    }
    const auto x = tps::Tps<double>::makeVariable(6, 0);
    nonlinear[0] += x * x;
    const tps::Map<double> quadratic(
            tps::MapContext(
                    plan.slices()[0].left.context().input(),
                    plan.slices()[0].left.context().output(), 2),
            nonlinear);
    EXPECT_THROW((void)transportSpaceChargeCovariance(quadratic, sigma()), OpalException);
    const auto good = transportSpaceCharge(plan, sigma(), Profile::GaussianRms, 1e-9);
    EXPECT_EQ(good.steps.size(), 1u);
}

TEST_F(StrangTransportTest, MpiRanksAgreeOnFullCovarianceAndComposedMap) {
    OpalBeamline lattice;
    PartData species(1, mass, momentum * mass);
    const auto source = sourcePass(lattice, species);
    const auto plan =
            buildSpaceChargeTransportPlan(lattice, species, *source.pass, dt, settings(2), "mpi");
    const auto result = transportSpaceCharge(plan, sigma(), Profile::GaussianRms, 1e-9);
    const auto j      = derivative(result.polynomial);
    std::array<double, 72> value{}, minimum{}, maximum{};
    for (unsigned i = 0; i < 6; ++i)
        for (unsigned k = 0; k < 6; ++k) {
            value[6 * i + k]      = j[i][k];
            value[36 + 6 * i + k] = result.covariance[i][k];
        }
    MPI_Allreduce(value.data(), minimum.data(), 72, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(value.data(), maximum.data(), 72, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    EXPECT_EQ(minimum, maximum);
}
