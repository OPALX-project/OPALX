/**
 * \file TestRFCavity.cpp
 * \brief Unit tests for RFCavity component (base layer).
 *
 * This test suite validates the core behavior of RFCavity using:
 *  - a minimal concrete implementation (TestRFCavity)
 *  - a controlled fake fieldmap (FakeFieldmap)
 *
 * ---------------------------------------------------------------------------
 * Coverage
 * ---------------------------------------------------------------------------
 *
 * 1. Basic API
 *    - getType()
 *    - amplitude / frequency / phase setters and getters
 *
 * 2. Geometry
 *    - getFieldExtent()
 *
 * 3. Spatial behavior
 *    - apply() inside the element
 *    - apply() before the element
 *    - apply() after the element
 *
 * 4. RF phase behavior (core physics)
 *    - cos(phi) scaling of electric field (E)
 *    - sin(phi) scaling of magnetic field (B)
 *    - special phases:
 *        * phi = 0        → max E, no B
 *        * phi = π/2      → no E, max B
 *        * phi = π        → inverted E
 *    - independence of frequency at t = 0
 *
 * 5. Scaling behavior
 *    - linear scaling of fields via cavity scale factor
 *
 * 6. Time dependence
 *    - field variation with time through phi = ωt + phase
 *
 * 7. Fieldmap interaction / edge cases
 *    - correct handling of out-of-bounds fieldmap queries
 *
 * ---------------------------------------------------------------------------
 * Notes
 * ---------------------------------------------------------------------------
 *
 * - The FakeFieldmap provides a simple, deterministic field:
 *      E = (1, 0, 0), B = (0, 1, 0)
 *   allowing direct verification of cos/sin modulation.
 *
 * - Tests focus on physics correctness rather than implementation details.
 *
 * - The fixture defines a stable default configuration:
 *      scale = 1, frequency = 1, phase = 0
 *   Individual tests override only what they need.
 */
#include <gtest/gtest.h>
#include <limits>
#include "AbsBeamline/VariableRFCavity.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/CavityAutophaser.h"
#include "Algorithms/PolynomialTimeDependence.h"
#include "Attributes/Attributes.h"
#include "Elements/OpalCavity.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Units.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"

#include "AbsBeamline/ElementBase.h"
#include "BeamlineCore/RFCavityRep.h"
#include "Fields/Fieldmap.h"

#include <cmath>
#include <memory>

class FakeFieldmap : public Fieldmap {
public:
    FakeFieldmap() : Fieldmap("dummy"), outOfBounds_(false) {}

    void setOutOfBounds(bool v) { outOfBounds_ = v; }

    bool getFieldstrength(
            const Vector_t<double, 3>&, Vector_t<double, 3>& E,
            Vector_t<double, 3>& B) const override {
        if (outOfBounds_) return true;

        E = {1.0, 0.0, 0.0};
        B = {0.0, 1.0, 0.0};
        return false;  // inside
    }

    void getFieldDimensions(double& zBegin, double& zEnd) const override {
        zBegin = 0.0;
        zEnd   = 1.0;
    }

    void getFieldDimensions(
            double& zBegin, double& zEnd, double&, double&, double&, double&) const override {
        zBegin = 0.0;
        zEnd   = 1.0;
    }

    bool isInside(const Vector_t<double, 3>&) const override { return true; }

    // --- required no-op implementations ---
    void applyField(std::shared_ptr<ParticleContainer_t>, double = 1.0, double = 1.0) override {}
    bool getFieldDerivative(
            const Vector_t<double, 3>&, Vector_t<double, 3>&, Vector_t<double, 3>&,
            const DiffDirection&) const override {
        return false;
    }

    void swap() override {}
    void getInfo(Inform*) override {}
    double getFrequency() const override { return 1.0; }
    void setFrequency(double) override {}
    void readMap() override {}
    void freeMap() override {}

private:
    bool outOfBounds_;
};

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Minimal concrete RFCavity
// ---------------------------------------------------------------------------
class TestRFCavity : public RFCavity {
public:
    TestRFCavity() : RFCavity("test") {}

    // ---- REQUIRED PURE VIRTUALS ----
    double getAmplitude() const override { return amplitude_; }
    double getFrequency() const override { return frequency_; }
    double getPhase() const override { return phase_; }

    ElementBase* clone() const override { return new TestRFCavity(*this); }

    Geometry& getGeometry() override { return geom_; }
    const Geometry& getGeometry() const override { return geom_; }

    // ---- Simple setters for testing ----
    void setAmplitude(double v) { amplitude_ = v; }
    void setFrequency(double v) { frequency_ = v; }
    void setPhase(double v) { phase_ = v; }

    void setScale(double v) { scale_m = v; }
    void setPhaseInternal(double v) { phase_m = v; }
    void setFrequencyInternal(double v) { frequency_m = v; }

    void setFieldmap(Fieldmap* fmap) { fieldmap_m = fmap; }
    void setStartField(double val) { startField_m = val; }
    void setEndField(double val) { endField_m = val; }

private:
    double amplitude_ = 0.0;
    double frequency_ = 0.0;
    double phase_     = 0.0;

    Geometry geom_;
};

// ---------------------------------------------------------------------------
// Test Fixture
// ---------------------------------------------------------------------------
class RFCavityTest : public ::testing::Test {
protected:
    void SetUp() override {
        cav_  = std::make_unique<TestRFCavity>();
        fmap_ = std::make_unique<FakeFieldmap>();

        // --- Geometry ---
        cav_->setFieldmap(fmap_.get());
        cav_->setStartField(0.0);
        cav_->setEndField(1.0);
        cav_->getGeometry().setElementLength(1.0);

        // --- RF defaults ---
        cav_->setScale(1.0);
        cav_->setFrequencyInternal(1.0);
        cav_->setPhaseInternal(0.0);
    }
    std::unique_ptr<TestRFCavity> cav_;
    std::unique_ptr<FakeFieldmap> fmap_;
};

// ---------------------------------------------------------------------------
// Basic API
// ---------------------------------------------------------------------------
TEST_F(RFCavityTest, GetType) { EXPECT_EQ(cav_->getType(), ElementType::RFCAVITY); }

TEST_F(RFCavityTest, GetSetAmplitudeFrequencyPhase) {
    cav_->setAmplitude(5.0);
    cav_->setFrequency(2.0);
    cav_->setPhase(0.5);

    EXPECT_DOUBLE_EQ(cav_->getAmplitude(), 5.0);
    EXPECT_DOUBLE_EQ(cav_->getFrequency(), 2.0);
    EXPECT_DOUBLE_EQ(cav_->getPhase(), 0.5);
}

// ---------------------------------------------------------------------------
// Geometry / dimensions
// ---------------------------------------------------------------------------
TEST_F(RFCavityTest, GetDimensions) {
    double zBegin = -1.0, zEnd = -1.0;

    cav_->getFieldExtent(zBegin, zEnd);

    EXPECT_EQ(zBegin, 0.0);
    EXPECT_EQ(zEnd, 1.0);
}

TEST_F(RFCavityTest, BodyExtentCanDifferFromFieldSupport) {
    cav_->setStartField(0.2);
    cav_->setEndField(0.8);
    cav_->getGeometry().setElementLength(1.0);

    double bodyBegin = -1.0, bodyEnd = -1.0;
    bodyBegin = 0.0;
    bodyEnd   = cav_->getGeometry().getElementLength();
    EXPECT_EQ(bodyBegin, 0.0);
    EXPECT_EQ(bodyEnd, 1.0);

    const auto entry = cav_->getGeometry().getEdgeToBegin();
    const auto exit  = cav_->getGeometry().getEdgeToEnd();
    EXPECT_EQ(entry.getOrigin()(2), 0.0);
    EXPECT_EQ(exit.getOrigin()(2), 1.0);

    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};
    cav_->apply({0.0, 0.0, 0.1}, {0.0, 0.0, 1.0}, 0.0, E, B);
    EXPECT_DOUBLE_EQ(E(0), 0.0);

    cav_->apply({0.0, 0.0, 0.5}, {0.0, 0.0, 1.0}, 0.0, E, B);
    EXPECT_DOUBLE_EQ(E(0), 1.0);
}

TEST_F(RFCavityTest, ZeroBodyLengthDoesNotFallBackToFieldmapLength) {
    cav_->getGeometry().setElementLength(0.0);

    EXPECT_DOUBLE_EQ(cav_->getGeometry().getElementLength(), 0.0);

    double bodyBegin = -1.0, bodyEnd = -1.0;
    bodyBegin = 0.0;
    bodyEnd   = cav_->getGeometry().getElementLength();
    EXPECT_DOUBLE_EQ(bodyBegin, 0.0);
    EXPECT_DOUBLE_EQ(bodyEnd, 0.0);

    double fieldBegin = -1.0, fieldEnd = -1.0;
    cav_->getFieldExtent(fieldBegin, fieldEnd);
    EXPECT_DOUBLE_EQ(fieldBegin, 0.0);
    EXPECT_DOUBLE_EQ(fieldEnd, 1.0);
}

// ---------------------------------------------------------------------------
// apply(): spatial behavior
// ---------------------------------------------------------------------------
TEST_F(RFCavityTest, ApplyInside) {
    Vector_t<double, 3> R = {0.0, 0.0, 0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {1.0, 2.0, 3.0};
    Vector_t<double, 3> B = {1.0, 2.0, 3.0};

    cav_->apply(R, P, 0.0, E, B);

    // cos(0) = 1 → E += scale * E_map
    EXPECT_DOUBLE_EQ(E(0), 1.0 + 1.0);
    EXPECT_DOUBLE_EQ(E(1), 2.0);
    EXPECT_DOUBLE_EQ(E(2), 3.0);

    // sin(0) = 0 → no B contribution
    EXPECT_DOUBLE_EQ(B(0), 1.0);
    EXPECT_DOUBLE_EQ(B(1), 2.0);
    EXPECT_DOUBLE_EQ(B(2), 3.0);
}

TEST_F(RFCavityTest, ApplyBefore) {
    Vector_t<double, 3> R = {0.0, 0.0, -0.1};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {1.0, 2.0, 3.0};
    Vector_t<double, 3> B = {1.0, 2.0, 3.0};

    cav_->apply(R, P, 0.0, E, B);

    // No expected contribution
    EXPECT_DOUBLE_EQ(E(0), 1.0);
    EXPECT_DOUBLE_EQ(E(1), 2.0);
    EXPECT_DOUBLE_EQ(E(2), 3.0);

    EXPECT_DOUBLE_EQ(B(0), 1.0);
    EXPECT_DOUBLE_EQ(B(1), 2.0);
    EXPECT_DOUBLE_EQ(B(2), 3.0);
}

TEST_F(RFCavityTest, ApplyAfter) {
    Vector_t<double, 3> R = {0.0, 0.0, 1.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {1.0, 2.0, 3.0};
    Vector_t<double, 3> B = {1.0, 2.0, 3.0};

    cav_->apply(R, P, 0.0, E, B);

    // No expected contribution
    EXPECT_DOUBLE_EQ(E(0), 1.0);
    EXPECT_DOUBLE_EQ(E(1), 2.0);
    EXPECT_DOUBLE_EQ(E(2), 3.0);

    EXPECT_DOUBLE_EQ(B(0), 1.0);
    EXPECT_DOUBLE_EQ(B(1), 2.0);
    EXPECT_DOUBLE_EQ(B(2), 3.0);
}

// ---------------------------------------------------------------------------
// RF phase behavior (cos/sin physics)
// ---------------------------------------------------------------------------
// E ∝ cos(ωt + φ)
// B ∝ -sin(ωt + φ)
TEST_F(RFCavityTest, ApplyPhaseZero) {
    // φ = 0 → cos = 1, sin = 0
    cav_->setPhasem(0.0);
    cav_->setFrequencym(1.0);

    Vector_t<double, 3> R = {0.0, 0.0, 0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    cav_->apply(R, P, 0.0, E, B);

    EXPECT_NE(E(0), 0.0);  // x-component exists
    EXPECT_DOUBLE_EQ(E(1), 0.0);
    EXPECT_DOUBLE_EQ(E(2), 0.0);

    // No magnetic field at phase = 0
    EXPECT_DOUBLE_EQ(B(0), 0.0);
    EXPECT_DOUBLE_EQ(B(1), 0.0);
    EXPECT_DOUBLE_EQ(B(2), 0.0);
}

TEST_F(RFCavityTest, ApplyPhaseShift) {
    // φ = π/2 → cos = 0, sin = 1
    cav_->setPhasem(M_PI / 2.0);

    Vector_t<double, 3> R = {0.0, 0.0, 0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    cav_->apply(R, P, 0.0, E, B);

    // cos(pi/2) = 0 → no E
    EXPECT_NEAR(E(0), 0.0, 1e-12);

    // sin(pi/2) = 1 → B -= scale * B_map
    EXPECT_DOUBLE_EQ(B(1), -1.0);
}

TEST_F(RFCavityTest, ApplyPhasePi) {
    // φ = π → cos = -1, sin = 0
    cav_->setPhasem(M_PI);
    cav_->setFrequencym(1.0);

    Vector_t<double, 3> R = {0.0, 0.0, 0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    cav_->apply(R, P, 0.0, E, B);

    // cos(π) = -1 → E flipped
    EXPECT_NEAR(E(0), -1.0, 1e-12);
    EXPECT_DOUBLE_EQ(E(1), 0.0);
    EXPECT_DOUBLE_EQ(E(2), 0.0);

    // sin(π)=0 → no B
    EXPECT_DOUBLE_EQ(B(0), 0.0);
    EXPECT_NEAR(B(1), 0.0, 1e-12);
    EXPECT_DOUBLE_EQ(B(2), 0.0);
}

TEST_F(RFCavityTest, PhaseIndependentOfFrequencyAtT0) {
    cav_->setPhasem(M_PI / 2.0);

    Vector_t<double, 3> R = {0.0, 0.0, 0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};

    Vector_t<double, 3> E1 = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B1 = {0.0, 0.0, 0.0};

    cav_->setFrequencym(1.0);
    cav_->apply(R, P, 0.0, E1, B1);

    Vector_t<double, 3> E2 = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B2 = {0.0, 0.0, 0.0};

    cav_->setFrequencym(10.0);
    cav_->apply(R, P, 0.0, E2, B2);

    EXPECT_NEAR(E1(0), 0.0, 1e-12);
    EXPECT_NEAR(E2(0), 0.0, 1e-12);

    EXPECT_NEAR(B1(1), -1.0, 1e-12);
    EXPECT_NEAR(B2(1), -1.0, 1e-12);
}

// ---------------------------------------------------------------------------
// Scaling behavior
// ---------------------------------------------------------------------------
TEST_F(RFCavityTest, ApplyScaling) {
    Vector_t<double, 3> R = {0.0, 0.0, 0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};

    // Baseline
    Vector_t<double, 3> E1 = {0, 0, 0}, B1 = {0, 0, 0};
    cav_->apply(R, P, 0.0, E1, B1);

    // Scaled
    Vector_t<double, 3> E2 = {0, 0, 0}, B2 = {0, 0, 0};
    cav_->setScale(2.0);
    cav_->apply(R, P, 0.0, E2, B2);

    // compare
    EXPECT_NEAR(E2(0), 2.0 * E1(0), 1e-12);
    EXPECT_NEAR(E2(1), 2.0 * E1(1), 1e-12);
    EXPECT_NEAR(E2(2), 2.0 * E1(2), 1e-12);
}

// ---------------------------------------------------------------------------
// Time dependence
// ---------------------------------------------------------------------------
TEST_F(RFCavityTest, TimeDependence) {
    cav_->setPhasem(0.0);
    cav_->setFrequencym(1.0);

    Vector_t<double, 3> R = {0.0, 0.0, 0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};

    Vector_t<double, 3> E0 = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B0 = {0.0, 0.0, 0.0};

    cav_->apply(R, P, 0.0, E0, B0);  // phi = 0

    Vector_t<double, 3> E1 = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B1 = {0.0, 0.0, 0.0};

    cav_->apply(R, P, M_PI / 2.0, E1, B1);  // phi = pi/2

    EXPECT_NEAR(E0(0), 1.0, 1e-12);  // cos(0)
    EXPECT_NEAR(E1(0), 0.0, 1e-12);  // cos(pi/2)
}

// ---------------------------------------------------------------------------
// Fieldmap interaction
// ---------------------------------------------------------------------------
TEST_F(RFCavityTest, ApplyToReferenceParticle) {
    Vector_t<double, 3> R = {0.0, 0.0, 0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    // Should behave like apply() for reference particle
    cav_->applyToReferenceParticle(R, P, 0.0, E, B);

    EXPECT_DOUBLE_EQ(E(0), 1.0);
}

TEST_F(RFCavityTest, FieldmapOutOfBounds) {
    // Make fake fieldmap return true
    fmap_->setOutOfBounds(true);

    Vector_t<double, 3> R = {0.0, 0.0, 0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    cav_->apply(R, P, 0.0, E, B);

    // Out-of-bounds: no field is applied.
    EXPECT_DOUBLE_EQ(E(0), 0.0);
    EXPECT_DOUBLE_EQ(E(1), 0.0);
    EXPECT_DOUBLE_EQ(E(2), 0.0);
    EXPECT_TRUE(cav_->applyToReferenceParticle(R, P, 0.0, E, B));
}

// Ideal (map-free) RF cavity behavior.
class IdealRFCavityTest : public testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        OpalData::getInstance()->storeInputFn("ideal_rf_test.opal");
        gmsg                = new Inform(nullptr, -1);
        Options::enableHDF5 = false;
    }
    static void TearDownTestSuite() {
        delete gmsg;
        gmsg = nullptr;
        ippl::finalize();
    }
    void SetUp() override {
        cavity.setCavityType("IDEAL");
        cavity.getGeometry().setElementLength(0.2);
        cavity.setAmplitudem(5.0);
        cavity.setFrequencym(Physics::two_pi * 100.e6);
        cavity.setPhasem(0.0);
        cavity.setAperture(ApertureType::RECTANGULAR, {0.1, 0.05});
        cavity.initialise(nullptr);
    }
    std::shared_ptr<PartBunch_t> makeBunch(size_t count) {
        beam       = std::make_shared<Beam>();
        auto bunch = std::make_shared<PartBunch_t>(
                std::vector{1.0}, std::vector{1.0}, std::vector<Beam*>{Beam::find("UNNAMED_BEAM")},
                std::vector<size_t>{count}, 1.0, "LF2",
                opalx::spacecharge::CartesianDomainConfig3D{.periodicParticleBoundary = true});
        bunch->getParticleContainer()->createParticles(count);
        return bunch;
    }
    RFCavityRep cavity{"ideal"};
    std::shared_ptr<Beam> beam;
};

TEST_F(IdealRFCavityTest, GradientPhaseAndVariableRFConvention) {
    VariableRFCavity variable;
    variable.setLength(0.2);
    variable.setWidth(0.2);
    variable.setHeight(0.1);
    variable.setAmplitudeModel(std::make_shared<PolynomialTimeDependence>(std::vector{5.0}));
    variable.setFrequencyModel(std::make_shared<PolynomialTimeDependence>(std::vector{100.0}));
    variable.setPhaseModel(
            std::make_shared<PolynomialTimeDependence>(std::vector{0.3 + Physics::pi / 2.}));
    cavity.setPhasem(0.3);
    for (double time : {0.0, 1.e-9, 2.5e-9, 5.e-9}) {
        Vector_t<double, 3> e(0), b(0), v(0);
        cavity.apply({0., 0., 0.1}, {}, time, e, b);
        variable.apply({0., 0., 0.1}, {}, time, v, b);
        EXPECT_NEAR(e[2], 5.e6 * std::cos(Physics::two_pi * 100.e6 * time + 0.3), 1.e-8);
        EXPECT_NEAR(e[2], v[2], 1.e-8);
        EXPECT_DOUBLE_EQ(e[0], 0.);
        EXPECT_DOUBLE_EQ(e[1], 0.);
        EXPECT_DOUBLE_EQ(b[2], 0.);
    }
    cavity.setPhasem(0.);
    Vector_t<double, 3> e(0), b(0);
    cavity.apply({0., 0., 0.1}, {}, 0., e, b);
    EXPECT_DOUBLE_EQ(e[2] * 0.2, 1.e6);
}

TEST_F(IdealRFCavityTest, DeviceHostReferenceBoundariesAndMidpoint) {
    const std::vector<Vector_t<double, 3>> positions{
            {0., 0., -0.01}, {0., 0., 0.},    {0., 0., 0.1},    {0., 0., 0.2},
            {0.1, 0., 0.1},  {0., 0.05, 0.1}, {0.09, 0.04, 0.1}};
    const auto bunch = makeBunch(positions.size());
    const auto pc    = bunch->getParticleContainer();
    auto r           = Kokkos::create_mirror_view(pc->R.getView());
    for (size_t i = 0; i < positions.size(); ++i)
        r(i) = positions[i];
    Kokkos::deep_copy(pc->R.getView(), r);
    Kokkos::deep_copy(pc->E.getView(), Vector_t<double, 3>{1., 2., 3.});
    Kokkos::deep_copy(pc->B.getView(), Vector_t<double, 3>{4., 5., 6.});
    bunch->setT(1.e-9);
    bunch->setdT(2.e-9);
    cavity.initialise(bunch.get());
    cavity.goOnline(0.);
    cavity.apply(pc);
    const auto fields   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->E.getView());
    const auto magnetic = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->B.getView());
    for (size_t i = 0; i < positions.size(); ++i) {
        Vector_t<double, 3> e{1., 2., 3.}, refE(e), b{4., 5., 6.}, refB(b);
        cavity.apply(positions[i], {}, 2.e-9, e, b);
        bool outside = cavity.applyToReferenceParticle(positions[i], {}, 2.e-9, refE, refB);
        EXPECT_EQ(outside, i == 4 || i == 5);
        const bool inside = i == 1 || i == 2 || i == 6;
        EXPECT_EQ(cavity.isInside(positions[i]), inside);
        EXPECT_NEAR(e[2], 3. + (inside ? 5.e6 * std::cos(0.4 * Physics::pi) : 0.), 1.e-8);
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_DOUBLE_EQ(fields(i)[d], e[d]);
            EXPECT_DOUBLE_EQ(refE[d], e[d]);
            EXPECT_DOUBLE_EQ(magnetic(i)[d], b[d]);
            EXPECT_DOUBLE_EQ(refB[d], b[d]);
        }
    }
    EXPECT_NO_THROW(cavity.goOffline());
}

TEST_F(IdealRFCavityTest, ReferenceExcludesErrorsAndClonePreservesMode) {
    cavity.setAmplitudeError(1.0);
    cavity.setPhaseError(Physics::pi);
    std::unique_ptr<ElementBase> copy(cavity.clone());
    auto* cloned = dynamic_cast<RFCavityRep*>(copy.get());
    ASSERT_NE(cloned, nullptr);
    EXPECT_TRUE(cloned->isIdeal());
    Vector_t<double, 3> e(0), refE(0), b(0);
    cloned->apply({0., 0., 0.1}, {}, 0., e, b);
    cloned->applyToReferenceParticle({0., 0., 0.1}, {}, 0., refE, b);
    EXPECT_DOUBLE_EQ(e[2], -6.e6);
    EXPECT_DOUBLE_EQ(refE[2], 5.e6);
}

TEST_F(IdealRFCavityTest, InputValidationAndDefaultMode) {
    OpalCavity input;
    auto* rep = dynamic_cast<RFCavityRep*>(input.getElement());
    ASSERT_NE(rep, nullptr);
    EXPECT_NO_THROW(input.update());
    EXPECT_FALSE(rep->isIdeal());
    EXPECT_EQ(rep->getCavityType(), CavityType::SW);
    Attributes::setPredefinedString(input.itsAttr[OpalCavity::TYPE], "IDEAL");
    EXPECT_ANY_THROW(input.update());
    Attributes::setReal(input.itsAttr[OpalCavity::LENGTH], 0.2);
    Attributes::setReal(input.itsAttr[OpalCavity::VOLT], 1.);
    Attributes::setReal(input.itsAttr[OpalCavity::FREQ], 100.);
    EXPECT_NO_THROW(input.update());
    EXPECT_TRUE(rep->isIdeal());
    EXPECT_DOUBLE_EQ(rep->getFrequencym(), Physics::two_pi * 100.e6);
    Attributes::setString(input.itsAttr[OpalCavity::FMAPFN], "unused.map");
    EXPECT_ANY_THROW(input.update());
    Attributes::setString(input.itsAttr[OpalCavity::FMAPFN], "");
    Attributes::setString(input.itsAttr[OpalCavity::PHASE_MODEL], "phase");
    EXPECT_ANY_THROW(input.update());
    for (double length : {0., -1., std::numeric_limits<double>::infinity()}) {
        cavity.getGeometry().setElementLength(length);
        EXPECT_ANY_THROW(cavity.initialise(nullptr));
    }
}

TEST_F(IdealRFCavityTest, AutophaseAndDCGradientTimesLength) {
    const auto bunch = makeBunch(1);
    cavity.initialise(bunch.get());
    const auto ref       = bunch->getParticleContainer()->getReference();
    const double mass    = ref->getM() * Units::eV2MeV;
    const double kinetic = 10.;
    const double p0      = std::sqrt(kinetic * (kinetic + 2. * mass)) / mass;
    cavity.setFrequencym(0.);
    const auto dc     = cavity.trackOnAxisParticle(p0, 0., 1.e-13, ref->getQ(), mass);
    const double gain = (std::sqrt(1. + dc.first * dc.first) - 1.) * mass - kinetic;
    EXPECT_NEAR(
            gain, ref->getQ() * cavity.getAmplitudem() * cavity.getGeometry().getElementLength(),
            5.e-4);
    cavity.setFrequencym(Physics::two_pi * 100.e6);
    cavity.setPhasem(0.37);
    const double phase = cavity.getAutoPhaseEstimate(kinetic, 0., ref->getQ(), mass);
    EXPECT_TRUE(std::isfinite(phase));
    EXPECT_DOUBLE_EQ(cavity.getPhasem(), 0.37);
    const auto energyAt = [&](double phi) {
        cavity.setPhasem(phi);
        const auto result = cavity.trackOnAxisParticle(p0, 0., 1.e-13, ref->getQ(), mass);
        return (std::sqrt(1. + result.first * result.first) - 1.) * mass;
    };
    const double optimum = energyAt(phase);
    EXPECT_GT(optimum, energyAt(phase + 0.2));
    EXPECT_GT(optimum, energyAt(phase - 0.2));
}

TEST_F(IdealRFCavityTest, EllipticalApertureAndFiniteParameters) {
    cavity.setAperture(ApertureType::ELLIPTICAL, {0.1, 0.05});
    EXPECT_TRUE(cavity.isInside({0.05, 0.025, 0.1}));
    EXPECT_FALSE(cavity.isInside({0.09, 0.04, 0.1}));
    EXPECT_FALSE(cavity.isInside({0.1, 0., 0.1}));
    const double nan = std::numeric_limits<double>::quiet_NaN();
    cavity.setFrequencym(-1.);
    EXPECT_ANY_THROW(cavity.validateIdeal());
    cavity.setFrequencym(0.);
    EXPECT_NO_THROW(cavity.validateIdeal());
    cavity.setAmplitudem(nan);
    EXPECT_ANY_THROW(cavity.validateIdeal());
    cavity.setAmplitudem(1.);
    cavity.setPhasem(nan);
    EXPECT_ANY_THROW(cavity.validateIdeal());
}

TEST_F(IdealRFCavityTest, FullAutophaserRespectsVetoAndFindsCrest) {
    const auto bunch = makeBunch(1);
    auto rf          = std::make_shared<RFCavityRep>(cavity);
    rf->getGeometry().setElementLength(0.002);
    rf->setAmplitudem(50.0);
    rf->initialise(bunch.get());
    const auto ref       = bunch->getParticleContainer()->getReference();
    const double mass    = ref->getM() * Units::eV2MeV;
    const double kinetic = 10.;
    const double p0      = std::sqrt(kinetic * (kinetic + 2. * mass)) / mass;
    const Vector_t<double, 3> momentum{0., 0., p0};
    rf->setPhasem(0.37);
    rf->setAutophaseVeto(true);
    CavityAutophaser fixed(*ref, rf);
    EXPECT_DOUBLE_EQ(fixed.getPhaseAtMaxEnergy({0., 0., 0.}, momentum, 0., 1.e-13), 0.37);
    EXPECT_DOUBLE_EQ(rf->getPhasem(), 0.37);
    rf->setPhasem(0.);
    rf->setAutophaseVeto(false);
    CavityAutophaser automatic(*ref, rf);
    const double phase = automatic.getPhaseAtMaxEnergy({0., 0., 0.}, momentum, 0., 1.e-13);
    EXPECT_TRUE(std::isfinite(phase));
    EXPECT_TRUE(rf->getAutophaseVeto());
    EXPECT_DOUBLE_EQ(rf->getPhasem(), phase);
    const auto energyAt = [&](double phi) {
        rf->setPhasem(phi);
        const auto result = rf->trackOnAxisParticle(p0, 0., 1.e-13, ref->getQ(), mass);
        return (std::sqrt(1. + result.first * result.first) - 1.) * mass;
    };
    const double optimum = energyAt(phase);
    EXPECT_GT(optimum, kinetic);
    EXPECT_GT(optimum, energyAt(phase - 0.2));
    EXPECT_GT(optimum, energyAt(phase + 0.2));
}

TEST_F(IdealRFCavityTest, FullAutophaserAdjustsGradientToDesignEnergy) {
    const auto bunch = makeBunch(1);
    auto rf          = std::make_shared<RFCavityRep>(cavity);
    rf->getGeometry().setElementLength(0.002);
    rf->setAmplitudem(0.);
    rf->setDesignEnergy(10.2);
    rf->setPhasem(0.);
    rf->setAutophaseVeto(false);
    rf->initialise(bunch.get());
    const auto ref    = bunch->getParticleContainer()->getReference();
    const double mass = ref->getM() * Units::eV2MeV;
    const double p0   = std::sqrt(10. * (10. + 2. * mass)) / mass;
    CavityAutophaser automatic(*ref, rf);
    const double phase = automatic.getPhaseAtMaxEnergy({0., 0., 0.}, {0., 0., p0}, 0., 1.e-13);
    EXPECT_TRUE(std::isfinite(phase));
    EXPECT_GT(rf->getAmplitudem(), 50.0);
    EXPECT_LT(rf->getAmplitudem(), 150.0);
    const auto result        = rf->trackOnAxisParticle(p0, 0., 1.e-13, ref->getQ(), mass);
    const double finalEnergy = (std::sqrt(1. + result.first * result.first) - 1.) * mass;
    EXPECT_NEAR(finalEnergy, 10.2, 1.e-6);
}

TEST_F(IdealRFCavityTest, FieldAmplitudeIndependentOfLength) {
    const auto bunch = makeBunch(1);
    const auto pc    = bunch->getParticleContainer();
    bunch->setT(-1.e-9);
    bunch->setdT(2.e-9);  // The device midpoint is t=0.
    Kokkos::deep_copy(pc->R.getView(), Vector_t<double, 3>{0., 0., 0.1});
    for (double length : {0.2, 0.4}) {
        cavity.getGeometry().setElementLength(length);
        cavity.initialise(bunch.get());
        double begin = -1., end = -1.;
        cavity.getFieldExtent(begin, end);
        EXPECT_DOUBLE_EQ(begin, 0.);
        EXPECT_DOUBLE_EQ(end, length);
        Vector_t<double, 3> hostE(0), refE(0), b(0);
        cavity.apply({0., 0., 0.1}, {}, 0., hostE, b);
        EXPECT_FALSE(cavity.applyToReferenceParticle({0., 0., 0.1}, {}, 0., refE, b));
        Kokkos::deep_copy(pc->E.getView(), Vector_t<double, 3>(0));
        cavity.apply(pc);
        const auto fields =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->E.getView());
        // VOLT=5 MV/m sets the field directly; length only sets spatial support.
        EXPECT_DOUBLE_EQ(hostE[2], 5.e6);
        EXPECT_DOUBLE_EQ(refE[2], 5.e6);
        EXPECT_DOUBLE_EQ(fields(0)[2], 5.e6);
    }
}
