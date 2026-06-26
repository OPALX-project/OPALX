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
 *    - bends()
 *    - amplitude / frequency / phase setters and getters
 *
 * 2. Geometry
 *    - getFieldExtend()
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
    void applyField(std::shared_ptr<ParticleContainer_t>, double = 1.0) override {}
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
// Dummy Geometry (fully concrete)
// ---------------------------------------------------------------------------
class DummyGeometry : public BGeometryBase {
public:
    double getArcLength() const override { return 0.0; }
    double getElementLength() const override { return length_m; }
    void setElementLength(double length) override { length_m = length; }

    Euclid3D getTransform(double, double) const override { return Euclid3D(); }

private:
    double length_m = 0.0;
};

// ---------------------------------------------------------------------------
// Dummy Field (fully concrete)
// ---------------------------------------------------------------------------
class DummyField : public EMField {
public:
    void scale(double) override {}
};

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

    BGeometryBase& getGeometry() override { return geom_; }
    const BGeometryBase& getGeometry() const override { return geom_; }

    EMField& getField() override { return field_; }
    const EMField& getField() const override { return field_; }

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

    DummyGeometry geom_;
    DummyField field_;
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
        cav_->setElementLength(1.0);

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

TEST_F(RFCavityTest, Bends) { EXPECT_FALSE(cav_->bends()); }

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

    cav_->getFieldExtend(zBegin, zEnd);

    EXPECT_EQ(zBegin, 0.0);
    EXPECT_EQ(zEnd, 1.0);
}

TEST_F(RFCavityTest, BodyExtentCanDifferFromFieldSupport) {
    cav_->setStartField(0.2);
    cav_->setEndField(0.8);
    cav_->setElementLength(1.0);

    double bodyBegin = -1.0, bodyEnd = -1.0;
    cav_->getElementDimensions(bodyBegin, bodyEnd);
    EXPECT_EQ(bodyBegin, 0.0);
    EXPECT_EQ(bodyEnd, 1.0);

    const auto entry = cav_->getEdgeToBegin();
    const auto exit  = cav_->getEdgeToEnd();
    EXPECT_EQ(entry.getOrigin()(2), 0.0);
    EXPECT_EQ(exit.getOrigin()(2), 1.0);

    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};
    EXPECT_FALSE(cav_->apply({0.0, 0.0, 0.1}, {0.0, 0.0, 1.0}, 0.0, E, B));
    EXPECT_DOUBLE_EQ(E(0), 0.0);

    EXPECT_FALSE(cav_->apply({0.0, 0.0, 0.5}, {0.0, 0.0, 1.0}, 0.0, E, B));
    EXPECT_DOUBLE_EQ(E(0), 1.0);
}

TEST_F(RFCavityTest, ZeroBodyLengthDoesNotFallBackToFieldmapLength) {
    cav_->setElementLength(0.0);

    EXPECT_DOUBLE_EQ(cav_->getElementLength(), 0.0);

    double bodyBegin = -1.0, bodyEnd = -1.0;
    cav_->getElementDimensions(bodyBegin, bodyEnd);
    EXPECT_DOUBLE_EQ(bodyBegin, 0.0);
    EXPECT_DOUBLE_EQ(bodyEnd, 0.0);

    double fieldBegin = -1.0, fieldEnd = -1.0;
    cav_->getFieldExtend(fieldBegin, fieldEnd);
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

    bool out = cav_->apply(R, P, 0.0, E, B);

    EXPECT_TRUE(out);
}
