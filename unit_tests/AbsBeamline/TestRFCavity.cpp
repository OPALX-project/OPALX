/**
 * \file TestRFCavity.cpp
 * \brief Unit tests for RFCavity component (base layer).
 *
 * Tests RFCavity core API and apply logic using a lightweight concrete subclass
 * and a fake Fieldmap.
 *
 * Covers:
 *  - getType, bends
 *  - amplitude / frequency / phase setters + getters
 *  - getDimensions
 *  - apply(R,P,t,E,B): inside, outside, boundaries
 *  - phase dependence (cos/sin behavior)
 */
#include <gtest/gtest.h>

#include "Fields/Fieldmap.h"
#include "AbsBeamline/ElementBase.h"
#include "BeamlineCore/RFCavityRep.h"

#include "Fields/EMField.h"
#include "BeamlineGeometry/Geometry.h"

#include <cmath>
#include <memory>
// #include <fstream> 

// class DummyGeometry : public BGeometryBase {
// public:
//     DummyGeometry() : BGeometryBase() {}
// };

// class DummyField : public EMField {
// public:
//     DummyField() : EMField() {}
// };

class FakeFieldmap : public Fieldmap {
public:
    FakeFieldmap() : Fieldmap("dummy") {}

    bool getFieldstrength(const Vector_t<double,3>&,
                          Vector_t<double,3>& E,
                          Vector_t<double,3>& B) const override {
        E = {1.0, 0.0, 0.0};
        B = {0.0, 1.0, 0.0};
        return false; // inside
    }

    void getFieldDimensions(double& zBegin, double& zEnd) const override {
        zBegin = 0.0;
        zEnd   = 1.0;
    }

    void getFieldDimensions(double& zBegin, double& zEnd,
                            double&, double&, double&, double&) const override {
        zBegin = 0.0;
        zEnd   = 1.0;
    }

    bool isInside(const Vector_t<double,3>&) const override {
        return true;
    }

    // --- required no-op implementations ---
    void applyField(std::shared_ptr<ParticleContainer_t>, double = 1.0) override {}
    bool getFieldDerivative(const Vector_t<double,3>&,
                            Vector_t<double,3>&,
                            Vector_t<double,3>&,
                            const DiffDirection&) const override { return false; }

    void swap() override {}
    void getInfo(Inform*) override {}
    double getFrequency() const override { return 1.0; }
    void setFrequency(double) override {}
    void readMap() override {}
    void freeMap() override {}
};

// ---------------------------------------------------------------------------
// Dummy Geometry (fully concrete)
// ---------------------------------------------------------------------------
class DummyGeometry : public BGeometryBase {
public:
    double getArcLength() const override { return 0.0; }
    double getElementLength() const override { return 0.0; }

    Euclid3D getTransform(double, double) const override {
        return Euclid3D();
    }
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

    ElementBase* clone() const override {
        return new TestRFCavity(*this);
    }

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
        cav_ = std::make_unique<TestRFCavity>();
        fmap_ = std::make_unique<FakeFieldmap>();

        cav_->setScale(2.0);
        cav_->setFrequencyInternal(1.0);
        cav_->setPhaseInternal(0.0);

        cav_->setFieldmap(fmap_.get());
        cav_->setStartField(0.0);
        cav_->setElementLength(1.0);
    }

    std::unique_ptr<TestRFCavity> cav_;
    std::unique_ptr<FakeFieldmap> fmap_;
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------
TEST_F(RFCavityTest, GetType) {
    EXPECT_EQ(cav_->getType(), ElementType::RFCAVITY);
}

TEST_F(RFCavityTest, Bends) {
    EXPECT_FALSE(cav_->bends());
}

TEST_F(RFCavityTest, GetSetAmplitudeFrequencyPhase) {
    cav_->setAmplitude(5.0);
    cav_->setFrequency(2.0);
    cav_->setPhase(0.5);

    EXPECT_DOUBLE_EQ(cav_->getAmplitude(), 5.0);
    EXPECT_DOUBLE_EQ(cav_->getFrequency(), 2.0);
    EXPECT_DOUBLE_EQ(cav_->getPhase(), 0.5);
}

TEST_F(RFCavityTest, GetDimensions) {
    double zBegin = -1.0, zEnd = -1.0;

    cav_->getDimensions(zBegin, zEnd);

    EXPECT_EQ(zBegin, 0.0);
    EXPECT_EQ(zEnd, 0.0);
}

TEST_F(RFCavityTest, ApplyInside) {
    Vector_t<double,3> R = {0.0, 0.0, 0.5};
    Vector_t<double,3> P = {0.0, 0.0, 1.0};
    Vector_t<double,3> E = {0.0, 0.0, 0.0};
    Vector_t<double,3> B = {0.0, 0.0, 0.0};

    cav_->apply(R, P, 0.0, E, B);

    // cos(0) = 1 → E += scale * E_map
    EXPECT_DOUBLE_EQ(E(0), 2.0); // 2 * 1
    EXPECT_DOUBLE_EQ(E(1), 0.0);
    EXPECT_DOUBLE_EQ(E(2), 0.0);

    // sin(0) = 0 → no B contribution
    EXPECT_DOUBLE_EQ(B(0), 0.0);
    EXPECT_DOUBLE_EQ(B(1), 0.0);
    EXPECT_DOUBLE_EQ(B(2), 0.0);
}

TEST_F(RFCavityTest, ApplyBefore) {
    Vector_t<double,3> R = {0.0, 0.0, -0.1};
    Vector_t<double,3> P = {0.0, 0.0, 1.0};
    Vector_t<double,3> E = {1.0, 1.0, 1.0};
    Vector_t<double,3> B = {1.0, 1.0, 1.0};

    cav_->apply(R, P, 0.0, E, B);

    EXPECT_DOUBLE_EQ(E(0), 1.0);
    EXPECT_DOUBLE_EQ(B(1), 1.0);
}

TEST_F(RFCavityTest, ApplyAfter) {
    Vector_t<double,3> R = {0.0, 0.0, 1.5};
    Vector_t<double,3> P = {0.0, 0.0, 1.0};
    Vector_t<double,3> E = {1.0, 1.0, 1.0};
    Vector_t<double,3> B = {1.0, 1.0, 1.0};

    cav_->apply(R, P, 0.0, E, B);

    EXPECT_DOUBLE_EQ(E(0), 1.0);
    EXPECT_DOUBLE_EQ(B(1), 1.0);
}

TEST_F(RFCavityTest, ApplyPhaseShift) {
    cav_->setPhasem(M_PI / 2.0); // 90 degrees

    Vector_t<double,3> R = {0.0, 0.0, 0.5};
    Vector_t<double,3> P = {0.0, 0.0, 1.0};
    Vector_t<double,3> E = {0.0, 0.0, 0.0};
    Vector_t<double,3> B = {0.0, 0.0, 0.0};

    cav_->apply(R, P, 0.0, E, B);

    // cos(pi/2) = 0 → no E
    EXPECT_NEAR(E(0), 0.0, 1e-12);

    // sin(pi/2) = 1 → B -= scale * B_map
    EXPECT_DOUBLE_EQ(B(1), -2.0);
}

TEST_F(RFCavityTest, ApplyToReferenceParticle) {
    Vector_t<double,3> R = {0.0, 0.0, 0.5};
    Vector_t<double,3> P = {0.0, 0.0, 1.0};
    Vector_t<double,3> E = {0.0, 0.0, 0.0};
    Vector_t<double,3> B = {0.0, 0.0, 0.0};

    cav_->applyToReferenceParticle(R, P, 0.0, E, B);

    EXPECT_DOUBLE_EQ(E(0), 2.0);
}
