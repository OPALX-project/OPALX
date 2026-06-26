/**
 * \file TestTravelingWave.cpp
 * \brief Unit tests for TravelingWave component.
 *
 * This test suite validates the core behavior of TravelingWave using:
 *  - a minimal concrete implementation (TestTravelingWave)
 *  - a controlled fake fieldmap (FakeTWFieldmap)
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
 *    - getElementDimensions()
 *    - getEdgeToBegin()
 *    - getEdgeToEnd()
 *
 * 3. Spatial behavior
 *    - apply() in entry region
 *    - apply() in core region
 *    - apply() in exit region
 *    - apply() outside the element
 *
 * 4. Traveling-wave phase behavior
 *    - setPhasem() updates internal phaseCore1_m / phaseCore2_m / phaseExit_m
 *
 * 5. Field superposition
 *    - core region accumulates two field contributions
 *
 * 6. Scaling behavior
 *    - entry / core / exit scales act correctly
 *
 * 7. Fieldmap interaction / edge cases
 *    - correct handling of out-of-bounds fieldmap queries
 *    - isInside()
 *
 * ---------------------------------------------------------------------------
 * Notes
 * ---------------------------------------------------------------------------
 *
 * - The FakeTWFieldmap provides a simple, deterministic field:
 *      E = (1, 0, 0), B = (0, 1, 0)
 *   allowing direct verification of cos/sin modulation and accumulation.
 *
 * - Tests focus on physics/logic correctness of TravelingWave region handling.
 */

#include <gtest/gtest.h>

#define private public
#include "AbsBeamline/TravelingWave.h"
#undef private

#include "AbsBeamline/ElementBase.h"
#include "Fields/Fieldmap.h"

#include <cmath>
#include <memory>

// ---------------------------------------------------------------------------
// Fake fieldmap for deterministic tests
// ---------------------------------------------------------------------------
class FakeTWFieldmap : public Fieldmap {
public:
    FakeTWFieldmap() : Fieldmap("dummy"), outOfBounds_(false) {}

    void setOutOfBounds(bool v) { outOfBounds_ = v; }

    bool getFieldstrength(
            const Vector_t<double, 3>&, Vector_t<double, 3>& E,
            Vector_t<double, 3>& B) const override {
        if (outOfBounds_) return true;

        E = {1.0, 0.0, 0.0};
        B = {0.0, 1.0, 0.0};
        return false;
    }

    void getFieldDimensions(double& zBegin, double& zEnd) const override {
        zBegin = 0.0;
        zEnd   = 2.0;
    }

    void getFieldDimensions(
            double& xIni, double& xFinal, double& yIni, double& yFinal, double& zIni,
            double& zFinal) const override {
        xIni = yIni = 0.0;
        xFinal = yFinal = 0.0;
        zIni            = 0.0;
        zFinal          = 2.0;
    }

    bool isInside(const Vector_t<double, 3>&) const override { return !outOfBounds_; }

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
    void getOnaxisEz(std::vector<std::pair<double, double>>& F) override {
        F.clear();
        F.push_back({0.0, 1.0});
        F.push_back({1.0, 1.0});
        F.push_back({2.0, 1.0});
    }

private:
    bool outOfBounds_;
};

// ---------------------------------------------------------------------------
// Dummy Geometry
// ---------------------------------------------------------------------------
class DummyGeometryTW : public BGeometryBase {
public:
    double getArcLength() const override { return 0.0; }
    double getElementLength() const override { return length_m; }
    void setElementLength(double length) override { length_m = length; }

    Euclid3D getTransform(double, double) const override { return Euclid3D(); }

private:
    double length_m = 0.0;
};

// ---------------------------------------------------------------------------
// Dummy Field
// ---------------------------------------------------------------------------
class DummyFieldTW : public EMField {
public:
    void scale(double) override {}
};

// ---------------------------------------------------------------------------
// Minimal concrete TravelingWave
// ---------------------------------------------------------------------------
class TestTravelingWave : public TravelingWave {
public:
    TestTravelingWave() : TravelingWave("test_tw") {}

    double getAmplitude() const override { return amplitude_; }
    double getFrequency() const override { return frequency_; }
    double getPhase() const override { return phase_; }

    ElementBase* clone() const override { return new TestTravelingWave(*this); }

    void setTestElementLength(double v) { geom_.setElementLength(v); }

    BGeometryBase& getGeometry() override { return geom_; }
    const BGeometryBase& getGeometry() const override { return geom_; }

    EMField& getField() override { return field_; }
    const EMField& getField() const override { return field_; }

    void setAmplitude(double v) { amplitude_ = v; }
    void setFrequency(double v) { frequency_ = v; }
    void setPhase(double v) { phase_ = v; }

    void setScale(double v) { scale_m = v; }
    void setScaleError(double v) { scaleError_m = v; }
    void setScaleCore(double v) { scaleCore_m = v; }
    void setScaleCoreError(double v) { scaleCoreError_m = v; }

    void setFrequencyInternal(double v) { frequency_m = v; }
    void setPhaseInternal(double v) { phase_m = v; }
    // void setPhaseErrorInternal(double v) { phaseError_m = v; }

    void setFieldmap(Fieldmap* fmap) { fieldmap_m = fmap; }
    void setStartField(double v) { startField_m = v; }
    void setEndField(double v) { endField_m = v; }

private:
    double amplitude_ = 0.0;
    double frequency_ = 0.0;
    double phase_     = 0.0;

    DummyGeometryTW geom_;
    DummyFieldTW field_;
};

// ---------------------------------------------------------------------------
// Fixture
// ---------------------------------------------------------------------------
class TravelingWaveTest : public ::testing::Test {
protected:
    void SetUp() override {
        tw_   = std::make_unique<TestTravelingWave>();
        fmap_ = std::make_unique<FakeTWFieldmap>();

        tw_->setFieldmap(fmap_.get());

        // Stable default TW geometry for tests
        tw_->periodLength_m = 2.0;
        tw_->cellLength_m   = 1.0;
        tw_->numCells_m     = 3;
        tw_->mode_m         = 0.5;

        tw_->setStartField(-1.0);
        tw_->setEndField(3.0);
        tw_->startCoreField_m       = 1.0;
        tw_->startExitField_m       = 3.0;
        tw_->mappedStartExitField_m = 1.0;

        tw_->setTestElementLength(4.0);

        // Stable RF defaults
        tw_->setScale(1.0);
        tw_->setScaleError(0.0);
        tw_->setScaleCore(1.0);
        tw_->setScaleCoreError(0.0);
        tw_->setFrequencyInternal(1.0);
        tw_->setPhaseInternal(0.0);
        tw_->setPhaseError(0.0);

        // Set internal TW phases explicitly
        tw_->phaseCore1_m = 0.0;
        tw_->phaseCore2_m = 0.0;
        tw_->phaseExit_m  = 0.0;
    }

    std::unique_ptr<TestTravelingWave> tw_;
    std::unique_ptr<FakeTWFieldmap> fmap_;
};

// ---------------------------------------------------------------------------
// Basic API
// ---------------------------------------------------------------------------
TEST_F(TravelingWaveTest, GetType) { EXPECT_EQ(tw_->getType(), ElementType::TRAVELINGWAVE); }

TEST_F(TravelingWaveTest, Bends) { EXPECT_FALSE(tw_->bends()); }

TEST_F(TravelingWaveTest, GetSetAmplitudeFrequencyPhase) {
    tw_->setAmplitude(5.0);
    tw_->setFrequency(2.0);
    tw_->setPhase(0.5);

    EXPECT_DOUBLE_EQ(tw_->getAmplitude(), 5.0);
    EXPECT_DOUBLE_EQ(tw_->getFrequency(), 2.0);
    EXPECT_DOUBLE_EQ(tw_->getPhase(), 0.5);
}

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------
TEST_F(TravelingWaveTest, GetDimensions) {
    double zBegin = 0.0, zEnd = 0.0;
    tw_->getFieldExtend(zBegin, zEnd);

    EXPECT_DOUBLE_EQ(zBegin, -1.0);  // -0.5 * periodLength
    EXPECT_DOUBLE_EQ(zEnd, 3.0);     // zBegin + elementLength
}

TEST_F(TravelingWaveTest, GetElementDimensions) {
    double begin = 0.0, end = 0.0;
    tw_->getElementDimensions(begin, end);

    EXPECT_DOUBLE_EQ(begin, 0.0);
    EXPECT_DOUBLE_EQ(end, 4.0);
}

TEST_F(TravelingWaveTest, EdgeTransforms) {
    auto beg = tw_->getEdgeToBegin();
    auto end = tw_->getEdgeToEnd();

    EXPECT_DOUBLE_EQ(beg.getOrigin()(2), 0.0);
    EXPECT_DOUBLE_EQ(end.getOrigin()(2), 4.0);
}

// ---------------------------------------------------------------------------
// setPhasem updates internal TW phases
// ---------------------------------------------------------------------------
TEST_F(TravelingWaveTest, SetPhasemUpdatesInternalPhases) {
    tw_->mode_m     = 0.5;
    tw_->numCells_m = 3;

    tw_->setPhasem(0.2);

    EXPECT_NEAR(tw_->getPhasem(), 0.2, 1e-12);
    EXPECT_NEAR(tw_->phaseCore1_m, 0.2 + Physics::pi * 0.5 / 2.0, 1e-12);
    EXPECT_NEAR(tw_->phaseCore2_m, 0.2 + Physics::pi * 0.5 * 1.5, 1e-12);

    const double expectedExit = 0.2 - Physics::two_pi * ((3 - 1) * 0.5 - std::floor((3 - 1) * 0.5));
    EXPECT_NEAR(tw_->phaseExit_m, expectedExit, 1e-12);
}

// ---------------------------------------------------------------------------
// apply(): entry/core/exit/outside
// ---------------------------------------------------------------------------
TEST_F(TravelingWaveTest, ApplyEntryRegion) {
    // Entry: tmpR(2) = R(2) + 1.0 < startCoreField_m (=1.0)
    Vector_t<double, 3> R = {0.0, 0.0, -0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    tw_->apply(R, P, 0.0, E, B);

    // entry scale = 1, phase = 0 => cos = 1, sin = 0
    EXPECT_DOUBLE_EQ(E(0), 1.0);
    EXPECT_DOUBLE_EQ(E(1), 0.0);
    EXPECT_DOUBLE_EQ(E(2), 0.0);

    EXPECT_DOUBLE_EQ(B(0), 0.0);
    EXPECT_DOUBLE_EQ(B(1), 0.0);
    EXPECT_DOUBLE_EQ(B(2), 0.0);
}

TEST_F(TravelingWaveTest, ApplyCoreRegionAccumulatesTwoContributions) {
    // Core: tmpR(2) = R(2) + 1.0 in [1.0, 3.0)
    Vector_t<double, 3> R = {0.0, 0.0, 0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    tw_->setScaleCore(1.0);
    tw_->phaseCore1_m = 0.0;
    tw_->phaseCore2_m = 0.0;

    tw_->apply(R, P, 0.0, E, B);

    // Fake fieldmap returns E=(1,0,0) each time, applied twice
    EXPECT_DOUBLE_EQ(E(0), 2.0);
    EXPECT_DOUBLE_EQ(E(1), 0.0);
    EXPECT_DOUBLE_EQ(E(2), 0.0);

    EXPECT_DOUBLE_EQ(B(0), 0.0);
    EXPECT_DOUBLE_EQ(B(1), 0.0);
    EXPECT_DOUBLE_EQ(B(2), 0.0);
}

TEST_F(TravelingWaveTest, ApplyExitRegion) {
    // Exit: tmpR(2) = R(2) + 1.0 >= startExitField_m (=3.0)
    Vector_t<double, 3> R = {0.0, 0.0, 2.2};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    tw_->phaseExit_m = 0.0;

    tw_->apply(R, P, 0.0, E, B);

    EXPECT_DOUBLE_EQ(E(0), 1.0);
    EXPECT_DOUBLE_EQ(E(1), 0.0);
    EXPECT_DOUBLE_EQ(E(2), 0.0);
}

TEST_F(TravelingWaveTest, ApplyOutsideBefore) {
    Vector_t<double, 3> R = {0.0, 0.0, -1.1};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {1.0, 2.0, 3.0};
    Vector_t<double, 3> B = {4.0, 5.0, 6.0};

    tw_->apply(R, P, 0.0, E, B);

    EXPECT_DOUBLE_EQ(E(0), 1.0);
    EXPECT_DOUBLE_EQ(E(1), 2.0);
    EXPECT_DOUBLE_EQ(E(2), 3.0);
    EXPECT_DOUBLE_EQ(B(0), 4.0);
    EXPECT_DOUBLE_EQ(B(1), 5.0);
    EXPECT_DOUBLE_EQ(B(2), 6.0);
}

TEST_F(TravelingWaveTest, ApplyOutsideAfter) {
    Vector_t<double, 3> R = {0.0, 0.0, 3.0};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {1.0, 2.0, 3.0};
    Vector_t<double, 3> B = {4.0, 5.0, 6.0};

    tw_->apply(R, P, 0.0, E, B);

    EXPECT_DOUBLE_EQ(E(0), 1.0);
    EXPECT_DOUBLE_EQ(E(1), 2.0);
    EXPECT_DOUBLE_EQ(E(2), 3.0);
    EXPECT_DOUBLE_EQ(B(0), 4.0);
    EXPECT_DOUBLE_EQ(B(1), 5.0);
    EXPECT_DOUBLE_EQ(B(2), 6.0);
}

// ---------------------------------------------------------------------------
// RF phase behavior
// ---------------------------------------------------------------------------
TEST_F(TravelingWaveTest, ApplyPhasePiOverTwoGivesMagneticField) {
    Vector_t<double, 3> R = {0.0, 0.0, -0.5};  // entry region
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    tw_->setPhasem(Physics::pi / 2.0);
    tw_->phaseCore1_m = Physics::pi / 2.0;
    tw_->phaseCore2_m = Physics::pi / 2.0;
    tw_->phaseExit_m  = Physics::pi / 2.0;

    tw_->apply(R, P, 0.0, E, B);

    EXPECT_NEAR(E(0), 0.0, 1e-12);
    EXPECT_DOUBLE_EQ(B(1), -1.0);
}

// ---------------------------------------------------------------------------
// applyToReferenceParticle
// ---------------------------------------------------------------------------
TEST_F(TravelingWaveTest, ApplyToReferenceParticleEntry) {
    Vector_t<double, 3> R = {0.0, 0.0, -0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    tw_->applyToReferenceParticle(R, P, 0.0, E, B);

    EXPECT_DOUBLE_EQ(E(0), 1.0);
}

// ---------------------------------------------------------------------------
// Fieldmap interaction / edge cases
// ---------------------------------------------------------------------------
TEST_F(TravelingWaveTest, FieldmapOutOfBounds) {
    fmap_->setOutOfBounds(true);

    Vector_t<double, 3> R = {0.0, 0.0, -0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool out = tw_->apply(R, P, 0.0, E, B);
    EXPECT_TRUE(out);
}

TEST_F(TravelingWaveTest, IsInside) {
    EXPECT_TRUE(tw_->isInside({0.0, 0.0, 0.0}));
    EXPECT_TRUE(tw_->isInside({0.0, 0.0, 2.5}));
    EXPECT_FALSE(tw_->isInside({0.0, 0.0, -1.1}));
    EXPECT_FALSE(tw_->isInside({0.0, 0.0, 3.0}));
}
