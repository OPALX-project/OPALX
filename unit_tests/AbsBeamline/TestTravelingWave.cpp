/**
 * \file TestTravelingWave.cpp
 * \brief Unit tests for TravelingWave component.
 *
 * This test suite validates the core behavior of TravelingWave using:
 *  - a subclass with test access to protected state (TestTravelingWave)
 *  - a controlled fake fieldmap (FakeTWFieldmap)
 *
 * ---------------------------------------------------------------------------
 * Coverage
 * ---------------------------------------------------------------------------
 *
 * 1. Basic API
 *    - getType()
 *
 * 2. Geometry
 *    - getFieldExtent()
 *    - getGeometry().getEdgeToBegin()
 *    - getGeometry().getEdgeToEnd()
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

#include "../TemporaryDirectory.h"
#include "Fields/Astra1DDynamic.h"
#include "Ippl.h"

#include <fstream>

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
    void applyField(
            std::shared_ptr<ParticleContainer_t>, double = 1.0, double = 1.0,
            const ElementRegions& = ElementRegions()) override {}
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
// ---------------------------------------------------------------------------
// TravelingWave with test access to protected state
// ---------------------------------------------------------------------------
class TestTravelingWave : public TravelingWave {
public:
    TestTravelingWave() : TravelingWave("test_tw") {}

    void setTestElementLength(double v) { getGeometry().setElementLength(v); }

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
};

// ---------------------------------------------------------------------------
// Fixture
// ---------------------------------------------------------------------------
class TravelingWaveTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }
    static void TearDownTestSuite() { ippl::finalize(); }

    void TearDown() override {
        Fieldmap::clearDictionary();
        directory_m.reset();
    }

    /**
     * Give the element an Astra map with a constant Ez over [0, 4] m, the map the kernel
     * formula Astra1DDynamic::computeTravelingWaveField() reads. Returns its Ez on axis [V/m].
     */
    double useConstantAstraMap() {
        directory_m.create("opalx-traveling-wave");
        const std::string file = (directory_m.path() / "tw.astra").string();
        {
            std::ofstream f(file);
            f << "AstraDynamic 8 TRUE\n1.0\n";
            for (int i = 0; i <= 40; ++i)
                f << 0.1 * i << " 1.0\n";
        }
        Fieldmap* map = Fieldmap::getFieldmap(file);
        Fieldmap::readMap(file);
        tw_->setFieldmap(map);
        Vector_t<double, 3> e(0.0), b(0.0);
        map->getFieldstrength({0.0, 0.0, 1.0}, e, b);
        return e(2);
    }

    opalx::test::TemporaryDirectory directory_m;

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

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------
TEST_F(TravelingWaveTest, GetDimensions) {
    double zBegin = 0.0, zEnd = 0.0;
    tw_->getFieldExtent(zBegin, zEnd);

    EXPECT_DOUBLE_EQ(zBegin, -1.0);  // -0.5 * periodLength
    EXPECT_DOUBLE_EQ(zEnd, 3.0);     // zBegin + elementLength
}

TEST_F(TravelingWaveTest, GetElementDimensions) {
    double begin = 0.0, end = 0.0;
    begin = 0.0;
    end   = tw_->getGeometry().getElementLength();

    EXPECT_DOUBLE_EQ(begin, 0.0);
    EXPECT_DOUBLE_EQ(end, 4.0);
}

TEST_F(TravelingWaveTest, EdgeTransforms) {
    auto beg = tw_->getGeometry().getEdgeToBegin();
    auto end = tw_->getGeometry().getEdgeToEnd();

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
    const double ez = useConstantAstraMap();
    // Entry: R(2) + periodLength/2 = 0.5 < startCoreField_m (= 1.0)
    Vector_t<double, 3> R = {0.0, 0.0, -0.5};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    tw_->applyOnHost(R, 0.0, E, B);

    // entry scale = 1, phase = 0 => cos = 1, sin = 0
    // The Fourier sum of the map is constant to round-off.
    EXPECT_DOUBLE_EQ(E(0), 0.0);
    EXPECT_DOUBLE_EQ(E(1), 0.0);
    EXPECT_NEAR(E(2), ez, 1e-9 * ez);

    EXPECT_DOUBLE_EQ(B(0), 0.0);
    EXPECT_DOUBLE_EQ(B(1), 0.0);
    EXPECT_DOUBLE_EQ(B(2), 0.0);
}

TEST_F(TravelingWaveTest, ApplyCoreRegionAccumulatesTwoContributions) {
    const double ez = useConstantAstraMap();
    // Core: R(2) + periodLength/2 = 1.5 in [1.0, 3.0)
    Vector_t<double, 3> R = {0.0, 0.0, 0.5};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    tw_->setScaleCore(1.0);
    tw_->phaseCore1_m = 0.0;
    tw_->phaseCore2_m = 0.0;

    tw_->applyOnHost(R, 0.0, E, B);

    // The constant map gives Ez at both core positions, added twice.
    EXPECT_NEAR(E(2), 2.0 * ez, 1e-9 * ez);
    EXPECT_DOUBLE_EQ(B(1), 0.0);
}

TEST_F(TravelingWaveTest, ApplyExitRegion) {
    const double ez = useConstantAstraMap();
    // Exit: R(2) + periodLength/2 = 3.2 >= startExitField_m (= 3.0); mapped to 3.2 - 1.0 = 2.2.
    Vector_t<double, 3> R = {0.0, 0.0, 2.2};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    tw_->phaseExit_m = 0.0;

    tw_->applyOnHost(R, 0.0, E, B);

    EXPECT_NEAR(E(2), ez, 1e-9 * ez);
}

TEST_F(TravelingWaveTest, ApplyOutsideBefore) {
    Vector_t<double, 3> R = {0.0, 0.0, -1.1};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {1.0, 2.0, 3.0};
    Vector_t<double, 3> B = {4.0, 5.0, 6.0};

    tw_->applyOnHost(R, 0.0, E, B);

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

    tw_->applyOnHost(R, 0.0, E, B);

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
    const double ez = useConstantAstraMap();
    // Entry region, 1 cm off axis, where the map has a magnetic field.
    Vector_t<double, 3> R = {0.01, 0.0, -0.5};
    Vector_t<double, 3> E0(0.0), B0(0.0), E(0.0), B(0.0);

    tw_->applyOnHost(R, 0.0, E0, B0);  // phase 0: cos = 1, sin = 0
    EXPECT_DOUBLE_EQ(B0(1), 0.0);

    tw_->setPhasem(Physics::pi / 2.0);
    tw_->phaseCore1_m = Physics::pi / 2.0;
    tw_->phaseCore2_m = Physics::pi / 2.0;
    tw_->phaseExit_m  = Physics::pi / 2.0;

    tw_->applyOnHost(R, 0.0, E, B);

    EXPECT_NEAR(E(2), 0.0, 1e-12 * ez);
    EXPECT_NE(B(1), 0.0);
}

// ---------------------------------------------------------------------------
// Fieldmap interaction / edge cases
// ---------------------------------------------------------------------------
TEST_F(TravelingWaveTest, PointMappedOutsideTheMapGetsNoField) {
    useConstantAstraMap();
    // Exit: 3.2 - mappedStartExitField_m (= -1.5) = 4.7 lies past the map end at 4.0.
    tw_->mappedStartExitField_m = -1.5;
    Vector_t<double, 3> R       = {0.0, 0.0, 2.2};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    tw_->applyOnHost(R, 0.0, E, B);

    // No field, and inside the aperture this is a field-free gap, not material.
    EXPECT_DOUBLE_EQ(E(0), 0.0);
    EXPECT_DOUBLE_EQ(E(1), 0.0);
    EXPECT_DOUBLE_EQ(E(2), 0.0);
    EXPECT_FALSE(tw_->getRegions().inLossRegion(R));
}

TEST_F(TravelingWaveTest, IsInside) {
    EXPECT_TRUE(tw_->isInside({0.0, 0.0, 0.0}));
    EXPECT_TRUE(tw_->isInside({0.0, 0.0, 2.5}));
    EXPECT_FALSE(tw_->isInside({0.0, 0.0, -1.1}));
    EXPECT_FALSE(tw_->isInside({0.0, 0.0, 3.0}));
}
