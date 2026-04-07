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

#include <cmath>
#include <memory>

// std::string writeRFCavityFieldmap(const std::string& path)
// {
//     std::ofstream f(path);

//     f << "2DDynamic XZ\n";
//     f << "0.0 1.0 1\n";
//     f << "100.0\n";
//     f << "0.0 0.0 1\n";

//     // minimal field (Ez only)
//     f << "1.0 0.0 0.0 0.0\n";

//     return path;
// }

// namespace {

// // ---------------------------------------------------------------------------
// // Minimal concrete RFCavity (since RFCavity is abstract)
// // ---------------------------------------------------------------------------
// class TestRFCavity : public RFCavity {
// public:
//     TestRFCavity() : RFCavity("TestRFCavity") {}

//     // Required RFCavity interface
//     double getAmplitude() const override { return scale_m; }
//     double getFrequency() const override { return frequency_m; }
//     double getPhase() const override { return phase_m; }

//     // ---- setters for protected members ----
//     void setFieldmap(Fieldmap* fmap) { fieldmap_m = fmap; }
//     void setStartField(double val) { startField_m = val; }

//     // ---- required abstract stubs ----
//     ElementBase* clone() const override {
//         return new TestRFCavity(*this);
//     }

//     BGeometryBase& getGeometry() override {
//         throw std::runtime_error("Not implemented");
//     }

//     const BGeometryBase& getGeometry() const override {
//         throw std::runtime_error("Not implemented");
//     }

//     EMField& getField() override {
//         throw std::runtime_error("Not implemented");
//     }

//     const EMField& getField() const override {
//         throw std::runtime_error("Not implemented");
//     }
// };

// ---------------------------------------------------------------------------
// Test fixture
// ---------------------------------------------------------------------------
// class RFCavityTest : public ::testing::Test {
// protected:
//     void SetUp() override {
//         cav_ = std::make_unique<TestRFCavity>();

//         cav_->setAmplitudem(2.0);
//         cav_->setFrequencym(1.0);
//         cav_->setPhasem(0.0);

//         auto file = writeRFCavityFieldmap("rf_test.dat");

//         cav_->setFieldMapFN(file);

//         cav_->setStartField(0.0);
//         cav_->setElementLength(1.0);
//     }

//     std::unique_ptr<TestRFCavity> cav_;
// };

// ---------------------------------------------------------------------------
// Type and geometry
// ---------------------------------------------------------------------------
// TEST_F(RFCavityTest, GetType) {
//     EXPECT_EQ(cav_->getType(), ElementType::RFCAVITY);
// }

// TEST_F(RFCavityTest, Bends) {
//     EXPECT_FALSE(cav_->bends());
// }

// TEST_F(RFCavityTest, GetSetAmplitudeFrequencyPhase)
// {
//     cav_->setAmplitudem(5.0);
//     cav_->setFrequencym(2.0);
//     cav_->setPhasem(0.5);

//     EXPECT_DOUBLE_EQ(cav_->getAmplitudem(), 5.0);
//     EXPECT_DOUBLE_EQ(cav_->getFrequencym(), 2.0);
//     EXPECT_DOUBLE_EQ(cav_->getPhasem(), 0.5);
// }

// TEST_F(RFCavityTest, GetDimensions)
// {
//     double zBegin = -1.0;
//     double zEnd   = -1.0;

//     cav_->initialise(nullptr, zBegin, zEnd);

//     cav_->getDimensions(zBegin, zEnd);

//     EXPECT_GT(zEnd, 0.0);
//     EXPECT_LE(zBegin, zEnd);
// }


// // ---------------------------------------------------------------------------
// // apply(R,P,t,E,B)
// // ---------------------------------------------------------------------------
// TEST_F(RFCavityTest, ApplyInside) {
//     Vector_t<double,3> R = {0.0, 0.0, 0.5};
//     Vector_t<double,3> P = {0.0, 0.0, 1.0};
//     Vector_t<double,3> E = {0.0, 0.0, 0.0};
//     Vector_t<double,3> B = {0.0, 0.0, 0.0};

//     cav_->apply(R, P, 0.0, E, B);

//     // cos(0)=1 → E scaled
//     EXPECT_DOUBLE_EQ(E(0), 2.0 * 1.0);
//     EXPECT_DOUBLE_EQ(E(1), 2.0 * 2.0);
//     EXPECT_DOUBLE_EQ(E(2), 2.0 * 3.0);

//     // sin(0)=0 → B unchanged
//     EXPECT_DOUBLE_EQ(B(0), 0.0);
//     EXPECT_DOUBLE_EQ(B(1), 0.0);
//     EXPECT_DOUBLE_EQ(B(2), 0.0);
// }

// TEST_F(RFCavityTest, ApplyOutsideBefore) {
//     Vector_t<double,3> R = {0.0, 0.0, -0.1};
//     Vector_t<double,3> P = {0.0, 0.0, 1.0};
//     Vector_t<double,3> E = {1.0, 1.0, 1.0};
//     Vector_t<double,3> B = {1.0, 1.0, 1.0};

//     cav_->apply(R, P, 0.0, E, B);

//     EXPECT_DOUBLE_EQ(E(0), 1.0);
//     EXPECT_DOUBLE_EQ(B(0), 1.0);
// }

// TEST_F(RFCavityTest, ApplyOutsideAfter) {
//     Vector_t<double,3> R = {0.0, 0.0, 1.5};
//     Vector_t<double,3> P = {0.0, 0.0, 1.0};
//     Vector_t<double,3> E = {1.0, 1.0, 1.0};
//     Vector_t<double,3> B = {1.0, 1.0, 1.0};

//     cav_->apply(R, P, 0.0, E, B);

//     EXPECT_DOUBLE_EQ(E(0), 1.0);
//     EXPECT_DOUBLE_EQ(B(0), 1.0);
// }

// TEST_F(RFCavityTest, ApplyAtStartBoundary) {
//     Vector_t<double,3> R = {0.0, 0.0, 0.0};
//     Vector_t<double,3> P = {0.0, 0.0, 1.0};
//     Vector_t<double,3> E = {0.0, 0.0, 0.0};
//     Vector_t<double,3> B = {0.0, 0.0, 0.0};

//     cav_->apply(R, P, 0.0, E, B);

//     EXPECT_DOUBLE_EQ(E(0), 2.0);
// }

// TEST_F(RFCavityTest, ApplyAtEndBoundaryExcluded) {
//     Vector_t<double,3> R = {0.0, 0.0, 1.0};
//     Vector_t<double,3> P = {0.0, 0.0, 1.0};
//     Vector_t<double,3> E = {0.0, 0.0, 0.0};
//     Vector_t<double,3> B = {0.0, 0.0, 0.0};

//     cav_->apply(R, P, 0.0, E, B);

//     // End is excluded: no effect
//     EXPECT_DOUBLE_EQ(E(0), 0.0);
// }

// // ---------------------------------------------------------------------------
// // Phase dependence
// // ---------------------------------------------------------------------------
// TEST_F(RFCavityTest, ApplyPhaseShiftPiOver2) {
//     cav_->setPhasem(M_PI / 2.0);

//     Vector_t<double,3> R = {0.0, 0.0, 0.5};
//     Vector_t<double,3> P = {0.0, 0.0, 1.0};
//     Vector_t<double,3> E = {0.0, 0.0, 0.0};
//     Vector_t<double,3> B = {0.0, 0.0, 0.0};

//     cav_->apply(R, P, 0.0, E, B);

//     // cos(pi/2)=0 → no E contribution
//     EXPECT_NEAR(E(0), 0.0, 1e-12);
//     EXPECT_NEAR(E(1), 0.0, 1e-12);
//     EXPECT_NEAR(E(2), 0.0, 1e-12);

//     // sin(pi/2)=1 → B affected
//     EXPECT_DOUBLE_EQ(B(0), -2.0 * 0.5);
//     EXPECT_DOUBLE_EQ(B(1), -2.0 * 1.0);
//     EXPECT_DOUBLE_EQ(B(2), -2.0 * 1.5);
// }



// } // namespace