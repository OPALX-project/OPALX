//
// Tests for MultipoleTBase
//
// Copyright (c) 2025, Jon Thompson, STFC Rutherford Appleton Laboratory, Didcot, UK
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//

#include <csignal>
#include "AbsBeamline/BeamlineVisitor.h"
#include "AbsBeamline/MultipoleT.h"
#include "Algorithms/SplineTimeDependence.h"
#include "gtest/gtest.h"

class TestMultipoleT : public testing::Test, public MultipoleT, public BeamlineVisitor {
public:
    TestMultipoleT() : MultipoleT("Magnet") {}

    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }
    static void TearDownTestSuite() { ippl::finalize(); }

    // Overrides of BeamlineVisitor
    void execute() override {}
    void visitBeamline(const Beamline&) override {}
    void visitComponent(const Component&) override {}
    void visitConstantEFieldCavity(const ConstantEFieldCavity&) override {}
    void visitDrift(const Drift&) override {}
    void visitFlaggedElmPtr(const FlaggedElmPtr&) override {}
    void visitMarker(const Marker&) override {}
    void visitMonitor(const Monitor&) override {}
    void visitMultipole(const Multipole&) override {}
    void visitMultipoleT(const MultipoleT&) override {}
    void visitOffset(const Offset&) override {}
    void visitRFCavity(const RFCavity&) override {}
    void visitScalingFFAMagnet(const ScalingFFAMagnet&) override {}
    void visitRing(const Ring&) override {}
    void visitSolenoid(const Solenoid&) override {}
    void visitTravelingWave(const TravelingWave&) override {}
    void visitVerticalFFAMagnet(const VerticalFFAMagnet&) override {}
    void visitProbe(const Probe&) override {}

    // Test helper functions
    double fieldAtT(const Vector_t<double, 3>& pos, const double t) {
        Vector_t<double, 3> E{}, B{};
        apply(pos, {}, t, E, B);
        return std::hypot(B[0], B[1], B[2]);
    }
};

// Does the clone API return an object with the same configuration
TEST_F(TestMultipoleT, Cloning) {
    // Set up the magnet
    setBendAngle(1, false);
    setElementLength(4);
    setAperture(3, 2);
    setFringeField(2, 1, 2);
    setRotation(0.4);
    setEntranceAngle(1);
    setMaxOrder(5, 10);
    setTransProfile({1, 2, 3});
    setBoundingBoxLength(7);
    setEntryOffset(3);
    setScalingName("Scaling");
    // Does the clone have the same configuration
    std::unique_ptr<MultipoleT> theClone;
    theClone.reset(dynamic_cast<MultipoleT*>(clone()));
    EXPECT_NE(theClone.get(), nullptr);
    EXPECT_EQ(theClone->getMaxFOrder(), getMaxFOrder());
    EXPECT_EQ(theClone->getMaxXOrder(), getMaxXOrder());
    EXPECT_EQ(theClone->getTransMaxOrder(), getTransMaxOrder());
    EXPECT_EQ(theClone->getTransProfile(), getTransProfile());
    EXPECT_EQ(theClone->getFringeField(), getFringeField());
    EXPECT_EQ(theClone->getEntryOffset(), getEntryOffset());
    EXPECT_EQ(theClone->getVariableRadius(), getVariableRadius());
    EXPECT_EQ(theClone->getBendAngle(), getBendAngle());
    EXPECT_EQ(theClone->getEntranceAngle(), getEntranceAngle());
    EXPECT_EQ(theClone->getLength(), getLength());
    EXPECT_EQ(theClone->getAperture(), getAperture());
    EXPECT_EQ(theClone->getRotation(), getRotation());
    EXPECT_EQ(theClone->getBoundingBoxLength(), getBoundingBoxLength());
    EXPECT_EQ(theClone->getScalingName(), getScalingName());
    EXPECT_EQ(theClone->getEntryOffset(), getEntryOffset());
}

// Attach a time dependency object and check it affects the field output appropriately
TEST_F(TestMultipoleT, TimeDependency) {
    // Set up a time dependency
    const auto timeDependence = std::make_shared<SplineTimeDependence>(
            1, std::vector<double>{0, 1, 1.01, 2}, std::vector<double>{1, 1, 2, 2});
    AbstractTimeDependence::setTimeDependence("STEPFUNCTION", timeDependence);
    // Set up the magnet
    setBendAngle(0, false);
    setElementLength(4);
    setAperture(0.3, 0.3);
    setFringeField(2, 0.2, 0.2);
    setMaxOrder(5, 10);
    setTransProfile({1});
    setScalingName("STEPFUNCTION");
    // This initialises the time dependency
    accept(*this);
    // Test the field at the magnet center at various times
    EXPECT_NEAR(fieldAtT({0.0, 0.0, 2.0}, 0.0), 1.0, 1e-6);
    EXPECT_NEAR(fieldAtT({0.0, 0.0, 2.0}, 1.0), 1.0, 1e-6);
    EXPECT_NEAR(fieldAtT({0.0, 0.0, 2.0}, 1.2), 2.0, 1e-6);
    EXPECT_NEAR(fieldAtT({0.0, 0.0, 2.0}, 2.0), 2.0, 1e-6);
    // Now clear the time dependence
    setScalingName("");
    // This initialises the time dependency
    accept(*this);
    // Test the field at the magnet center at various times
    EXPECT_NEAR(fieldAtT({0.0, 0.0, 2.0}, 0.0), 1.0, 1e-6);
    EXPECT_NEAR(fieldAtT({0.0, 0.0, 2.0}, 1.0), 1.0, 1e-6);
    EXPECT_NEAR(fieldAtT({0.0, 0.0, 2.0}, 1.2), 1.0, 1e-6);
    EXPECT_NEAR(fieldAtT({0.0, 0.0, 2.0}, 2.0), 1.0, 1e-6);
}

// Does the bends API return the correct value
TEST_F(TestMultipoleT, Bends) {
    setBendAngle(0, false);
    EXPECT_FALSE(bends());
    setBendAngle(1, false);
    EXPECT_TRUE(bends());
}

// Check that an exception if thrown for configuration that is not supported.
TEST_F(TestMultipoleT, ConfigurationValidation) {
    // Set up the magnet
    setBendAngle(0, false);
    setElementLength(4);
    setAperture(0.3, 0.3);
    setFringeField(2, 0.2, 0.2);
    setTransProfile({1});
    // Illegal F order
    setMaxOrder(10, 10);
    EXPECT_THROW(fieldAtT({0.0, 0.0, 2.0}, 0.0), OpalException);
}

// Tests for the few remaining API functions are collected here
TEST_F(TestMultipoleT, OddApis) {
    // The field object API which currently always returns a dummy object
    MultipoleT* magnet = this;
    auto* field        = &magnet->getField();
    EXPECT_NE(field, nullptr);
    auto* constField = &const_cast<const MultipoleT*>(magnet)->getField();
    EXPECT_NE(constField, nullptr);
    EXPECT_EQ(field, constField);
    // Just make sure these functions do not throw
    EXPECT_NO_THROW(finalise());
    double a, b;
    EXPECT_NO_THROW(getDimensions(a, b));
}
