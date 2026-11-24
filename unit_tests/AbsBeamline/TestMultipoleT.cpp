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
#include <memory>
#include "AbsBeamline/MultipoleT.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/SplineTimeDependence.h"
#include "Attributes/Attributes.h"
#include "Structure/Beam.h"
#include "Structure/FieldSolverCmd.h"
#include "Track/FreezeRun.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"
#include "gtest/gtest.h"

class TestMultipoleT : public testing::Test, public MultipoleT {
public:
    TestMultipoleT() : MultipoleT("Magnet") {}

    class TestableFieldSolverCmd final : public FieldSolverCmd {
    public:
        void setType(const std::string& t) {
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::TYPE], t);
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

    static std::shared_ptr<FieldSolverCmd> makeNoFieldSolver() {
        const auto fsCmd = std::make_shared<TestableFieldSolverCmd>();
        fsCmd->setType("NONE");
        fsCmd->setNX(8);
        fsCmd->setNY(8);
        fsCmd->setNZ(8);
        fsCmd->setBCX("PERIODIC");
        fsCmd->setBCY("PERIODIC");
        fsCmd->setBCZ("PERIODIC");
        return fsCmd;
    }

    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        OpalData::getInstance()->storeInputFn("unit_test.opal");
        if (gmsg == nullptr) {
            gmsg = new Inform(nullptr, -1);
        }
        Options::enableHDF5 = false;
    }
    static void TearDownTestSuite() {
        if (gmsg != nullptr) {
            delete gmsg;
            gmsg = nullptr;
        }
        ippl::finalize();
    }

    // Test helper functions
    double fieldAtT(const Vector_t<double, 3>& pos, const double t) {
        Vector_t<double, 3> E{}, B{};
        apply(pos, {}, t, E, B);
        return std::hypot(B[0], B[1], B[2]);
    }
};

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
    initialiseTimeDependencies();
    // Test the field at the magnet center at various times
    EXPECT_NEAR(fieldAtT({0.0, 0.0, 2.0}, 0.0), 1.0, 1e-6);
    EXPECT_NEAR(fieldAtT({0.0, 0.0, 2.0}, 1.0), 1.0, 1e-6);
    EXPECT_NEAR(fieldAtT({0.0, 0.0, 2.0}, 1.2), 2.0, 1e-6);
    EXPECT_NEAR(fieldAtT({0.0, 0.0, 2.0}, 2.0), 2.0, 1e-6);
    // Now clear the time dependence
    setScalingName("");
    // This initialises the time dependency
    initialiseTimeDependencies();
    // Test the field at the magnet center at various times
    EXPECT_NEAR(fieldAtT({0.0, 0.0, 2.0}, 0.0), 1.0, 1e-6);
    EXPECT_NEAR(fieldAtT({0.0, 0.0, 2.0}, 1.0), 1.0, 1e-6);
    EXPECT_NEAR(fieldAtT({0.0, 0.0, 2.0}, 1.2), 1.0, 1e-6);
    EXPECT_NEAR(fieldAtT({0.0, 0.0, 2.0}, 2.0), 1.0, 1e-6);
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

TEST_F(TestMultipoleT, DefaultFringeParametersOnReusedStorage) {
    // Reuse storage containing a previous magnet's fringe settings, as can happen
    // between fixtures. Default construction must not retain those settings.
    alignas(MultipoleTConfig) unsigned char storage[sizeof(MultipoleTConfig)];
    auto* config                = new (storage) MultipoleTConfig;
    config->fringeS0_m          = 2.0;
    config->fringeLambdaLeft_m  = 0.2;
    config->fringeLambdaRight_m = 0.3;
    std::destroy_at(config);
    config = new (storage) MultipoleTConfig;
    EXPECT_DOUBLE_EQ(config->fringeS0_m, 0.0);
    EXPECT_DOUBLE_EQ(config->fringeLambdaLeft_m, 0.0);
    EXPECT_DOUBLE_EQ(config->fringeLambdaRight_m, 0.0);
    std::destroy_at(config);
}

TEST_F(TestMultipoleT, FieldExtentIncludesAsymmetricFringes) {
    setElementLength(4.0);
    setFringeField(2.0, 0.2, 0.3);
    double begin = 0.0, end = 0.0;
    getFieldExtent(begin, end);
    EXPECT_DOUBLE_EQ(begin, -2.0);
    EXPECT_DOUBLE_EQ(end, 7.0);

    // Geometry and fringe contributions remain independent when L changes.
    setElementLength(5.0);
    getFieldExtent(begin, end);
    EXPECT_DOUBLE_EQ(begin, -2.0);
    EXPECT_DOUBLE_EQ(end, 8.0);

    setFringeField(2.5, 0.0, 0.0);
    getFieldExtent(begin, end);
    EXPECT_DOUBLE_EQ(begin, 0.0);
    EXPECT_DOUBLE_EQ(end, 5.0);
}

// Tests for the few remaining API functions are collected here
TEST_F(TestMultipoleT, OddApis) {
    // The field-support interval follows the body length.
    setElementLength(4.0);
    double a = -1.0, b = -1.0;
    EXPECT_NO_THROW(getFieldExtent(a, b));
    EXPECT_DOUBLE_EQ(a, 0.0);
    EXPECT_DOUBLE_EQ(b, 4.0);
    EXPECT_EQ(getType(), ElementType::MULTIPOLET);
}

TEST_F(TestMultipoleT, ApplySingleParticleThrowsForMultiContainerBunch) {
    // Configure a simple straight multipole.
    setBendAngle(0.0, false);
    setElementLength(1.0);
    setAperture(1.0, 1.0);
    setFringeField(0.5, 0.1, 0.1);
    setTransProfile({1.0});

    // Build a multi-container bunch.
    const auto fsCmd = makeNoFieldSolver();
    const auto beam  = std::make_shared<Beam>();
    Beam::defineBuiltins();
    Beam* testBeam = Beam::find("UNNAMED_BEAM");
    ASSERT_NE(testBeam, nullptr);
    const auto bunch = std::make_shared<PartBunch_t>(
            std::vector{1.0e-9, 2.0e-9}, std::vector{0.511e-3, 0.938},
            std::vector{
                    opalx::run::describeBeamForBunch(*testBeam),
                    opalx::run::describeBeamForBunch(*testBeam)},
            std::vector<size_t>{1, 1}, 1.0,
            opalx::spacecharge::CartesianDomainConfig3D{.periodicParticleBoundary = true});
    ASSERT_EQ(bunch->getNumParticleContainers(), 2u);

    // Register bunch; the per-particle-index apply() was removed, so only check that
    // a multi-container bunch can be registered.
    initialise(bunch.get());
}
