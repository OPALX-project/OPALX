//
// Unit tests for analytic SBEND/RBEND element definitions
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//

#include "Attributes/Attributes.h"
#include "BeamlineCore/RBendRep.h"
#include "BeamlineCore/SBendRep.h"
#include "Elements/OpalRBend.h"
#include "Elements/OpalSBend.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "ValueDefinitions/RealVariable.h"
#include "gtest/gtest.h"

extern Inform* gmsg;

class TestOpalBend : public testing::Test {
public:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        gmsg = new Inform(nullptr, -1);

        // Construct the standard REAL variable environment used by parser-
        // driven element updates. The analytic bend runtime normalization no
        // longer depends on parser-global P0, but the surrounding attribute
        // infrastructure still expects the REAL variable registry to exist.
        realVariableFactory_m = std::make_unique<RealVariable>();
    }

    static void TearDownTestSuite() {
        realVariableFactory_m.reset();
        delete gmsg;
        gmsg = nullptr;
        ippl::finalize();
    }

protected:
    static std::unique_ptr<RealVariable> realVariableFactory_m;
};

std::unique_ptr<RealVariable> TestOpalBend::realVariableFactory_m;

TEST_F(TestOpalBend, SBendRejectsK0OnlyDefinition) {
    OpalSBend ui;
    std::unique_ptr<OpalSBend> instance(ui.clone("SB1"));
    Attributes::setReal(instance->itsAttr[OpalSBend::LENGTH], 2.0);
    Attributes::setReal(instance->itsAttr[OpalSBend::K0], 0.2);

    EXPECT_THROW(instance->update(), OpalException);
}

TEST_F(TestOpalBend, SBendRejectsInconsistentK0) {
    OpalSBend ui;
    std::unique_ptr<OpalSBend> instance(ui.clone("SB1"));
    Attributes::setReal(instance->itsAttr[OpalSBend::LENGTH], 2.0);
    Attributes::setReal(instance->itsAttr[OpalSBend::ANGLE], 0.4);
    Attributes::setReal(instance->itsAttr[OpalSBend::K0], 0.25);

    EXPECT_THROW(instance->update(), OpalException);
}

TEST_F(TestOpalBend, SBendAcceptsConsistentK0ButUsesAngleAsPrimaryDefinition) {
    OpalSBend ui;
    std::unique_ptr<OpalSBend> instance(ui.clone("SB1"));
    Attributes::setReal(instance->itsAttr[OpalSBend::LENGTH], 2.0);
    Attributes::setReal(instance->itsAttr[OpalSBend::ANGLE], 0.4);
    Attributes::setReal(instance->itsAttr[OpalSBend::K0], 0.2);

    EXPECT_NO_THROW(instance->update());

    const auto* bend = dynamic_cast<const SBendRep*>(instance->getElement());
    ASSERT_NE(bend, nullptr);
    EXPECT_NEAR(bend->getGeometry().getElementLength(), 2.0, 1.0e-12);
    EXPECT_NEAR(bend->getGeometry().getBendAngle(), 0.4, 1.0e-12);
    EXPECT_NEAR(bend->getFieldAmplitude(), 0.2, 1.0e-12);
}

TEST_F(TestOpalBend, SBendFringeSupportExtendsFieldExtentAndRenormalizesDipoleCoefficient) {
    OpalSBend ui;
    std::unique_ptr<OpalSBend> instance(ui.clone("SB4"));
    Attributes::setReal(instance->itsAttr[OpalSBend::LENGTH], 1.0);
    Attributes::setReal(instance->itsAttr[OpalSBend::ANGLE], 0.4);
    Attributes::setReal(instance->itsAttr[OpalSBend::HGAP], 0.02);
    Attributes::setReal(instance->itsAttr[OpalSBend::FINT], 0.5);

    EXPECT_NO_THROW(instance->update());

    const auto* bend = dynamic_cast<const SBendRep*>(instance->getElement());
    ASSERT_NE(bend, nullptr);

    double zBegin = 0.0;
    double zEnd   = 0.0;
    bend->getFieldExtend(zBegin, zEnd);
    EXPECT_NEAR(zBegin, -0.01, 1.0e-12);
    EXPECT_NEAR(zEnd, 1.01, 1.0e-12);
    EXPECT_NEAR(bend->getEffectiveFieldLength(), 1.01, 1.0e-12);
    EXPECT_NEAR(bend->getFieldAmplitude(), 0.4 / 1.01, 1.0e-12);
}

TEST_F(TestOpalBend, RBendRequiresAngleAndUsesAngleOverK0) {
    OpalRBend ui;
    std::unique_ptr<OpalRBend> invalidInstance(ui.clone("RB1"));
    Attributes::setReal(invalidInstance->itsAttr[OpalRBend::LENGTH], 1.5);
    Attributes::setReal(invalidInstance->itsAttr[OpalRBend::K0], 0.3);
    EXPECT_THROW(invalidInstance->update(), OpalException);

    std::unique_ptr<OpalRBend> validInstance(ui.clone("RB2"));
    const double length       = 1.5;
    const double angle        = 0.3;
    const double consistentK0 = angle / (length * (angle / 2.0) / std::sin(angle / 2.0));
    Attributes::setReal(validInstance->itsAttr[OpalRBend::LENGTH], 1.5);
    Attributes::setReal(validInstance->itsAttr[OpalRBend::ANGLE], 0.3);
    Attributes::setReal(validInstance->itsAttr[OpalRBend::K0], consistentK0);

    EXPECT_NO_THROW(validInstance->update());

    const auto* bend = dynamic_cast<const RBendRep*>(validInstance->getElement());
    ASSERT_NE(bend, nullptr);
    EXPECT_NEAR(bend->getGeometry().getElementLength(), 1.5, 1.0e-12);
    EXPECT_NEAR(bend->getGeometry().getBendAngle(), 0.3, 1.0e-12);
    EXPECT_NEAR(bend->getFieldAmplitude(), consistentK0, 1.0e-12);
}

TEST_F(TestOpalBend, RBendFringeSupportUsesGapFallbackForHalfGap) {
    OpalRBend ui;
    std::unique_ptr<OpalRBend> instance(ui.clone("RB4"));
    Attributes::setReal(instance->itsAttr[OpalRBend::LENGTH], 1.2);
    Attributes::setReal(instance->itsAttr[OpalRBend::ANGLE], 0.3);
    Attributes::setReal(instance->itsAttr[OpalRBend::GAP], 0.06);
    Attributes::setReal(instance->itsAttr[OpalRBend::FINT], 0.5);

    EXPECT_NO_THROW(instance->update());

    const auto* bend = dynamic_cast<const RBendRep*>(instance->getElement());
    ASSERT_NE(bend, nullptr);

    double zBegin = 0.0;
    double zEnd   = 0.0;
    bend->getFieldExtend(zBegin, zEnd);
    EXPECT_LT(zBegin, 0.0);
    EXPECT_GT(zEnd, 1.2);
    EXPECT_NEAR(bend->getFringeHalfGap(), 0.03, 1.0e-12);
    EXPECT_NEAR(bend->getFieldAmplitude(), 0.3 / bend->getEffectiveFieldLength(), 1.0e-12);
}

TEST_F(TestOpalBend, ExemplarUpdateDoesNotEnforceInstanceOnlyAngleChecks) {
    OpalSBend sbendExemplar;
    OpalRBend rbendExemplar;

    EXPECT_NO_THROW(sbendExemplar.update());
    EXPECT_NO_THROW(rbendExemplar.update());
}

TEST_F(TestOpalBend, SBendWithoutHapertKeepsNonzeroFallbackAperture) {
    OpalSBend ui;
    std::unique_ptr<OpalSBend> instance(ui.clone("SB3"));
    Attributes::setReal(instance->itsAttr[OpalSBend::LENGTH], 0.5);
    Attributes::setReal(instance->itsAttr[OpalSBend::ANGLE], Physics::pi / 8.0);

    EXPECT_NO_THROW(instance->update());

    const auto* bend = dynamic_cast<const SBendRep*>(instance->getElement());
    ASSERT_NE(bend, nullptr);
    const auto aperture = bend->getAperture();
    EXPECT_GT(aperture.second.size(), 1u);
    EXPECT_GT(aperture.second[0], 0.0);
    EXPECT_GT(aperture.second[1], 0.0);
}

TEST_F(TestOpalBend, RBendWithoutHapertKeepsNonzeroFallbackAperture) {
    OpalRBend ui;
    std::unique_ptr<OpalRBend> instance(ui.clone("RB3"));
    Attributes::setReal(instance->itsAttr[OpalRBend::LENGTH], 0.5);
    Attributes::setReal(instance->itsAttr[OpalRBend::ANGLE], Physics::pi / 8.0);

    EXPECT_NO_THROW(instance->update());

    const auto* bend = dynamic_cast<const RBendRep*>(instance->getElement());
    ASSERT_NE(bend, nullptr);
    const auto aperture = bend->getAperture();
    EXPECT_GT(aperture.second.size(), 1u);
    EXPECT_GT(aperture.second[0], 0.0);
    EXPECT_GT(aperture.second[1], 0.0);
}
