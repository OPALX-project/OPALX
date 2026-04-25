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

        // Analytic bend field coefficients are scaled by the OPAL reference
        // momentum P0. Construct the standard REAL variable environment so P0
        // exists during parser/update tests.
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
    Attributes::setReal(ui.itsAttr[OpalSBend::LENGTH], 2.0);
    Attributes::setReal(ui.itsAttr[OpalSBend::K0], 0.2);

    EXPECT_THROW(ui.update(), OpalException);
}

TEST_F(TestOpalBend, SBendRejectsInconsistentK0) {
    OpalSBend ui;
    Attributes::setReal(ui.itsAttr[OpalSBend::LENGTH], 2.0);
    Attributes::setReal(ui.itsAttr[OpalSBend::ANGLE], 0.4);
    Attributes::setReal(ui.itsAttr[OpalSBend::K0], 0.25);

    EXPECT_THROW(ui.update(), OpalException);
}

TEST_F(TestOpalBend, SBendAcceptsConsistentK0ButUsesAngleAsPrimaryDefinition) {
    OpalSBend ui;
    Attributes::setReal(ui.itsAttr[OpalSBend::LENGTH], 2.0);
    Attributes::setReal(ui.itsAttr[OpalSBend::ANGLE], 0.4);
    Attributes::setReal(ui.itsAttr[OpalSBend::K0], 0.2);

    EXPECT_NO_THROW(ui.update());

    const auto* bend = dynamic_cast<const SBendRep*>(ui.getElement());
    ASSERT_NE(bend, nullptr);
    EXPECT_NEAR(bend->getGeometry().getElementLength(), 2.0, 1.0e-12);
    EXPECT_NEAR(bend->getGeometry().getBendAngle(), 0.4, 1.0e-12);
    EXPECT_NEAR(bend->getFieldAmplitude(), 0.2, 1.0e-12);
}

TEST_F(TestOpalBend, RBendRequiresAngleAndUsesAngleOverK0) {
    OpalRBend invalidUi;
    Attributes::setReal(invalidUi.itsAttr[OpalRBend::LENGTH], 1.5);
    Attributes::setReal(invalidUi.itsAttr[OpalRBend::K0], 0.3);
    EXPECT_THROW(invalidUi.update(), OpalException);

    OpalRBend validUi;
    Attributes::setReal(validUi.itsAttr[OpalRBend::LENGTH], 1.5);
    Attributes::setReal(validUi.itsAttr[OpalRBend::ANGLE], 0.3);
    Attributes::setReal(validUi.itsAttr[OpalRBend::K0], 0.2);

    EXPECT_NO_THROW(validUi.update());

    const auto* bend = dynamic_cast<const RBendRep*>(validUi.getElement());
    ASSERT_NE(bend, nullptr);
    EXPECT_NEAR(bend->getGeometry().getElementLength(), 1.5, 1.0e-12);
    EXPECT_NEAR(bend->getGeometry().getBendAngle(), 0.3, 1.0e-12);
    EXPECT_NEAR(bend->getFieldAmplitude(), 0.2, 1.0e-12);
}
