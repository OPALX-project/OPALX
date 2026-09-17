//
// Unit tests for class OpalEnge
//
// Copyright (c) 2017-2026, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
// All rights reserved.
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.

#include "Attributes/Attributes.h"
#include "gtest/gtest.h"
#include "AbsBeamline/EndFieldModel/EndFieldModelManager.h"
#include "Elements/OpalAsymmetricEnge.h"
#include "AbsBeamline/EndFieldModel/AsymmetricEnge.h"

TEST(TestOpalAsymmetricEnge, TestSetup) {
    // Make the UI
    OpalAsymmetricEnge ui;
    // Set the attributes
    Attributes::setReal(ui.itsAttr[OpalAsymmetricEnge::X0_START], 7);
    Attributes::setReal(ui.itsAttr[OpalAsymmetricEnge::LAMBDA_START], 8);
    Attributes::setRealArray(ui.itsAttr[OpalAsymmetricEnge::COEFFICIENTS_START], {101.0, 3.0, 4.0});
    Attributes::setReal(ui.itsAttr[OpalAsymmetricEnge::X0_END], 9);
    Attributes::setReal(ui.itsAttr[OpalAsymmetricEnge::LAMBDA_END], 11);
    Attributes::setRealArray(ui.itsAttr[OpalAsymmetricEnge::COEFFICIENTS_END], {12.0, 17.0, 21.0});
    ui.update();
    auto efmMan = endfieldmodel::EndFieldModelManager::getEFMManager();
    EXPECT_NO_THROW(efmMan->getEndFieldModel("ASYMMETRIC_ENGE"));
    auto efm = efmMan->getEndFieldModel("ASYMMETRIC_ENGE");
    auto enge = std::dynamic_pointer_cast<endfieldmodel::AsymmetricEnge,
                                          endfieldmodel::EndFieldModel>(efm);
    endfieldmodel::AsymmetricEngeConfig config = enge->getConfig();
    EXPECT_EQ(config.engeStart_m.x0_m, 7.0);
    EXPECT_EQ(config.engeStart_m.lambda_m, 8.0);
    EXPECT_EQ(config.engeStart_m.a_m, std::vector<double>({101.0, 3.0, 4.0}));

    EXPECT_EQ(config.engeEnd_m.x0_m, 9.0);
    EXPECT_EQ(config.engeEnd_m.lambda_m, 11.0);
    EXPECT_EQ(config.engeEnd_m.a_m, std::vector<double>({12.0, 17.0, 21.0}));

}

