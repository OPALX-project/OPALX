//
// Unit tests for class OpalTanh
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
#include "Elements/OpalTanh.h"
#include "AbsBeamline/EndFieldModel/Tanh.h"

TEST(TestOpalTanh, TestSetup) {
    // Make the UI
    OpalTanh ui;
    // Set the attributes
    Attributes::setReal(ui.itsAttr[OpalTanh::X0], 4);
    Attributes::setReal(ui.itsAttr[OpalTanh::LAMBDA], 2);
    ui.update();
    auto efmMan = endfieldmodel::EndFieldModelManager::getEFMManager();
    EXPECT_NO_THROW(efmMan->getEndFieldModel("TANH"));
    auto efm = efmMan->getEndFieldModel("TANH");
    EXPECT_EQ(efm->getCentreLength(), 8.0);
    EXPECT_EQ(efm->getEndLength(), 2.0);
}

