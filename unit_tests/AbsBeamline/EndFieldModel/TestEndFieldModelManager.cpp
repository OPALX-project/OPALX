//
// Unit tests for class EndFieldModelManager
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

#include "AbsBeamline/EndFieldModel/EndFieldModelManager.h"
#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "gtest/gtest.h"

TEST(TestEndFieldModelManager, TestEFMSetGet) {
    auto efmMan = endfieldmodel::EndFieldModelManager::getEFMManager();
    EXPECT_ANY_THROW(efmMan->getEndFieldModel("test"));

    auto efm    = std::make_shared<endfieldmodel::Tanh>();
    auto efmBad = std::make_shared<endfieldmodel::Tanh>();
    efmMan->setEndFieldModel("test", efm);
    EXPECT_EQ(efmMan->getEndFieldModel("test"), efm);
    EXPECT_EQ(efmMan->getName(efm), "test");
    EXPECT_ANY_THROW(efmMan->getName(efmBad));

    efmMan->clearEFMManager();

    EXPECT_ANY_THROW(efmMan->getEndFieldModel("test"));
}
