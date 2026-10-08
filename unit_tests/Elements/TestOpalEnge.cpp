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
#include "Elements/OpalEnge.h"
#include "AbsBeamline/EndFieldModel/Enge.h"


class TestOpalEnge : public testing::Test {
public:
    TestOpalEnge() {
    }

    static void SetUpTestSuite() {
        Kokkos::initialize();
    }
    static void TearDownTestSuite() {
        endfieldmodel::EndFieldModelManager::clearEFMManager();
        Kokkos::fence();
        Kokkos::finalize();
    }
};


TEST_F(TestOpalEnge, TestSetup) {
    // Make the UI
    OpalEnge ui;
    // Set the attributes
    Attributes::setReal(ui.itsAttr[OpalEnge::LENGTH], 4);
    Attributes::setReal(ui.itsAttr[OpalEnge::FRINGE_LENGTH], 2);
    Attributes::setRealArray(ui.itsAttr[OpalEnge::COEFFICIENTS], {0.0, 1.0, 2.0});
    ui.update();
    auto efmMan = endfieldmodel::EndFieldModelManager::getEFMManager();
    EXPECT_NO_THROW(efmMan->getEndFieldModel("ENGE"));
    auto efm = efmMan->getEndFieldModel("ENGE");
    auto enge = std::dynamic_pointer_cast<endfieldmodel::Enge, endfieldmodel::EndFieldModel>(efm);
    EXPECT_EQ(efm->getCentreLength(), 4.0);
    EXPECT_EQ(efm->getEndLength(), 2.0);
    EXPECT_EQ(enge->getCoefficients(), std::vector<double>({0.0, 1.0, 2.0}));
}

