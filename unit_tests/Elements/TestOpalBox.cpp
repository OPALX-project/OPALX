//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#include "AbsBeamline/Box.h"
#include "Attributes/Attributes.h"
#include "Elements/OpalBox.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

#include <memory>
#include <optional>
#include <string>
#include "../KeywordObjects.h"
#include "AbstractObjects/OpalData.h"

namespace {

    class TestOpalBox : public testing::Test {
    public:
        TestOpalBox() = default;

        /**
         * @brief A BOX as an input file would define it.
         *
         * Made from the keyword table, as the parser does for "B1: BOX;". Returns the runtime
         * element a RUN builds for it.
         */
        std::shared_ptr<ElementBase> makeBox(
                double length, double width, double height,
                const std::optional<std::string>& aperture = std::nullopt,
                const std::optional<double>& reach         = std::nullopt) {
            auto box = opalx::test::makeKeywordObject<OpalBox>("BOX", "B1");
            Attributes::setReal(*box->findAttribute("L"), length);
            Attributes::setReal(*box->findAttribute("WIDTH"), width);
            Attributes::setReal(*box->findAttribute("HEIGHT"), height);
            Attributes::setReal(*box->findAttribute("ELEMEDGE"), 0.5);
            if (aperture.has_value()) {
                Attributes::setString(*box->findAttribute("APERTURE"), aperture.value());
            }
            if (reach.has_value()) {
                Attributes::setReal(*box->findAttribute("REACH"), reach.value());
            }
            return opalx::test::buildFrom(*box);
        }
    };

}  // namespace

TEST_F(TestOpalBox, UpdateSetsTypeLengthAndHalfSize) {
    auto box             = makeBox(0.1, 0.02, 0.04);
    ElementBase* element = box.get();
    ASSERT_NE(element, nullptr);
    EXPECT_EQ(element->getType(), ElementType::BOX);
    EXPECT_DOUBLE_EQ(element->getGeometry().getElementLength(), 0.1);

    // Input values are full widths; the element stores half sizes.
    const Box* boxElement = dynamic_cast<const Box*>(element);
    ASSERT_NE(boxElement, nullptr);
    EXPECT_DOUBLE_EQ(boxElement->getHalfWidth(), 0.01);
    EXPECT_DOUBLE_EQ(boxElement->getHalfHeight(), 0.02);

    // The reach is the larger half-size, so the element's box in the grid covers the block.
    ASSERT_TRUE(element->getReach().has_value());
    EXPECT_DOUBLE_EQ(*element->getReach(), 0.02);
}

TEST_F(TestOpalBox, NonPositiveSizeThrows) {
    EXPECT_THROW(makeBox(0.0, 0.02, 0.04), OpalException);
    EXPECT_THROW(makeBox(0.1, 0.0, 0.04), OpalException);
    EXPECT_THROW(makeBox(0.1, 0.02, -0.04), OpalException);
}

TEST_F(TestOpalBox, ApertureOrReachThrows) {
    // WIDTH and HEIGHT set the size of a BOX.
    EXPECT_THROW(makeBox(0.1, 0.02, 0.04, "CIRCLE(0.02)"), OpalException);
    EXPECT_THROW(makeBox(0.1, 0.02, 0.04, std::nullopt, 1.0), OpalException);
}

TEST_F(TestOpalBox, KeywordIsNotUpdatedWithTheInputObjects) {
    // The BOX keyword has none of the attributes an input object must set. It is kept in the
    // keyword table, apart from the input objects, so updating all input objects leaves it alone.
    OpalData* opal = OpalData::getInstance();
    opalx::test::keyword<OpalBox>("BOX");
    opal->makeDirty(nullptr);
    EXPECT_NO_THROW(opal->update());
}
