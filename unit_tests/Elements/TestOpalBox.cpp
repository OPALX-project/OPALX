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

namespace {

    class TestOpalBox : public testing::Test {
    public:
        TestOpalBox() = default;

        /**
         * @brief A BOX as an input file would define it.
         *
         * The tests must run against a *clone*, not the exemplar: OpalData
         * calls update() on the builtin prototype too, and the attribute
         * validation is skipped there.
         */
        std::unique_ptr<OpalBox> makeBox(
                double length, double width, double height,
                const std::optional<std::string>& aperture = std::nullopt) {
            std::unique_ptr<OpalBox> box(exemplar_m.clone("B1"));
            Attributes::setReal(*box->findAttribute("L"), length);
            Attributes::setReal(*box->findAttribute("WIDTH"), width);
            Attributes::setReal(*box->findAttribute("HEIGHT"), height);
            Attributes::setReal(*box->findAttribute("ELEMEDGE"), 0.5);
            if (aperture.has_value()) {
                Attributes::setString(*box->findAttribute("APERTURE"), aperture.value());
            }
            box->update();
            return box;
        }

    private:
        /// Outlives every clone made from it.
        OpalBox exemplar_m;
    };

}  // namespace

TEST_F(TestOpalBox, UpdateSetsTypeLengthAndHalfSize) {
    auto box             = makeBox(0.1, 0.02, 0.04);
    ElementBase* element = box->getElement();
    ASSERT_NE(element, nullptr);
    EXPECT_EQ(element->getType(), ElementType::BOX);
    EXPECT_DOUBLE_EQ(element->getGeometry().getElementLength(), 0.1);

    // Input values are full widths; the element stores half sizes.
    const Box* boxElement = dynamic_cast<const Box*>(element);
    ASSERT_NE(boxElement, nullptr);
    EXPECT_DOUBLE_EQ(boxElement->getHalfWidth(), 0.01);
    EXPECT_DOUBLE_EQ(boxElement->getHalfHeight(), 0.02);
}

TEST_F(TestOpalBox, NonPositiveSizeThrows) {
    EXPECT_THROW(makeBox(0.0, 0.02, 0.04), OpalException);
    EXPECT_THROW(makeBox(0.1, 0.0, 0.04), OpalException);
    EXPECT_THROW(makeBox(0.1, 0.02, -0.04), OpalException);
}

TEST_F(TestOpalBox, ApertureThrows) {
    // The aperture decides whether the element is registered along the
    // reference path, so a BOX must keep the wide-open default.
    EXPECT_THROW(makeBox(0.1, 0.02, 0.04, "CIRCLE(0.02)"), OpalException);
}

TEST_F(TestOpalBox, BuiltinExemplarPassesValidation) {
    OpalBox exemplar;
    EXPECT_NO_THROW(exemplar.update());
}
