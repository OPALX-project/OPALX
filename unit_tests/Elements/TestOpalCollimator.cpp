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
#include "AbsBeamline/ElementBase.h"
#include "Attributes/Attributes.h"
#include "Elements/OpalCollimator.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

#include <memory>
#include <optional>
#include <string>

namespace {

    class TestOpalCollimator : public testing::Test {
    public:
        TestOpalCollimator() = default;

        /**
         * @brief A COLLIMATOR as a input would define it.
         *
         * The tests must run against a *clone*, not the exemplar: OpalData
         * calls update() on the builtin prototype too, and the attribute
         * validation is skipped there (same rule as
         * OpalElement::validatePlacement).
         */
        std::unique_ptr<OpalCollimator> makeCollimator(
                double length, const std::optional<std::string>& aperture,
                const std::optional<bool>& deleteOnTransverseExit = std::nullopt) {
            std::unique_ptr<OpalCollimator> coll(exemplar_m.clone("C1"));
            Attributes::setReal(*coll->findAttribute("L"), length);
            Attributes::setReal(*coll->findAttribute("ELEMEDGE"), 0.5);
            if (aperture.has_value()) {
                Attributes::setString(*coll->findAttribute("APERTURE"), aperture.value());
            }
            if (deleteOnTransverseExit.has_value()) {
                Attributes::setBool(
                        *coll->findAttribute("DELETEONTRANSVERSEEXIT"),
                        deleteOnTransverseExit.value());
            }
            coll->update();
            return coll;
        }

    private:
        /// Outlives every clone made from it.
        OpalCollimator exemplar_m;
    };

}  // namespace

TEST_F(TestOpalCollimator, ElementTypeIsCollimator) {
    auto coll            = makeCollimator(0.1, "CIRCLE(0.02)");
    ElementBase* element = coll->getElement();
    ASSERT_NE(element, nullptr);
    EXPECT_EQ(element->getType(), ElementType::COLLIMATOR);
}

TEST_F(TestOpalCollimator, UpdateSetsGeometryLengthAndAperture) {
    auto coll            = makeCollimator(0.1, "CIRCLE(0.02)");
    ElementBase* element = coll->getElement();
    ASSERT_NE(element, nullptr);

    // length is taken as-is (no minimum clamp like MONITOR applies)
    EXPECT_DOUBLE_EQ(element->getGeometry().getElementLength(), 0.1);

    const auto aperture = element->getAperture();
    EXPECT_EQ(aperture.first, ApertureType::ELLIPTICAL);
    ASSERT_EQ(aperture.second.size(), 2u);
    EXPECT_DOUBLE_EQ(aperture.second[0], 0.01);  // diameter -> half-aperture
    EXPECT_DOUBLE_EQ(aperture.second[1], 0.01);
}

TEST_F(TestOpalCollimator, RectangularApertureStringReachesElement) {
    auto coll            = makeCollimator(0.05, "RECTANGLE(0.02,0.04)");
    ElementBase* element = coll->getElement();
    ASSERT_NE(element, nullptr);

    const auto aperture = element->getAperture();
    EXPECT_EQ(aperture.first, ApertureType::RECTANGULAR);
    ASSERT_EQ(aperture.second.size(), 2u);
    EXPECT_DOUBLE_EQ(aperture.second[0], 0.01);
    EXPECT_DOUBLE_EQ(aperture.second[1], 0.02);
}

TEST_F(TestOpalCollimator, UpdateSetsDeleteOnTransverseExitDefaultTrue) {
    auto withDefault = makeCollimator(0.1, "CIRCLE(0.02)");
    ASSERT_NE(withDefault->getElement(), nullptr);
    EXPECT_TRUE(withDefault->getElement()->getFlagDeleteOnTransverseExit());

    auto disabled = makeCollimator(0.1, "CIRCLE(0.02)", false);
    ASSERT_NE(disabled->getElement(), nullptr);
    EXPECT_FALSE(disabled->getElement()->getFlagDeleteOnTransverseExit());
}

TEST_F(TestOpalCollimator, MissingApertureThrows) {
    // Without an aperture the element would inherit the wide-open default and
    // transmit everything, which is never what a collimator is for.
    EXPECT_THROW(makeCollimator(0.1, std::nullopt), OpalException);
    EXPECT_THROW(makeCollimator(0.1, ""), OpalException);
}

TEST_F(TestOpalCollimator, BuiltinExemplarPassesValidation) {
    // OpalData::update() runs update() on the prototype, whose attributes are
    // all unset. Validating it would make every inputfile containing a COLLIMATOR
    // abort during parsing.
    OpalCollimator exemplar;
    EXPECT_NO_THROW(exemplar.update());
}

TEST_F(TestOpalCollimator, NonPositiveLengthThrows) {
    // Scraping only happens for z in [0, L).
    EXPECT_THROW(makeCollimator(0.0, "CIRCLE(0.02)"), OpalException);
    EXPECT_THROW(makeCollimator(-0.1, "CIRCLE(0.02)"), OpalException);
}
