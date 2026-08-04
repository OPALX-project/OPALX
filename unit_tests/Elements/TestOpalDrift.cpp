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
#include "Elements/OpalDrift.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace {

    class TestOpalDrift : public testing::Test {
    public:
        TestOpalDrift() = default;
    };

    std::unique_ptr<OpalDrift> makeDrift(const std::optional<std::string>& aperture) {
        auto drift = std::make_unique<OpalDrift>();
        Attributes::setReal(*drift->findAttribute("L"), 2.0);
        Attributes::setReal(*drift->findAttribute("ELEMEDGE"), 0.0);
        if (aperture.has_value()) {
            Attributes::setString(*drift->findAttribute("APERTURE"), aperture.value());
        }
        drift->update();
        return drift;
    }

}  // namespace

TEST_F(TestOpalDrift, CircleDefaultsMatchDefaultApertureBehaviour) {
    auto noAperture = makeDrift(std::nullopt);
    auto circle     = makeDrift("CIRCLE(1)");

    ElementBase* noApertureElement = noAperture->getElement();
    ElementBase* circleElement     = circle->getElement();

    ASSERT_NE(noApertureElement, nullptr);
    ASSERT_NE(circleElement, nullptr);

    const std::vector<Vector_t<double, 3>> probes = {
            {0.00, 0.00, 0.20}, {0.20, 0.10, 1.00},  {0.49, 0.00, 0.20}, {0.49, 0.00, 1.00},
            {0.49, 0.00, 1.80}, {0.50, 0.00, 1.00},  {0.00, 0.50, 1.00}, {0.36, 0.36, 1.00},
            {0.40, 0.40, 1.00}, {0.00, 0.00, -0.10}, {0.00, 0.00, 2.00}};

    for (const Vector_t<double, 3>& r : probes) {
        SCOPED_TRACE(
                ::testing::Message()
                << "Probe point (" << r[0] << ", " << r[1] << ", " << r[2] << ")");
        EXPECT_EQ(noApertureElement->isInside(r), circleElement->isInside(r));
    }
}

TEST_F(TestOpalDrift, SquareAndRectangleEquivalentBehaviour) {
    auto rectangle = makeDrift("RECTANGLE(1,1)");
    auto square    = makeDrift("SQUARE(1)");

    ElementBase* rectangleElement = rectangle->getElement();
    ElementBase* squareElement    = square->getElement();

    ASSERT_NE(rectangleElement, nullptr);
    ASSERT_NE(squareElement, nullptr);

    const std::vector<Vector_t<double, 3>> probes = {
            {0.00, 0.00, 0.20},  {0.49, 0.49, 1.00}, {0.49, 0.10, 1.80}, {0.10, 0.49, 0.20},
            {0.50, 0.10, 1.00},  {0.10, 0.50, 1.00}, {0.60, 0.40, 1.00}, {0.40, 0.60, 1.00},
            {0.00, 0.00, -0.10}, {0.00, 0.00, 2.00}};

    for (const Vector_t<double, 3>& r : probes) {
        SCOPED_TRACE(
                ::testing::Message()
                << "Probe point (" << r[0] << ", " << r[1] << ", " << r[2] << ")");
        EXPECT_EQ(rectangleElement->isInside(r), squareElement->isInside(r));
    }
}

TEST_F(TestOpalDrift, CircleConstantAlongLengthAndLongitudinalBounds) {
    auto circle                = makeDrift("CIRCLE(1)");
    ElementBase* circleElement = circle->getElement();

    ASSERT_NE(circleElement, nullptr);

    const std::vector<Vector_t<double, 3>> insideTransverseInsideLength = {
            {0.49, 0.00, 0.10}, {0.49, 0.00, 1.00}, {0.49, 0.00, 1.90}};
    for (const Vector_t<double, 3>& r : insideTransverseInsideLength) {
        EXPECT_TRUE(circleElement->isInside(r));
    }

    const std::vector<Vector_t<double, 3>> outsideTransverseInsideLength = {
            {0.51, 0.00, 0.10}, {0.51, 0.00, 1.00}, {0.51, 0.00, 1.90}};
    for (const Vector_t<double, 3>& r : outsideTransverseInsideLength) {
        EXPECT_FALSE(circleElement->isInside(r));
    }

    EXPECT_FALSE(circleElement->isInside(Vector_t<double, 3>({0.00, 0.00, -0.10})));
    EXPECT_FALSE(circleElement->isInside(Vector_t<double, 3>({0.00, 0.00, 2.00})));
}

TEST_F(TestOpalDrift, ConicApertureStringsThrow) {
    // Conic (tapered) apertures are disabled: any extra scale argument is rejected.
    EXPECT_THROW(makeDrift("CIRCLE(1,0.5)"), OpalException);
    EXPECT_THROW(makeDrift("CIRCLE(1,2)"), OpalException);
    EXPECT_THROW(makeDrift("SQUARE(1,1)"), OpalException);
    EXPECT_THROW(makeDrift("ELLIPSE(2,1,2)"), OpalException);
    EXPECT_THROW(makeDrift("RECTANGLE(1,1,-1)"), OpalException);
    EXPECT_THROW(makeDrift("ELLIPSE(1,1,nan)"), OpalException);
}

TEST_F(TestOpalDrift, MalformedRectangleArgumentsThrow) {
    EXPECT_THROW(makeDrift("RECTANGLE(1,)"), OpalException);
    EXPECT_THROW(makeDrift("RECTANGLE(,1)"), OpalException);
}

namespace {

    void expectAperture(
            const std::unique_ptr<OpalDrift>& drift, ApertureType type,
            const std::vector<double>& args) {
        ElementBase* element = drift->getElement();
        ASSERT_NE(element, nullptr);
        const auto aperture = element->getAperture();
        EXPECT_EQ(aperture.first, type);
        ASSERT_EQ(aperture.second.size(), args.size());
        for (size_t i = 0; i < args.size(); ++i) {
            EXPECT_DOUBLE_EQ(aperture.second[i], args[i]) << "aperture arg " << i;
        }
    }

}  // namespace

TEST_F(TestOpalDrift, GetApertureStoresHalfWidths) {
    // The APERTURE grammar takes full widths/diameters; the element stores
    // half-apertures (factor 0.5).
    expectAperture(makeDrift("CIRCLE(1)"), ApertureType::ELLIPTICAL, {0.5, 0.5});
    expectAperture(makeDrift("ELLIPSE(2,1)"), ApertureType::ELLIPTICAL, {1.0, 0.5});
    expectAperture(makeDrift("RECTANGLE(2,4)"), ApertureType::RECTANGULAR, {1.0, 2.0});
    expectAperture(makeDrift("SQUARE(3)"), ApertureType::RECTANGULAR, {1.5, 1.5});
}

TEST_F(TestOpalDrift, ApertureKeywordIsCaseInsensitive) {
    expectAperture(makeDrift("circle(1)"), ApertureType::ELLIPTICAL, {0.5, 0.5});
    expectAperture(makeDrift("Circle (1)"), ApertureType::ELLIPTICAL, {0.5, 0.5});
    expectAperture(makeDrift("rectangle(2,4)"), ApertureType::RECTANGULAR, {1.0, 2.0});
}

TEST_F(TestOpalDrift, UnknownApertureStringThrows) {
    EXPECT_THROW(makeDrift("TRIANGLE(1)"), OpalException);
    EXPECT_THROW(makeDrift("CIRCLE"), OpalException);
}

TEST_F(TestOpalDrift, NoApertureDefaultsToHalfMeterEllipse) {
    expectAperture(makeDrift(std::nullopt), ApertureType::ELLIPTICAL, {0.5, 0.5});
}
