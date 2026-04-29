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
    auto noAperture  = makeDrift(std::nullopt);
    auto circle      = makeDrift("CIRCLE(1)");
    auto conicCircle = makeDrift("CIRCLE(1,1)");

    ElementBase* noApertureElement  = noAperture->getElement();
    ElementBase* circleElement      = circle->getElement();
    ElementBase* conicCircleElement = conicCircle->getElement();

    ASSERT_NE(noApertureElement, nullptr);
    ASSERT_NE(circleElement, nullptr);
    ASSERT_NE(conicCircleElement, nullptr);

    const std::vector<Vector_t<double, 3>> probes = {
            {0.00, 0.00, 0.20}, {0.20, 0.10, 1.00},  {0.49, 0.00, 0.20}, {0.49, 0.00, 1.00},
            {0.49, 0.00, 1.80}, {0.50, 0.00, 1.00},  {0.00, 0.50, 1.00}, {0.36, 0.36, 1.00},
            {0.40, 0.40, 1.00}, {0.00, 0.00, -0.10}, {0.00, 0.00, 2.00}};

    for (const Vector_t<double, 3>& r : probes) {
        SCOPED_TRACE(
                ::testing::Message()
                << "Probe point (" << r[0] << ", " << r[1] << ", " << r[2] << ")");
        EXPECT_EQ(noApertureElement->isInside(r), circleElement->isInside(r));
        EXPECT_EQ(noApertureElement->isInside(r), conicCircleElement->isInside(r));
    }
}

TEST_F(TestOpalDrift, SquareAndRectangleEquivalentBehaviour) {
    auto rectangle      = makeDrift("RECTANGLE(1,1)");
    auto conicRectangle = makeDrift("RECTANGLE(1,1,1)");
    auto square         = makeDrift("SQUARE(1)");
    auto conicSquare    = makeDrift("SQUARE(1,1)");

    ElementBase* rectangleElement      = rectangle->getElement();
    ElementBase* conicRectangleElement = conicRectangle->getElement();
    ElementBase* squareElement         = square->getElement();
    ElementBase* conicSquareElement    = conicSquare->getElement();

    ASSERT_NE(rectangleElement, nullptr);
    ASSERT_NE(conicRectangleElement, nullptr);
    ASSERT_NE(squareElement, nullptr);
    ASSERT_NE(conicSquareElement, nullptr);

    const std::vector<Vector_t<double, 3>> probes = {
            {0.00, 0.00, 0.20},  {0.49, 0.49, 1.00}, {0.49, 0.10, 1.80}, {0.10, 0.49, 0.20},
            {0.50, 0.10, 1.00},  {0.10, 0.50, 1.00}, {0.60, 0.40, 1.00}, {0.40, 0.60, 1.00},
            {0.00, 0.00, -0.10}, {0.00, 0.00, 2.00}};

    for (const Vector_t<double, 3>& r : probes) {
        SCOPED_TRACE(
                ::testing::Message()
                << "Probe point (" << r[0] << ", " << r[1] << ", " << r[2] << ")");
        const bool expected = rectangleElement->isInside(r);
        EXPECT_EQ(expected, conicRectangleElement->isInside(r));
        EXPECT_EQ(expected, squareElement->isInside(r));
        EXPECT_EQ(expected, conicSquareElement->isInside(r));
    }
}

TEST_F(TestOpalDrift, ConicCircleOpeningBehaviour) {
    auto conicCircle                = makeDrift("CIRCLE(1,2)");
    ElementBase* conicCircleElement = conicCircle->getElement();

    ASSERT_NE(conicCircleElement, nullptr);

    const Vector_t<double, 3> centerline = {0.00, 0.00, 1.00};
    EXPECT_TRUE(conicCircleElement->isInside(centerline));

    const Vector_t<double, 3> startProbe  = {0.75, 0.00, 0.10};
    const Vector_t<double, 3> middleProbe = {0.75, 0.00, 1.00};
    const Vector_t<double, 3> endProbe    = {0.75, 0.00, 1.90};

    EXPECT_FALSE(conicCircleElement->isInside(startProbe));
    EXPECT_FALSE(conicCircleElement->isInside(middleProbe));
    EXPECT_TRUE(conicCircleElement->isInside(endProbe));
}

TEST_F(TestOpalDrift, ConicCircleClosingBehaviour) {
    auto conicCircle                = makeDrift("CIRCLE(1,0.5)");
    ElementBase* conicCircleElement = conicCircle->getElement();

    ASSERT_NE(conicCircleElement, nullptr);

    const Vector_t<double, 3> startProbe  = {0.40, 0.00, 0.10};
    const Vector_t<double, 3> middleProbe = {0.40, 0.00, 1.00};
    const Vector_t<double, 3> endProbe    = {0.40, 0.00, 1.90};

    EXPECT_TRUE(conicCircleElement->isInside(startProbe));
    EXPECT_FALSE(conicCircleElement->isInside(middleProbe));
    EXPECT_FALSE(conicCircleElement->isInside(endProbe));
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