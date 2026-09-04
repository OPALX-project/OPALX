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
#include <mpi.h>

#include "Ippl.h"

#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/MultipoleTFieldModel.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MultipoleTRep.h"
#include "BeamlineGeometry/Geometry.h"
#include "Elements/OpalMultipoleT.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace {

    /// What the deck writes. An unset field is left at its attribute default, which is
    /// what update() distinguishes through Attribute::defaultUsed().
    struct Spec {
        double length = 1.0;
        std::optional<double> angle;
        std::optional<std::vector<double>> profile;
        std::optional<double> lfringe;
        std::optional<double> rfringe;
        std::optional<double> maxForder;
        std::optional<double> hapert;
        std::optional<double> vapert;
        std::optional<std::string> aperture;
        std::optional<std::string> scalingModel;
    };

    std::unique_ptr<OpalMultipoleT> makeMagnet(const Spec& spec) {
        auto magnet = std::make_unique<OpalMultipoleT>();

        Attributes::setReal(*magnet->findAttribute("L"), spec.length);
        Attributes::setReal(*magnet->findAttribute("ELEMEDGE"), 0.0);
        if (spec.angle.has_value()) {
            Attributes::setReal(*magnet->findAttribute("ANGLE"), spec.angle.value());
        }
        if (spec.profile.has_value()) {
            Attributes::setRealArray(*magnet->findAttribute("TP"), spec.profile.value());
        }
        if (spec.lfringe.has_value()) {
            Attributes::setReal(*magnet->findAttribute("LFRINGE"), spec.lfringe.value());
        }
        if (spec.rfringe.has_value()) {
            Attributes::setReal(*magnet->findAttribute("RFRINGE"), spec.rfringe.value());
        }
        if (spec.maxForder.has_value()) {
            Attributes::setReal(*magnet->findAttribute("MAXFORDER"), spec.maxForder.value());
        }
        if (spec.hapert.has_value()) {
            Attributes::setReal(*magnet->findAttribute("HAPERT"), spec.hapert.value());
        }
        if (spec.vapert.has_value()) {
            Attributes::setReal(*magnet->findAttribute("VAPERT"), spec.vapert.value());
        }
        if (spec.aperture.has_value()) {
            Attributes::setString(*magnet->findAttribute("APERTURE"), spec.aperture.value());
        }
        if (spec.scalingModel.has_value()) {
            Attributes::setUpperCaseString(
                    *magnet->findAttribute("SCALING_MODEL"), spec.scalingModel.value());
        }
        return magnet;
    }

    MultipoleTRep* elementOf(const std::unique_ptr<OpalMultipoleT>& magnet) {
        return dynamic_cast<MultipoleTRep*>(magnet->getElement());
    }

    class TestOpalMultipoleT : public ::testing::Test {
    protected:
        // The element allocates a Kokkos view for its tanh table, so the runtime has to be
        // up before the first element is built.
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }

        static void TearDownTestSuite() { ippl::finalize(); }
    };

}  // namespace

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

// Without a bend angle the body is straight and L is its length.
TEST_F(TestOpalMultipoleT, NoAngleGivesAStraightBody) {
    auto magnet = makeMagnet({.length = 2.5});
    ASSERT_NO_THROW(magnet->update());

    const Geometry& geometry = elementOf(magnet)->getGeometry();
    EXPECT_EQ(geometry.kind(), GeometryKind::Straight);
    EXPECT_FALSE(geometry.isBend());
    EXPECT_DOUBLE_EQ(geometry.getElementLength(), 2.5);
    EXPECT_DOUBLE_EQ(geometry.getCurvature(), 0.0);
}

// A bend angle makes the body a planar arc, with L the arc length: the same
// convention as SBEND, so the placement treats it as a bend.
TEST_F(TestOpalMultipoleT, AngleGivesASectorBody) {
    const double angle = Physics::pi / 6.0;
    auto magnet        = makeMagnet({.length = 1.5, .angle = angle});
    ASSERT_NO_THROW(magnet->update());

    const Geometry& geometry = elementOf(magnet)->getGeometry();
    EXPECT_EQ(geometry.kind(), GeometryKind::SBend);
    EXPECT_TRUE(geometry.isBend());
    EXPECT_FALSE(geometry.isRectangularBend());
    EXPECT_DOUBLE_EQ(geometry.getElementLength(), 1.5);
    EXPECT_DOUBLE_EQ(geometry.getArcLength(), 1.5);
    EXPECT_NEAR(geometry.getBendAngle(), angle, 1.0e-12);
    EXPECT_NEAR(geometry.getCurvature(), angle / 1.5, 1.0e-12);
}

// A negative bend angle bends the other way: the geometry keeps the sign and the
// curvature goes with it, so the field frame follows the design orbit either way.
// (SBEND flips the sign into a 180 degree roll instead, but only once the element
// has been placed; before that both keep the sign the deck gave.)
TEST_F(TestOpalMultipoleT, NegativeAngleBendsTheOtherWay) {
    const double angle = -Physics::pi / 6.0;
    auto magnet        = makeMagnet({.length = 1.0, .angle = angle});
    ASSERT_NO_THROW(magnet->update());

    const Geometry& geometry = elementOf(magnet)->getGeometry();
    EXPECT_TRUE(geometry.isBend());
    EXPECT_NEAR(geometry.getBendAngle(), angle, 1.0e-12);
    EXPECT_NEAR(geometry.getCurvature(), angle, 1.0e-12);
}

// ---------------------------------------------------------------------------
// The field profile
// ---------------------------------------------------------------------------

// TP is taken as the physical mid-plane field, zero-padded to the number of poles
// the field model carries.
TEST_F(TestOpalMultipoleT, TransverseProfileIsZeroPadded) {
    auto magnet = makeMagnet({.profile = std::vector<double>{-0.35, 1.25}});
    ASSERT_NO_THROW(magnet->update());

    const auto& profile = elementOf(magnet)->getTransverseProfile();
    EXPECT_DOUBLE_EQ(profile[0], -0.35);
    EXPECT_DOUBLE_EQ(profile[1], 1.25);
    for (unsigned int i = 2; i < MultipoleTFieldModel::NumPoles; ++i) {
        EXPECT_DOUBLE_EQ(profile[i], 0.0) << "coefficient " << i;
    }
    // TP[0] is the dipole coefficient, which is also the "BY" channel.
    EXPECT_DOUBLE_EQ(elementOf(magnet)->getB(), -0.35);
}

// More coefficients than there are poles is a deck error.
TEST_F(TestOpalMultipoleT, TooManyProfileCoefficientsThrows) {
    auto magnet = makeMagnet(
            {.profile = std::vector<double>(MultipoleTFieldModel::NumPoles + 1, 1.0)});
    EXPECT_THROW(magnet->update(), OpalException);
}

// The fringe lengths are stored as given; zero is a hard edge.
TEST_F(TestOpalMultipoleT, FringeLengths) {
    auto magnet = makeMagnet({.lfringe = 0.05, .rfringe = 0.07});
    ASSERT_NO_THROW(magnet->update());
    EXPECT_DOUBLE_EQ(elementOf(magnet)->getFringeLeft(), 0.05);
    EXPECT_DOUBLE_EQ(elementOf(magnet)->getFringeRight(), 0.07);

    auto hardEdge = makeMagnet({.length = 1.0});
    ASSERT_NO_THROW(hardEdge->update());
    double zBegin = 0.0;
    double zEnd   = 0.0;
    elementOf(hardEdge)->getFieldExtent(zBegin, zEnd);
    EXPECT_DOUBLE_EQ(zBegin, 0.0);
    EXPECT_DOUBLE_EQ(zEnd, 1.0);
}

// A negative fringe length has no meaning.
TEST_F(TestOpalMultipoleT, NegativeFringeLengthThrows) {
    EXPECT_THROW(makeMagnet({.lfringe = -0.01})->update(), OpalException);
    EXPECT_THROW(makeMagnet({.rfringe = -0.01})->update(), OpalException);
}

// The expansion order has to stay inside what the derivative arrays hold.
TEST_F(TestOpalMultipoleT, MaxFOrderRange) {
    EXPECT_THROW(makeMagnet({.maxForder = 0.0})->update(), OpalException);
    EXPECT_THROW(
            makeMagnet({.maxForder = MultipoleTFieldModel::MaxFOrder + 1.0})->update(),
            OpalException);

    auto lowest = makeMagnet({.maxForder = 1.0});
    ASSERT_NO_THROW(lowest->update());
    EXPECT_EQ(elementOf(lowest)->getMaxFOrder(), 1u);

    auto highest = makeMagnet({.maxForder = MultipoleTFieldModel::MaxFOrder});
    ASSERT_NO_THROW(highest->update());
    EXPECT_EQ(elementOf(highest)->getMaxFOrder(), MultipoleTFieldModel::MaxFOrder);

    // The default is used when the deck says nothing.
    auto plain = makeMagnet({});
    ASSERT_NO_THROW(plain->update());
    EXPECT_EQ(elementOf(plain)->getMaxFOrder(), 3u);
}

// ---------------------------------------------------------------------------
// Aperture
// ---------------------------------------------------------------------------

// Without any aperture attribute the element keeps the default.
TEST_F(TestOpalMultipoleT, NoApertureAttributesKeepsTheDefault) {
    auto magnet = makeMagnet({});
    ASSERT_NO_THROW(magnet->update());

    const auto aperture = elementOf(magnet)->getAperture();
    EXPECT_EQ(aperture.first, ApertureType::ELLIPTICAL);
    EXPECT_DOUBLE_EQ(aperture.second[0], 1e6);
    EXPECT_DOUBLE_EQ(aperture.second[1], 1e6);
}

// HAPERT and VAPERT are the full width and height of a rectangular chamber; the
// element stores half widths.
TEST_F(TestOpalMultipoleT, HapertAndVapertInstallARectangle) {
    auto magnet = makeMagnet({.hapert = 0.3, .vapert = 0.1});
    ASSERT_NO_THROW(magnet->update());

    const auto aperture = elementOf(magnet)->getAperture();
    EXPECT_EQ(aperture.first, ApertureType::RECTANGULAR);
    EXPECT_DOUBLE_EQ(aperture.second[0], 0.15);
    EXPECT_DOUBLE_EQ(aperture.second[1], 0.05);
}

// The generic APERTURE keeps its own shape, and its widths are halved the same way.
TEST_F(TestOpalMultipoleT, ApertureKeepsItsOwnShape) {
    auto magnet = makeMagnet({.aperture = std::string("ELLIPSE(0.3, 0.1)")});
    ASSERT_NO_THROW(magnet->update());

    const auto aperture = elementOf(magnet)->getAperture();
    EXPECT_EQ(aperture.first, ApertureType::ELLIPTICAL);
    EXPECT_DOUBLE_EQ(aperture.second[0], 0.15);
    EXPECT_DOUBLE_EQ(aperture.second[1], 0.05);
}

// The two ways of giving an aperture are alternatives, and HAPERT/VAPERT go together.
TEST_F(TestOpalMultipoleT, ConflictingApertureAttributesThrow) {
    EXPECT_THROW(
            makeMagnet({.hapert = 0.3, .aperture = std::string("ELLIPSE(0.3, 0.1)")})->update(),
            OpalException);
    EXPECT_THROW(
            makeMagnet({.vapert = 0.1, .aperture = std::string("ELLIPSE(0.3, 0.1)")})->update(),
            OpalException);
    EXPECT_THROW(makeMagnet({.hapert = 0.3})->update(), OpalException);
    EXPECT_THROW(makeMagnet({.vapert = 0.1})->update(), OpalException);
    EXPECT_THROW(makeMagnet({.hapert = 0.0, .vapert = 0.1})->update(), OpalException);
    EXPECT_THROW(makeMagnet({.hapert = 0.3, .vapert = -0.1})->update(), OpalException);
}

// ---------------------------------------------------------------------------
// The scaling model
// ---------------------------------------------------------------------------

// The name of the time dependence is stored upper case, like every element name.
TEST_F(TestOpalMultipoleT, ScalingModelNameIsUpperCased) {
    auto magnet = makeMagnet({.scalingModel = std::string("scaling")});
    ASSERT_NO_THROW(magnet->update());
    EXPECT_EQ(elementOf(magnet)->getScalingName(), "SCALING");

    auto plain = makeMagnet({});
    ASSERT_NO_THROW(plain->update());
    EXPECT_EQ(elementOf(plain)->getScalingName(), "");
}

// ---------------------------------------------------------------------------
// Attributes and housekeeping
// ---------------------------------------------------------------------------

// The element carries only the attributes the field needs plus the generic ones.
// The old field-frame attributes are gone: a rolled magnet uses PSI, and the
// aperture window is the aperture.
TEST_F(TestOpalMultipoleT, DroppedAttributesAreNoLongerDeclared) {
    OpalMultipoleT magnet;
    for (const char* name :
         {"ROTATION", "EANGLE", "BBLENGTH", "MAXXORDER", "VARRADIUS", "ENTRYOFFSET"}) {
        EXPECT_EQ(magnet.findAttribute(name), nullptr) << name << " should be gone";
    }
    for (const char* name :
         {"TP", "LFRINGE", "RFRINGE", "MAXFORDER", "ANGLE", "HAPERT", "VAPERT", "SCALING_MODEL",
          "L", "ELEMEDGE", "APERTURE", "PSI"}) {
        EXPECT_NE(magnet.findAttribute(name), nullptr) << name << " should be declared";
    }
}

// What the element reports about itself.
TEST_F(TestOpalMultipoleT, TypeAndPrint) {
    auto magnet = makeMagnet({});
    ASSERT_NO_THROW(magnet->update());
    EXPECT_EQ(elementOf(magnet)->getType(), ElementType::MULTIPOLET);
    EXPECT_EQ(elementOf(magnet)->getTypeString(), "MULTIPOLET");

    std::stringstream stream;
    magnet->print(stream);
    EXPECT_EQ(stream.str(), "MULTIPOLET;\n");
}

// Features that are not wired up are rejected rather than ignored.
TEST_F(TestOpalMultipoleT, UnsupportedFeaturesThrow) {
    auto wake = makeMagnet({});
    Attributes::setString(*wake->findAttribute("WAKEF"), "SOMEWAKE");
    EXPECT_THROW(wake->update(), OpalException);

    auto interaction = makeMagnet({});
    Attributes::setString(
            *interaction->findAttribute("PARTICLEMATTERINTERACTION"), "SOMEINTERACTION");
    EXPECT_THROW(interaction->update(), OpalException);
}

// A cloned element is a separate magnet with its own state; it has to be placed,
// like every other non-builtin element.
TEST_F(TestOpalMultipoleT, Clone) {
    OpalMultipoleT exemplar;
    std::unique_ptr<OpalMultipoleT> clone(exemplar.clone("MYMAGNET"));
    ASSERT_NE(clone, nullptr);

    Attributes::setReal(*clone->findAttribute("L"), 0.8);
    Attributes::setReal(*clone->findAttribute("ELEMEDGE"), 0.0);
    Attributes::setReal(*clone->findAttribute("ANGLE"), 0.25);
    Attributes::setRealArray(*clone->findAttribute("TP"), std::vector<double>{-0.4});
    Attributes::setReal(*clone->findAttribute("LFRINGE"), 0.02);
    Attributes::setReal(*clone->findAttribute("RFRINGE"), 0.02);
    ASSERT_NO_THROW(clone->update());

    auto* element = dynamic_cast<MultipoleTRep*>(clone->getElement());
    ASSERT_NE(element, nullptr);
    EXPECT_DOUBLE_EQ(element->getTransverseProfile()[0], -0.4);
    EXPECT_NEAR(element->getGeometry().getBendAngle(), 0.25, 1.0e-12);
    EXPECT_DOUBLE_EQ(element->getGeometry().getElementLength(), 0.8);
}
