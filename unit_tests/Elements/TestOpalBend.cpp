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
#include "Attributes/Attributes.h"
#include "Elements/OpalElement.h"
#include "Elements/OpalRBend.h"
#include "Elements/OpalSBend.h"
#include "Utilities/OpalException.h"
#include "ValueDefinitions/RealVariable.h"
#include "gtest/gtest.h"

#include <memory>
#include <optional>
#include <string>
#include <vector>

// The transverse aperture of a bend is bounded vertically by the pole gap (HGAP) and
// horizontally by HAPERT or by APERTURE. HAPERT and APERTURE are mutually exclusive, and
// neither works without a positive HGAP.

namespace {

    enum class BendKind { SBend, RBend };

    /// What the deck writes. An unset field is left at its attribute default, which is
    /// what update() distinguishes through Attribute::defaultUsed().
    struct BendSpec {
        std::optional<double> hgap;
        std::optional<double> hapert;
        std::optional<std::string> aperture;
    };

    std::unique_ptr<OpalElement> makeBend(BendKind kind, const BendSpec& spec) {
        std::unique_ptr<OpalElement> bend;
        if (kind == BendKind::SBend) {
            bend = std::make_unique<OpalSBend>();
        } else {
            bend = std::make_unique<OpalRBend>();
        }

        Attributes::setReal(*bend->findAttribute("L"), 1.0);
        Attributes::setReal(*bend->findAttribute("ELEMEDGE"), 0.0);
        if (spec.hgap.has_value()) {
            Attributes::setReal(*bend->findAttribute("HGAP"), spec.hgap.value());
        }
        if (spec.hapert.has_value()) {
            Attributes::setReal(*bend->findAttribute("HAPERT"), spec.hapert.value());
        }
        if (spec.aperture.has_value()) {
            Attributes::setString(*bend->findAttribute("APERTURE"), spec.aperture.value());
        }
        bend->update();
        return bend;
    }

    void expectAperture(
            const std::unique_ptr<OpalElement>& bend, ApertureType type,
            const std::vector<double>& args) {
        ElementBase* element = bend->getElement();
        ASSERT_NE(element, nullptr);
        const auto aperture = element->getAperture();
        EXPECT_EQ(aperture.first, type);
        ASSERT_EQ(aperture.second.size(), args.size());
        for (size_t i = 0; i < args.size(); ++i) {
            EXPECT_DOUBLE_EQ(aperture.second[i], args[i]) << "aperture arg " << i;
        }
    }

    class TestOpalBend : public testing::TestWithParam<BendKind> {
    protected:
        // The bends allocate Kokkos views for their multipole coefficients, so the runtime
        // has to be up before the first update().
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }

        static void TearDownTestSuite() { ippl::finalize(); }

        void SetUp() override {
            // Both bends read OpalData::getP0() when scaling their multipole coefficients.
            // In a full run the RealVariable exemplar creates and registers the P0 variable;
            // without it getP0() dereferences a null pointer. Created once for the suite and
            // intentionally never deleted, since OpalData keeps referring to it.
            [[maybe_unused]] static RealVariable* p0Owner = new RealVariable();
        }
    };

}  // namespace

TEST_P(TestOpalBend, NoApertureAttributesKeepsTheDefault) {
    expectAperture(makeBend(GetParam(), {}), ApertureType::ELLIPTICAL, {1e6, 1e6});
}

TEST_P(TestOpalBend, HgapAloneInstallsTheVerticalWall) {
    // The poles are the vertical limit even when the deck asks for nothing else; the
    // horizontal half-width stays at the "no limit" default.
    BendSpec spec;
    spec.hgap = 0.02;
    expectAperture(makeBend(GetParam(), spec), ApertureType::RECTANGULAR, {1e6, 0.02});
}

TEST_P(TestOpalBend, HapertSetsTheHorizontalHalfWidthOnly) {
    BendSpec spec;
    spec.hgap   = 0.02;
    spec.hapert = 0.15;
    // Not {0.15, 0.15}: HAPERT is the bend-plane half-aperture, not a square.
    expectAperture(makeBend(GetParam(), spec), ApertureType::RECTANGULAR, {0.15, 0.02});
}

TEST_P(TestOpalBend, HapertIsAHalfWidthWhileApertureTakesFullWidths) {
    BendSpec fromHapert;
    fromHapert.hgap   = 0.02;
    fromHapert.hapert = 0.15;

    BendSpec fromAperture;
    fromAperture.hgap     = 0.02;
    fromAperture.aperture = "RECTANGLE(0.3,0.04)";

    // 0.15 written as HAPERT is the same wall as 0.3 written as an APERTURE full width.
    expectAperture(makeBend(GetParam(), fromHapert), ApertureType::RECTANGULAR, {0.15, 0.02});
    expectAperture(makeBend(GetParam(), fromAperture), ApertureType::RECTANGULAR, {0.15, 0.02});
}

TEST_P(TestOpalBend, ApertureKeepsItsOwnShape) {
    BendSpec spec;
    spec.hgap     = 0.02;
    spec.aperture = "ELLIPSE(0.3,0.03)";
    expectAperture(makeBend(GetParam(), spec), ApertureType::ELLIPTICAL, {0.15, 0.015});
}

TEST_P(TestOpalBend, ApertureFlushWithThePoleFacesIsAccepted) {
    // Vertical half-aperture exactly HGAP: a chamber flush with the poles is legal.
    BendSpec spec;
    spec.hgap     = 0.02;
    spec.aperture = "RECTANGLE(0.3,0.04)";
    EXPECT_NO_THROW(makeBend(GetParam(), spec));
}

TEST_P(TestOpalBend, ApertureTallerThanTheGapThrows) {
    BendSpec spec;
    spec.hgap     = 0.02;
    spec.aperture = "RECTANGLE(0.3,0.05)";  // vertical half-width 0.025 > HGAP
    EXPECT_THROW(makeBend(GetParam(), spec), OpalException);
}

TEST_P(TestOpalBend, ApertureAndHapertTogetherThrow) {
    BendSpec spec;
    spec.hgap     = 0.02;
    spec.hapert   = 0.15;
    spec.aperture = "RECTANGLE(0.3,0.04)";
    EXPECT_THROW(makeBend(GetParam(), spec), OpalException);
}

TEST_P(TestOpalBend, ScrapingWithoutAPositiveHgapThrows) {
    BendSpec hapertNoHgap;
    hapertNoHgap.hapert = 0.15;
    EXPECT_THROW(makeBend(GetParam(), hapertNoHgap), OpalException);

    BendSpec hapertZeroHgap;
    hapertZeroHgap.hgap   = 0.0;
    hapertZeroHgap.hapert = 0.15;
    EXPECT_THROW(makeBend(GetParam(), hapertZeroHgap), OpalException);

    BendSpec apertureNoHgap;
    apertureNoHgap.aperture = "RECTANGLE(0.3,0.04)";
    EXPECT_THROW(makeBend(GetParam(), apertureNoHgap), OpalException);

    BendSpec apertureZeroHgap;
    apertureZeroHgap.hgap     = 0.0;
    apertureZeroHgap.aperture = "RECTANGLE(0.3,0.04)";
    EXPECT_THROW(makeBend(GetParam(), apertureZeroHgap), OpalException);
}

TEST_P(TestOpalBend, ExplicitZeroHapertThrows) {
    // A zero half-aperture would put every particle outside the element.
    BendSpec spec;
    spec.hgap   = 0.02;
    spec.hapert = 0.0;
    EXPECT_THROW(makeBend(GetParam(), spec), OpalException);
}

TEST_P(TestOpalBend, HardEdgeWithoutScrapingStillWorks) {
    // HGAP = 0 means "no fringe modelled". It is only an error together with HAPERT or
    // APERTURE; on its own it must leave the bend unrestricted.
    BendSpec spec;
    spec.hgap = 0.0;
    expectAperture(makeBend(GetParam(), spec), ApertureType::ELLIPTICAL, {1e6, 1e6});
}

INSTANTIATE_TEST_SUITE_P(
        Bends, TestOpalBend, testing::Values(BendKind::SBend, BendKind::RBend),
        [](const testing::TestParamInfo<BendKind>& info) {
            return info.param == BendKind::SBend ? "SBEND" : "RBEND";
        });
