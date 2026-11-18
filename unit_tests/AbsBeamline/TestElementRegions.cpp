/**
 * \file TestElementRegions.cpp
 * \brief Unit tests for ElementRegions and ElementBase::getRegions() / checkRegions().
 *
 * An element applies its field inside its field window and aperture, and deletes particles
 * inside its body window and reach but outside its aperture. Beyond the reach it does nothing.
 */

#include <gtest/gtest.h>

#include "Ippl.h"

#include "AbsBeamline/ElementBase.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/MultipoleRep.h"
#include "BeamlineCore/SBendRep.h"
#include "BeamlineGeometry/Geometry.h"
#include "Utilities/GeneralOpalException.h"
#include "Utilities/Options.h"

#include <cmath>
#include <optional>

namespace {
    using Vector = Vector_t<double, 3>;

    /// A drift whose field model covers an off-centre rectangle, like a 3D field map.
    class OffCentreRangeDrift final : public DriftRep {
    public:
        OffCentreRangeDrift(double centreX, double halfX)
            : DriftRep("OFF_CENTRE"), centreX_m(centreX), halfX_m(halfX) {
            getGeometry().setElementLength(1.0);
        }
        ElementBase* clone() const override { return new OffCentreRangeDrift(*this); }
        std::optional<ApertureShape> getFieldRange() const override {
            ApertureShape range;
            range.centreX = centreX_m;
            range.halfX   = halfX_m;
            range.halfY   = 0.05;
            return range;
        }

    private:
        double centreX_m;
        double halfX_m;
    };

    /// A drift with an off-centre field map, like a 3D G4beamline grid.
    class OffCentreMapDrift final : public DriftRep {
    public:
        OffCentreMapDrift() : DriftRep("MAP") { getGeometry().setElementLength(1.0); }
        ElementBase* clone() const override { return new OffCentreMapDrift(*this); }
        std::optional<ApertureShape> getFieldMapRange() const override {
            ApertureShape range;
            range.centreX = 0.5;  // x = 0.2 .. 0.8 m, wider than the 0.5 m reach
            range.halfX   = 0.3;
            range.halfY   = 0.1;
            return range;
        }
    };

    class ElementRegionsTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }
        static void TearDownTestSuite() { ippl::finalize(); }

        void SetUp() override { Options::elementMappingReach = 0.5; }

        static MultipoleRep makeQuadrupole(const double apertureHalfWidth) {
            MultipoleRep quad("Q");
            quad.getGeometry().setElementLength(1.0);
            if (apertureHalfWidth > 0.0)
                quad.setAperture(ApertureType::RECTANGULAR, {apertureHalfWidth, apertureHalfWidth});
            return quad;
        }
    };

    TEST_F(ElementRegionsTest, DefaultAppliesEverywhereAndDeletesNothing) {
        const ElementRegions regions;
        EXPECT_TRUE(regions.inFieldRegion(Vector({100.0, -100.0, 1.0e6})));
        EXPECT_FALSE(regions.inLossRegion(Vector({100.0, -100.0, 0.0})));
    }

    TEST_F(ElementRegionsTest, FieldInsideApertureLossBetweenApertureAndReach) {
        const ElementRegions regions = makeQuadrupole(0.05).getRegions();

        // Field: in [0, L) and inside the aperture.
        EXPECT_TRUE(regions.inFieldRegion(Vector({0.01, 0.0, 0.5})));
        EXPECT_TRUE(regions.inFieldRegion(Vector({0.0, 0.0, 0.0})));
        EXPECT_FALSE(regions.inFieldRegion(Vector({0.0, 0.0, 1.0})));

        // Lost: in the body, inside the reach, outside the aperture.
        EXPECT_FALSE(regions.inLossRegion(Vector({0.01, 0.0, 0.5})));
        EXPECT_TRUE(regions.inLossRegion(Vector({0.1, 0.0, 0.5})));
        EXPECT_TRUE(regions.inLossRegion(Vector({0.0, -0.4, 0.5})));

        // Beyond the reach, or outside the body: nothing.
        EXPECT_FALSE(regions.inFieldRegion(Vector({0.6, 0.0, 0.5})));
        EXPECT_FALSE(regions.inLossRegion(Vector({0.6, 0.0, 0.5})));
        EXPECT_FALSE(regions.inLossRegion(Vector({0.1, 0.0, 1.2})));
    }

    TEST_F(ElementRegionsTest, WithoutApertureTheReachIsTheAperture) {
        MultipoleRep quad            = makeQuadrupole(0.0);
        const ElementRegions regions = quad.getRegions();
        EXPECT_TRUE(regions.inFieldRegion(Vector({0.4, 0.0, 0.5})));
        EXPECT_FALSE(regions.inFieldRegion(Vector({0.6, 0.0, 0.5})));
        EXPECT_FALSE(regions.inLossRegion(Vector({0.4, 0.0, 0.5})));
        EXPECT_FALSE(regions.inLossRegion(Vector({0.6, 0.0, 0.5})));

        quad.setReach(1.0);
        EXPECT_TRUE(quad.getRegions().inFieldRegion(Vector({0.6, 0.0, 0.5})));
    }

    TEST_F(ElementRegionsTest, SectorBendUsesArcCoordinates) {
        SBendRep bend("B");
        constexpr double curvature = 0.5;
        bend.getGeometry()         = Geometry::makeSBend(1.0, curvature);
        bend.getGeometry().setElementLength(1.0);
        bend.getGeometry().setBendAngle(curvature);
        bend.setAperture(ApertureType::RECTANGULAR, {0.05, 0.05});
        const ElementRegions regions = bend.getRegions();

        // A point at arc length s = 0.9 m, offset radially by dx from the design arc.
        const auto onArc = [](const double s, const double dx) {
            const double phi = curvature * s;
            return Vector(
                    {(std::cos(phi) - 1.0) / curvature + dx * std::cos(phi), 0.0,
                     std::sin(phi) / curvature + dx * std::sin(phi)});
        };
        EXPECT_TRUE(regions.inFieldRegion(onArc(0.9, 0.01)));
        EXPECT_TRUE(regions.inLossRegion(onArc(0.9, 0.1)));
        // The straight entrance frame would put this point 0.17 m off axis.
        EXPECT_FALSE(regions.inLossRegion(onArc(0.9, 0.0)));
    }

    TEST_F(ElementRegionsTest, OffCentreFieldRangeIsTheDefaultAperture) {
        const OffCentreRangeDrift drift(0.1, 0.05);
        const ElementRegions regions = drift.getRegions();
        EXPECT_TRUE(regions.inFieldRegion(Vector({0.12, 0.0, 0.5})));
        EXPECT_FALSE(regions.inFieldRegion(Vector({0.0, 0.0, 0.5})));
        EXPECT_TRUE(regions.inLossRegion(Vector({0.0, 0.0, 0.5})));
        EXPECT_NO_THROW(drift.checkRegions());
    }

    TEST_F(ElementRegionsTest, FieldMapWithoutApertureIsItsMapAndDeletesNothing) {
        OffCentreMapDrift drift;
        const ElementRegions regions = drift.getRegions();
        EXPECT_TRUE(regions.inFieldRegion(Vector({0.7, 0.0, 0.5})));
        EXPECT_FALSE(regions.inFieldRegion(Vector({0.0, 0.0, 0.5})));
        EXPECT_FALSE(regions.inLossRegion(Vector({0.0, 0.0, 0.5})));
        EXPECT_FALSE(regions.inLossRegion(Vector({0.7, 0.3, 0.5})));
        EXPECT_NO_THROW(drift.checkRegions());

        const auto [lower, upper] = drift.getLocalSupportBox(0.5).getCorners();
        EXPECT_NEAR(lower[0], 0.2, 1e-12);
        EXPECT_NEAR(upper[0], 0.8, 1e-12);
        EXPECT_NEAR(upper[1], 0.1, 1e-12);

        // With an APERTURE the element is material out to its reach again.
        drift.setAperture(ApertureType::RECTANGULAR, {0.05, 0.05});
        EXPECT_TRUE(drift.getRegions().inLossRegion(Vector({0.3, 0.0, 0.5})));
    }

    TEST_F(ElementRegionsTest, ApertureWiderThanTheReachIsAnError) {
        MultipoleRep quad = makeQuadrupole(0.8);
        EXPECT_THROW(quad.checkRegions(), GeneralOpalException);
        quad.setReach(1.0);
        EXPECT_NO_THROW(quad.checkRegions());
    }

    TEST_F(ElementRegionsTest, FieldRangeWiderThanTheReachIsAnError) {
        // The range reaches |x| = 0.4 + 0.2 = 0.6 m, beyond the 0.5 m reach.
        OffCentreRangeDrift drift(0.4, 0.2);
        EXPECT_THROW(drift.checkRegions(), GeneralOpalException);
        // An APERTURE inside the reach replaces the range as the aperture.
        drift.setAperture(ApertureType::ELLIPTICAL, {0.3, 0.3});
        EXPECT_NO_THROW(drift.checkRegions());
    }
}  // namespace
