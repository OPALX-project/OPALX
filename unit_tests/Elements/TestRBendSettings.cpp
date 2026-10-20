// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// SPDX-License-Identifier: GPL-3.0-or-later
#include <filesystem>
#include <fstream>
#include "AbsBeamline/RBendProfile.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/RBendRep.h"
#include "BeamlineCore/RFCavityRep.h"
#include "Elements/OpalRBend.h"
#include "Elements/PlacementResolver.h"
#include "Ippl.h"
#include "Physics/Physics.h"
#include "Steppers/BorisPusher.h"
#include "Utilities/OpalException.h"
#include "ValueDefinitions/RealVariable.h"
#include "gtest/gtest.h"

namespace {
    using Vec                     = Vector_t<double, 3>;
    constexpr double electronMass = 510998.95;
    constexpr double energy       = 42.411e6;
    constexpr double angle        = 20.0 * Physics::pi / 180.0;

    struct ProfileFieldFunctor {
        Kokkos::View<Vec*> values;
        RBendProfile::Inputs inputs;

        KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
            values(i) = RBendProfile::field(Vec(0, 0.001, i == 0 ? 0.0 : 0.3), inputs);
        }
    };

    class RBendSettingsTest : public testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
            gmsg                             = new Inform(nullptr, -1);
            [[maybe_unused]] static auto* p0 = new RealVariable();
        }
        static void TearDownTestSuite() {
            delete gmsg;
            gmsg = nullptr;
            ippl::finalize();
        }
        std::string filename;
        void SetUp() override {
            filename = "rbend-settings-"
                       + std::string(testing::UnitTest::GetInstance()->current_test_info()->name())
                       + ".map";
            std::ofstream(filename) << "1DProfile1 1 1 5\n-25 0 15 0\n-15 0 25 0\n0\n2\n0\n2\n";
        }
        void TearDown() override { std::filesystem::remove(filename); }
        std::unique_ptr<OpalRBend> element(double theta = angle, double e1 = 0.0) {
            auto el = std::make_unique<OpalRBend>();
            for (const auto& [name, value] : std::vector<std::pair<std::string, double>>{
                         {"L", 0.3},
                         {"ANGLE", theta},
                         {"E1", e1},
                         {"DESIGNENERGY", 42.411},
                         {"ELEMEDGE", 1.0}})
                Attributes::setReal(*el->findAttribute(name), value);
            Attributes::setString(*el->findAttribute("FMAPFN"), filename);
            return el;
        }
        Vec field(RBend& bend, double x, double y, double z) {
            Vec e(0.0), b(0.0);
            bend.apply(Vec(x, y, z), Vec(0.0), 0.0, e, b);
            return b;
        }
        void awaProfile() {
            std::ofstream(filename) << "1DProfile1 5 5 5\n-25 0 15.2314 0\n-15.2314 0 25 0\n"
                                    << ".478959\n1.91129\n-1.18595\n1.63055\n-1.08266\n.318111\n"
                                    << ".478959\n1.91129\n-1.18595\n1.63055\n-1.08266\n.31811\n";
        }
    };
}  // namespace

TEST_F(RBendSettingsTest, ParserConvertsCentimetresAndUsesMapGap) {
    auto el = element();
    el->update();
    auto& b = dynamic_cast<RBend&>(*el->getElement());
    double begin, end;
    b.getFieldExtent(begin, end);
    EXPECT_DOUBLE_EQ(begin, -0.25);
    EXPECT_DOUBLE_EQ(end, 0.55);
    EXPECT_DOUBLE_EQ(b.getFullGap(), 0.05);
    EXPECT_DOUBLE_EQ(b.getDesignEnergy(), energy);
    EXPECT_DOUBLE_EQ(b.getAperture().second[1], 0.025);
    EXPECT_NEAR(b.getGeometry().getChordLength(), 0.3 / std::cos(angle / 2), 1e-14);
    EXPECT_NEAR(b.getGeometry().getArcLength(), 0.3 * angle / std::sin(angle), 1e-14);
    auto clone = std::unique_ptr<ElementBase>(b.clone());
    EXPECT_EQ(dynamic_cast<RBend&>(*clone).hasFieldMap(), true);
    EXPECT_DOUBLE_EQ(clone->getGeometry().getEntranceAngle(), 0.0);
    clone->setDesignEnergy(1e6);
    EXPECT_DOUBLE_EQ(clone->getDesignEnergy(), energy);
}

TEST_F(RBendSettingsTest, ReferencePathMatchesLegacyFringePlacementAndClones) {
    awaProfile();
    for (double sign : {-1.0, 1.0})
        for (double e1 : {0.0, sign * angle}) {
            auto el = element(sign * angle, e1);
            el->update();
            auto& bend = dynamic_cast<RBend&>(*el->getElement());
            bend.initialiseField(electronMass, -1.0);
            EXPECT_THROW(bend.getPlacementPathLength(), OpalException);
            const auto idealLength = bend.getGeometry().getArcLength();
            const auto idealChord  = bend.getGeometry().getChordLength();
            bend.initialiseReferencePath(electronMass, -1.0, 1e-13);
            // Independent OPAL 2024.1 AWA benchmark, reconstructed from the
            // printed D2 and TDC positions, at the same initialization step.
            const double expected = e1 == 0.0 ? .30594858 : .30590490;
            EXPECT_NEAR(bend.getPlacementPathLength(), expected, 5e-8);
            EXPECT_DOUBLE_EQ(bend.getGeometry().getArcLength(), idealLength);
            EXPECT_DOUBLE_EQ(bend.getGeometry().getChordLength(), idealChord);
            auto clone = std::unique_ptr<ElementBase>(bend.clone());
            EXPECT_DOUBLE_EQ(
                    dynamic_cast<RBend&>(*clone).getPlacementPathLength(),
                    bend.getPlacementPathLength());
            const double first = bend.getPlacementPathLength();
            bend.initialiseReferencePath(electronMass, -1.0, 1e-13);
            EXPECT_DOUBLE_EQ(bend.getPlacementPathLength(), first);
            bend.initialiseField(electronMass, -1.0);
            EXPECT_THROW(bend.getPlacementPathLength(), OpalException);
        }
}

TEST_F(RBendSettingsTest, DoglegPlacesCavityAtLegacyPositionWithoutChangingPhase) {
    awaProfile();
    auto first = element(-angle, 0.0), second = element(angle, angle);
    first->update();
    second->update();
    std::shared_ptr<ElementBase> b1(first->getElement()->clone()),
            b2(second->getElement()->clone());
    b1->setElementPosition(.623686);
    b2->setElementPosition(2.431686);
    for (auto bend : {b1, b2}) {
        auto& b = dynamic_cast<RBend&>(*bend);
        b.initialiseField(electronMass, -1.0);
        b.initialiseReferencePath(electronMass, -1.0, 1e-13);
    }
    auto cavity = std::make_shared<RFCavityRep>("TDC");
    cavity->setElementPosition(3.335485);
    cavity->setPhasem(1.5514735997026097 + Physics::pi / 2.0);
    const double phase   = cavity->getPhasem();
    ElementList elements = {b1, b2, cavity};
    PlacementResolver::resolve(elements, CoordinateSystemTrafo());
    const auto origin = cavity->getCSTrafoGlobal2Local().getOrigin();
    // Actual OPAL benchmark entrance, not ELEMEDGE: both bend-placement
    // passes must use the effective path to reproduce these lab coordinates.
    EXPECT_NEAR(origin(0), .6195280296, 5e-8);
    EXPECT_NEAR(origin(2), 3.233046737, 1e-7);
    const auto direction = cavity->getCSTrafoGlobal2Local().rotateFrom(Vec(0, 0, 1));
    EXPECT_NEAR(direction(0), 0.0, 1e-14);
    EXPECT_NEAR(direction(2), 1.0, 1e-14);
    EXPECT_DOUBLE_EQ(cavity->getPhasem(), phase);
}

TEST_F(RBendSettingsTest, NativePlacementLengthAndInvalidLegacyTimeStep) {
    RBendRep native;
    native.getGeometry() = Geometry::makeRBend(.3, angle);
    EXPECT_DOUBLE_EQ(native.getPlacementPathLength(), native.getGeometry().getArcLength());
    auto el = element();
    el->update();
    auto& bend = dynamic_cast<RBend&>(*el->getElement());
    bend.initialiseField(electronMass, -1.0);
    EXPECT_THROW(bend.initialiseReferencePath(electronMass, -1.0, 0.0), OpalException);
    EXPECT_THROW(bend.initialiseReferencePath(electronMass, -1.0, -1e-13), OpalException);
}

TEST_F(RBendSettingsTest, GapAliasKeepsNativeModelAndItsAperture) {
    OpalRBend a, b;
    for (auto* p : {&a, &b}) {
        Attributes::setReal(*p->findAttribute("L"), .3);
        Attributes::setReal(*p->findAttribute("ELEMEDGE"), 0);
        Attributes::setReal(*p->findAttribute("ANGLE"), angle);
    }
    Attributes::setReal(*a.findAttribute("GAP"), .05);
    Attributes::setReal(*b.findAttribute("HGAP"), .025);
    a.update();
    b.update();
    auto& ar = dynamic_cast<RBend&>(*a.getElement());
    auto& br = dynamic_cast<RBend&>(*b.getElement());
    EXPECT_FALSE(ar.hasFieldMap());
    EXPECT_FALSE(ar.getGeometry().hasPoleFaceAngles());
    EXPECT_EQ(ar.getAperture(), br.getAperture());
    for (double z : {-.05, 0., .15, .3, .35})
        for (int d = 0; d < 3; ++d)
            EXPECT_DOUBLE_EQ(field(ar, .01, .001, z)(d), field(br, .01, .001, z)(d));
}

TEST_F(RBendSettingsTest, CalibrationUsesSpeciesChargeAndDesignEnergy) {
    auto el = element();
    el->update();
    auto& b = dynamic_cast<RBend&>(*el->getElement());
    b.initialiseField(electronMass, -1);
    double electronField = field(b, 0, 0, .1)(1);
    b.initialiseField(electronMass, 1);
    EXPECT_DOUBLE_EQ(field(b, 0, 0, .1)(1), -electronField);
    constexpr double protonMass = 938272088.16;
    b.initialiseField(protonMass, 1);
    EXPECT_NEAR(
            field(b, 0, 0, .1)(1) / (-electronField),
            std::sqrt((energy + 2 * protonMass) / (energy + 2 * electronMass)), 1e-12);
}

TEST_F(RBendSettingsTest, GapAndHgapScaleTheSameProfile) {
    auto a = element(), b = element();
    Attributes::setReal(*a->findAttribute("GAP"), 0.10);
    Attributes::setReal(*b->findAttribute("HGAP"), 0.05);
    a->update();
    b->update();
    for (auto* el : {a.get(), b.get()}) {
        auto& bend = dynamic_cast<RBend&>(*el->getElement());
        bend.initialiseField(electronMass, -1);
        double start, end;
        bend.getFieldExtent(start, end);
        EXPECT_DOUBLE_EQ(start, -0.5);
        EXPECT_DOUBLE_EQ(end, 0.8);
    }
    EXPECT_DOUBLE_EQ(
            field(dynamic_cast<RBend&>(*a->getElement()), 0, 0, .1)(1),
            field(dynamic_cast<RBend&>(*b->getElement()), 0, 0, .1)(1));
    Attributes::setReal(*a->findAttribute("HGAP"), 0.05);
    EXPECT_THROW(a->update(), OpalException);
}

TEST_F(RBendSettingsTest, RejectsInconsistentFacesAndUnsupportedProfileOptions) {
    auto el = element();
    Attributes::setReal(*el->findAttribute("E2"), 0.0);
    EXPECT_THROW(el->update(), OpalException);
    Attributes::setReal(*el->findAttribute("E2"), angle);
    EXPECT_NO_THROW(el->update());
    Attributes::setReal(*el->findAttribute("FINT"), 0.5);
    EXPECT_THROW(el->update(), OpalException);
}

TEST_F(RBendSettingsTest, RejectsMalformedOrUnsupportedMaps) {
    for (const std::string contents :
         {"3DDynamic FALSE\n", "1DProfile1 99 1 5\n",
          "1DProfile1 1 1 -5\n-25 0 15 0\n-15 0 25 0\n0\n2\n0\n2\n",
          "1DProfile1 1 1 5\n-25 0 15 0\n-15 0 25 0\n0\n2\n"}) {
        std::ofstream(filename) << contents;
        EXPECT_THROW(element()->update(), OpalException);
    }
}

TEST_F(RBendSettingsTest, ProfileDerivativesAndDeviceEvaluation) {
    RBendProfile::Inputs in;
    in.enabled        = true;
    in.length         = 0.3;
    in.gap            = 0.05;
    in.amplitude      = 2;
    in.entrance.order = in.exit.order = 1;
    for (auto* p : {&in.entrance, &in.exit}) {
        p->coefficients[0] = 0;
        p->coefficients[1] = 2;
    }
    auto v = RBendProfile::profile(0, in.gap, in.entrance);
    EXPECT_DOUBLE_EQ(v.f, 0.5);
    EXPECT_DOUBLE_EQ(v.first, -10.0);
    EXPECT_DOUBLE_EQ(v.second, 0.0);
    // At the entrance Bz = -B0 F' y, at exit Bz = B0 F' y.
    Kokkos::View<Vec*> values("fringe", 2);
    Kokkos::parallel_for("RBEND profile field test", 2, ProfileFieldFunctor{values, in});
    auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), values);
    EXPECT_DOUBLE_EQ(host(0)(1), 1.0);
    EXPECT_NEAR(host(0)(2), 0.02, 1e-15);
    EXPECT_NEAR(host(1)(2), -0.02, 1e-15);
    RBendProfile::Profile p;
    constexpr double z = .013, h = 1e-6;
    auto a = RBendProfile::profile(z - h, .05, p);
    auto b = RBendProfile::profile(z, .05, p);
    auto c = RBendProfile::profile(z + h, .05, p);
    EXPECT_NEAR(b.first, (c.f - a.f) / (2 * h), 1e-7);
    EXPECT_NEAR(b.second, (c.first - a.first) / (2 * h), 1e-5);
}

TEST_F(RBendSettingsTest, CalibratedFieldTracksBothBendSignsAndBothEntranceAngles) {
    for (double sign : {-1.0, 1.0})
        for (double e1 : {0.0, sign * angle}) {
            auto el = element(sign * angle, e1);
            el->update();
            auto& bend = dynamic_cast<RBend&>(*el->getElement());
            bend.initialiseField(electronMass, -1);
            const double betaGamma = std::sqrt(energy * (energy + 2 * electronMass)) / electronMass;
            Vec p                  = betaGamma * Vec(std::sin(e1), 0, std::cos(e1));
            Vec r(-.25 * std::tan(e1), 0, -.25);
            constexpr double dt    = 1e-13;
            const double halfDrift = .5 * Physics::c * dt / std::sqrt(1 + betaGamma * betaGamma);
            int steps              = 0;
            while (r(2) < .55 && steps++ < 100000) {
                r += halfDrift * p;
                const auto b = field(bend, r(0), 0, r(2));
                BorisPusher().kick(r, p, Vec(0), b, dt, electronMass, -1);
                r += halfDrift * p;
            }
            EXPECT_LT(steps, 100000);
            EXPECT_NEAR(e1 - std::atan2(p(0), p(2)), sign * angle, 2e-7);
            EXPECT_NEAR(std::sqrt(dot(p, p)), betaGamma, 1e-9);
        }
}

TEST_F(RBendSettingsTest, DefaultProfileAndRepeatedInitialisation) {
    auto el = element();
    Attributes::setString(*el->findAttribute("FMAPFN"), "1DPROFILE1-DEFAULT");
    el->update();
    auto& b = dynamic_cast<RBend&>(*el->getElement());
    b.initialiseField(electronMass, -1);
    auto first = field(b, 0, .001, .02);
    b.initialiseField(electronMass, -1);
    auto second = field(b, 0, .001, .02);
    for (int d = 0; d < 3; ++d)
        EXPECT_DOUBLE_EQ(first(d), second(d));
}

TEST_F(RBendSettingsTest, PlacementUsesE1AndPreservesDownstreamReferenceTangent) {
    for (double sign : {-1.0, 1.0})
        for (double e1 : {0.0, sign * angle}) {
            auto bend = std::make_shared<RBendRep>("B");
            auto& g   = bend->getGeometry();
            g         = Geometry::makeRBend(.3, sign * angle);
            g.setPoleFaceAnglesEnabled(true);
            g.setEntranceAngle(e1);
            bend->setElementPosition(1.0);
            auto drift = std::make_shared<DriftRep>("D");
            drift->setElementPosition(1.0 + g.getArcLength() + .4);
            ElementList elements = {bend, drift};
            PlacementResolver::resolve(elements, CoordinateSystemTrafo());
            Vec entryDirection = bend->getCSTrafoGlobal2Local().rotateFrom(Vec(0, 0, 1));
            EXPECT_NEAR(entryDirection(0), -std::sin(e1), 1e-13);
            Vec outgoing = drift->getCSTrafoGlobal2Local().rotateFrom(Vec(0, 0, 1));
            EXPECT_NEAR(outgoing(0), -std::sin(sign * angle), 1e-13);
            auto origin = drift->getCSTrafoGlobal2Local().getOrigin();
            EXPECT_NEAR(
                    origin(0),
                    -g.getChordLength() * std::sin(sign * angle / 2) - .4 * std::sin(sign * angle),
                    1e-13);
            EXPECT_NEAR(
                    origin(2), 1 + g.getChordLength() * std::cos(angle / 2) + .4 * std::cos(angle),
                    1e-13);
        }
}

TEST_F(RBendSettingsTest, PoleFaceAnglesDoNotSelectAProfileOrRequireDesignEnergy) {
    for (bool useExit : {false, true}) {
        OpalRBend el;
        Attributes::setReal(*el.findAttribute("L"), .3);
        Attributes::setReal(*el.findAttribute("ANGLE"), angle);
        Attributes::setReal(*el.findAttribute("ELEMEDGE"), 0.0);
        Attributes::setReal(*el.findAttribute(useExit ? "E2" : "E1"), useExit ? angle : 0.0);
        Attributes::setReal(*el.findAttribute("K1"), .2);
        Attributes::setReal(*el.findAttribute("FINT"), .4);
        ASSERT_NO_THROW(el.update());
        auto& bend = dynamic_cast<RBend&>(*el.getElement());
        EXPECT_FALSE(bend.hasFieldMap());
        EXPECT_DOUBLE_EQ(bend.getDesignEnergy(), 0.0);
        EXPECT_TRUE(bend.getGeometry().hasPoleFaceAngles());
        EXPECT_DOUBLE_EQ(bend.getGeometry().getEntranceAngle(), 0.0);
        EXPECT_DOUBLE_EQ(bend.getGeometry().getExitAngle(), angle);
        EXPECT_NEAR(bend.getGeometry().getArcLength(), .3 * angle / std::sin(angle), 1e-14);
        EXPECT_DOUBLE_EQ(bend.getPlacementPathLength(), bend.getGeometry().getArcLength());
    }
}

TEST_F(RBendSettingsTest, ZeroFintKeepsMaxwellFocusingWithoutAdditionalEdgeField) {
    OpalRBend el;
    Attributes::setReal(*el.findAttribute("L"), .3);
    Attributes::setReal(*el.findAttribute("ANGLE"), angle);
    Attributes::setReal(*el.findAttribute("E1"), 0.0);
    Attributes::setReal(*el.findAttribute("GAP"), .05);
    Attributes::setReal(*el.findAttribute("FINT"), 0.0);
    Attributes::setReal(*el.findAttribute("ELEMEDGE"), 0.0);
    el.update();
    auto& bend = dynamic_cast<RBend&>(*el.getElement());
    bend.setB(1.0);
    for (double z : {0.0, .3}) {
        const auto b = field(bend, 0.0, .001, z);
        EXPECT_NEAR(b(0), 0.0, 1e-15);
        EXPECT_GT(std::abs(b(2)), 1e-5);
    }
}

TEST_F(RBendSettingsTest, DesignEnergyOnlyCalibratesTheStandardField) {
    OpalRBend el;
    Attributes::setReal(*el.findAttribute("L"), .3);
    Attributes::setReal(*el.findAttribute("ANGLE"), angle);
    Attributes::setReal(*el.findAttribute("ELEMEDGE"), 0.0);
    Attributes::setReal(*el.findAttribute("DESIGNENERGY"), 42.411);
    el.update();
    auto& bend = dynamic_cast<RBend&>(*el.getElement());
    EXPECT_FALSE(bend.hasFieldMap());
    EXPECT_FALSE(bend.getGeometry().hasPoleFaceAngles());
    EXPECT_DOUBLE_EQ(bend.getGeometry().getChordLength(), .3);
    bend.setReferenceMomentum(100e6);
    bend.setB(1.0);
    bend.initialiseField(electronMass, 1.0);
    const double expected = std::sqrt(energy * (energy + 2 * electronMass)) / 100e6;
    EXPECT_NEAR(field(bend, 0.0, 0.0, .15)(1), expected, 1e-14);
    bend.initialiseField(electronMass, 1.0);
    EXPECT_NEAR(field(bend, 0.0, 0.0, .15)(1), expected, 1e-14);
    auto clone = std::unique_ptr<ElementBase>(bend.clone());
    EXPECT_NEAR(field(dynamic_cast<RBend&>(*clone), 0.0, 0.0, .15)(1), expected, 1e-14);
}

TEST_F(RBendSettingsTest, ProfileWithoutDesignEnergyUsesReferenceMomentum) {
    OpalRBend el;
    Attributes::setReal(*el.findAttribute("L"), .3);
    Attributes::setReal(*el.findAttribute("ANGLE"), angle);
    Attributes::setReal(*el.findAttribute("E1"), 0.0);
    Attributes::setReal(*el.findAttribute("ELEMEDGE"), 0.0);
    Attributes::setString(*el.findAttribute("FMAPFN"), filename);
    ASSERT_NO_THROW(el.update());
    auto& bend = dynamic_cast<RBend&>(*el.getElement());
    EXPECT_TRUE(bend.hasFieldMap());
    EXPECT_DOUBLE_EQ(bend.getDesignEnergy(), 0.0);
    bend.setReferenceMomentum(std::sqrt(energy * (energy + 2 * electronMass)));
    bend.initialiseField(electronMass, -1.0);
    bend.initialiseReferencePath(electronMass, -1.0, 1e-13);
    auto explicitEnergy = element();
    explicitEnergy->update();
    auto& reference = dynamic_cast<RBend&>(*explicitEnergy->getElement());
    reference.initialiseField(electronMass, -1.0);
    reference.initialiseReferencePath(electronMass, -1.0, 1e-13);
    EXPECT_DOUBLE_EQ(field(bend, 0.0, .001, .1)(1), field(reference, 0.0, .001, .1)(1));
    EXPECT_DOUBLE_EQ(bend.getPlacementPathLength(), reference.getPlacementPathLength());
}

TEST_F(RBendSettingsTest, ProfileDoesNotChangeDefaultPoleFaceGeometry) {
    OpalRBend el;
    Attributes::setReal(*el.findAttribute("L"), .3);
    Attributes::setReal(*el.findAttribute("ANGLE"), angle);
    Attributes::setReal(*el.findAttribute("ELEMEDGE"), 0.0);
    Attributes::setString(*el.findAttribute("FMAPFN"), filename);
    el.update();
    auto& bend = dynamic_cast<RBend&>(*el.getElement());
    EXPECT_TRUE(bend.hasFieldMap());
    EXPECT_FALSE(bend.getGeometry().hasPoleFaceAngles());
    EXPECT_DOUBLE_EQ(bend.getGeometry().getChordLength(), .3);
    EXPECT_DOUBLE_EQ(bend.getGeometry().getEntranceAngle(), .5 * angle);
    EXPECT_DOUBLE_EQ(bend.getGeometry().getExitAngle(), .5 * angle);
    const auto chord = bend.getGeometry().getChordLength();
    bend.setFieldMap("");
    EXPECT_FALSE(bend.hasFieldMap());
    EXPECT_DOUBLE_EQ(bend.getGeometry().getChordLength(), chord);
}

TEST_F(RBendSettingsTest, ProfilePreservesPositionedNegativeBendConvention) {
    OpalRBend standard, mapped;
    for (auto* el : {&standard, &mapped}) {
        Attributes::setReal(*el->findAttribute("L"), .3);
        Attributes::setReal(*el->findAttribute("ANGLE"), -angle);
        Attributes::setReal(*el->findAttribute("ELEMEDGE"), 0.0);
    }
    Attributes::setString(*mapped.findAttribute("FMAPFN"), filename);
    standard.update();
    mapped.update();
    standard.getElement()->fixPosition();
    mapped.getElement()->fixPosition();
    standard.update();
    mapped.update();
    const auto& a = standard.getElement()->getGeometry();
    const auto& b = mapped.getElement()->getGeometry();
    EXPECT_DOUBLE_EQ(a.getBendAngle(), b.getBendAngle());
    EXPECT_DOUBLE_EQ(a.getEntranceAngle(), b.getEntranceAngle());
    EXPECT_DOUBLE_EQ(a.getExitAngle(), b.getExitAngle());
    EXPECT_DOUBLE_EQ(a.getChordLength(), b.getChordLength());
    const auto da = standard.getElement()->getCSTrafoGlobal2Local().rotateFrom(Vec(1, 0, 0));
    const auto db = mapped.getElement()->getCSTrafoGlobal2Local().rotateFrom(Vec(1, 0, 0));
    for (int d = 0; d < 3; ++d)
        EXPECT_DOUBLE_EQ(da(d), db(d));
}
