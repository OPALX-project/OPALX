/**
 * @file TestBendRep.cpp
 * @brief Regression tests for the OPALX-native SBEND/RBEND representation layer.
 *
 * The tests cover the first port of the analytic bend representation:
 * geometry ownership, dipole-field channels, and the curved placement hooks
 * used by the element-placement bridge.
 */

#include <gtest/gtest.h>

#include "BeamlineCore/RBendRep.h"
#include "BeamlineCore/SBendRep.h"
#include "Channels/Channel.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"

#include <cmath>
#include <memory>
#include <stdexcept>

namespace {

    struct ReferenceState {
        Vector_t<double, 3> r;
        Vector_t<double, 3> p;
    };

    /**
     * @brief Compute the dimensionless momentum \f$\beta\gamma\f$ from kinetic energy.
     *
     * The relativistic momentum variable used by OPALX tracking is
     * \f[
     * \mathbf{p} = \beta\gamma \hat{\mathbf{u}},
     * \f]
     * where \f$\gamma = 1 + E_\mathrm{kin}/m c^2\f$.
     */
    double betaGammaFromKineticEnergy(const double kineticEnergyEV, const double massEV) {
        const double gamma = 1.0 + kineticEnergyEV / massEV;
        return std::sqrt(gamma * gamma - 1.0);
    }

    /**
     * @brief Relativistic magnetic rigidity relation used for the bend field.
     *
     * For motion in a pure dipole field the design radius satisfies
     * \f[
     * \rho = \frac{\beta\gamma m}{|q| c |B|},
     * \f]
     * so the field magnitude required for a target radius is
     * \f[
     * |B| = \frac{\beta\gamma m}{|q| c \rho}.
     * \f]
     */
    double dipoleFieldMagnitude(
            const double betaGamma, const double massEV, const double charge, const double radius) {
        return betaGamma * massEV / (Physics::c * std::abs(charge) * radius);
    }

    ReferenceState evaluateDerivative(
            BendBase& bend, const ReferenceState& state, const double massEV, const double charge) {
        Vector_t<double, 3> E(0.0);
        Vector_t<double, 3> B(0.0);
        if (bend.applyToReferenceParticle(state.r, state.p, 0.0, E, B)) {
            throw std::runtime_error(
                    "Reference particle left bend field support during regression.");
        }

        const double gamma = std::sqrt(1.0 + dot(state.p, state.p));
        ReferenceState derivative{};
        derivative.r = Physics::c * state.p / gamma;
        derivative.p = charge * Physics::c / massEV * E
                       + charge * Physics::c * Physics::c / (gamma * massEV) * cross(state.p, B);
        return derivative;
    }

    ReferenceState integrateReferenceParticle(
            BendBase& bend, const double massEV, const double charge, const double kineticEnergyEV,
            const double integrationTime, const std::size_t steps, const double initialZ = 0.0) {
        ReferenceState state{};
        state.r = Vector_t<double, 3>(0.0, 0.0, initialZ);
        state.p =
                Vector_t<double, 3>(0.0, 0.0, betaGammaFromKineticEnergy(kineticEnergyEV, massEV));

        const double dt = integrationTime / static_cast<double>(steps);
        for (std::size_t i = 0; i < steps; ++i) {
            const ReferenceState k1 = evaluateDerivative(bend, state, massEV, charge);

            ReferenceState tmp = state;
            tmp.r += 0.5 * dt * k1.r;
            tmp.p += 0.5 * dt * k1.p;
            const ReferenceState k2 = evaluateDerivative(bend, tmp, massEV, charge);

            tmp = state;
            tmp.r += 0.5 * dt * k2.r;
            tmp.p += 0.5 * dt * k2.p;
            const ReferenceState k3 = evaluateDerivative(bend, tmp, massEV, charge);

            tmp = state;
            tmp.r += dt * k3.r;
            tmp.p += dt * k3.p;
            const ReferenceState k4 = evaluateDerivative(bend, tmp, massEV, charge);

            state.r += dt / 6.0 * (k1.r + 2.0 * k2.r + 2.0 * k3.r + k4.r);
            state.p += dt / 6.0 * (k1.p + 2.0 * k2.p + 2.0 * k3.p + k4.p);
        }

        return state;
    }

    ReferenceState integrateReferenceParticle(
            BendBase& bend, const double massEV, const double charge,
            const ReferenceState& initialState, const double integrationTime,
            const std::size_t steps) {
        ReferenceState state = initialState;

        const double dt = integrationTime / static_cast<double>(steps);
        for (std::size_t i = 0; i < steps; ++i) {
            const ReferenceState k1 = evaluateDerivative(bend, state, massEV, charge);

            ReferenceState tmp = state;
            tmp.r += 0.5 * dt * k1.r;
            tmp.p += 0.5 * dt * k1.p;
            const ReferenceState k2 = evaluateDerivative(bend, tmp, massEV, charge);

            tmp = state;
            tmp.r += 0.5 * dt * k2.r;
            tmp.p += 0.5 * dt * k2.p;
            const ReferenceState k3 = evaluateDerivative(bend, tmp, massEV, charge);

            tmp = state;
            tmp.r += dt * k3.r;
            tmp.p += dt * k3.p;
            const ReferenceState k4 = evaluateDerivative(bend, tmp, massEV, charge);

            state.r += dt / 6.0 * (k1.r + 2.0 * k2.r + 2.0 * k3.r + k4.r);
            state.p += dt / 6.0 * (k1.p + 2.0 * k2.p + 2.0 * k3.p + k4.p);
        }

        return state;
    }

    TEST(SBendRep, GeometryAndDesignPathFollowPlanarArc) {
        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(2.0, M_PI / 12.0);
        bend.setElementLength(2.0);
        bend.setBendAngle(M_PI / 6.0);
        bend.setExitAngle(M_PI / 18.0);

        EXPECT_EQ(bend.getType(), ElementType::SBEND);
        EXPECT_NEAR(bend.getGeometry().getElementLength(), 2.0, 1e-12);
        EXPECT_NEAR(bend.getBendAngle(), M_PI / 6.0, 1e-12);

        const double radius = 2.0 / (M_PI / 6.0);
        const double chord  = 2.0 * radius * std::sin(M_PI / 12.0);
        EXPECT_NEAR(bend.getChordLength(), chord, 1e-12);

        const auto path = bend.getDesignPath(7);
        ASSERT_GE(path.size(), 7u);

        const auto entry = bend.getEdgeToBegin().getOrigin();
        const auto exit  = bend.getEdgeToEnd().getOrigin();
        EXPECT_NEAR(path.front()(0), entry(0), 1e-12);
        EXPECT_NEAR(path.front()(1), entry(1), 1e-12);
        EXPECT_NEAR(path.front()(2), entry(2), 1e-12);
        EXPECT_NEAR(path.back()(0), exit(0), 1e-12);
        EXPECT_NEAR(path.back()(1), exit(1), 1e-12);
        EXPECT_NEAR(path.back()(2), exit(2), 1e-12);
    }

    TEST(SBendRep, DipoleChannelUsesCurrentZeroBasedMultipoleConvention) {
        SBendRep bend("SBEND");

        Channel* by = bend.getChannel("BY", false);
        ASSERT_NE(by, nullptr);
        EXPECT_TRUE(by->set(1.75));
        EXPECT_DOUBLE_EQ(bend.getB(), 1.75);
        EXPECT_DOUBLE_EQ(bend.getField().getNormalComponent(0), 1.75);
        delete by;
    }

    TEST(RBendRep, ExitAngleTracksRectangularBendConvention) {
        RBendRep bend("RBEND");
        bend.getGeometry().setElementLength(1.2);
        bend.getGeometry().setBendAngle(M_PI / 9.0);
        bend.setElementLength(1.2);
        bend.setBendAngle(M_PI / 9.0);
        bend.setEntranceAngle(M_PI / 18.0);

        EXPECT_EQ(bend.getType(), ElementType::RBEND);
        EXPECT_NEAR(bend.getExitAngle(), M_PI / 18.0, 1e-12);
        EXPECT_NEAR(bend.getChordLength(), 1.2, 1e-12);

        const auto path = bend.getDesignPath(5);
        ASSERT_GE(path.size(), 5u);
        const auto entry = bend.getEdgeToBegin().getOrigin();
        const auto exit  = bend.getEdgeToEnd().getOrigin();
        EXPECT_NEAR(path.front()(0), entry(0), 1e-12);
        EXPECT_NEAR(path.back()(0), exit(0), 1e-12);
    }

    TEST(RBendRep, ClonePreservesGeometryAndField) {
        RBendRep bend("RBEND");
        bend.getGeometry().setElementLength(0.8);
        bend.getGeometry().setBendAngle(M_PI / 8.0);
        bend.setElementLength(0.8);
        bend.setBendAngle(M_PI / 8.0);
        bend.setB(0.9);

        std::unique_ptr<ElementBase> clone(bend.clone());
        auto* copy = dynamic_cast<RBendRep*>(clone.get());
        ASSERT_NE(copy, nullptr);
        EXPECT_NEAR(copy->getElementLength(), 0.8, 1e-12);
        EXPECT_NEAR(copy->getBendAngle(), M_PI / 8.0, 1e-12);
        EXPECT_NEAR(copy->getB(), 0.9, 1e-12);
    }

    TEST(SBendRep, ProtonOnAxisIntegratesToFortyFiveDegreeDeflection) {
        constexpr double angle           = Physics::pi / 4.0;
        constexpr double arcLength       = 1.0;
        constexpr double kineticEnergyEV = 590.0 * Units::MeV2eV;
        constexpr double massEV          = Physics::m_p * Units::GeV2eV;
        constexpr double charge          = 1.0;
        constexpr std::size_t steps      = 20000;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(arcLength, angle / arcLength);
        bend.setElementLength(arcLength);
        bend.setBendAngle(angle);

        const double radius    = arcLength / angle;
        const double betaGamma = betaGammaFromKineticEnergy(kineticEnergyEV, massEV);
        const double by        = -dipoleFieldMagnitude(betaGamma, massEV, charge, radius);
        bend.setB(by);
        bend.setFieldAmplitude(by, 0.0);

        const double beta            = betaGamma / std::sqrt(1.0 + betaGamma * betaGamma);
        const double integrationTime = arcLength / (beta * Physics::c);
        const ReferenceState state   = integrateReferenceParticle(
                bend, massEV, charge, kineticEnergyEV, integrationTime, steps);

        const double deflection = std::atan2(state.p(0), state.p(2));
        EXPECT_NEAR(deflection, angle, 5.0e-5);
        EXPECT_NEAR(state.r(0), radius * (1.0 - std::cos(angle)), 1.0e-6);
        EXPECT_NEAR(state.r(2), radius * std::sin(angle), 1.0e-6);
    }

    TEST(RBendRep, ElectronOnAxisIntegratesToFortyFiveDegreeDeflection) {
        constexpr double angle           = Physics::pi / 4.0;
        constexpr double chordLength     = 1.0;
        constexpr double kineticEnergyEV = 1.0 * Units::GeV2eV;
        constexpr double massEV          = Physics::m_e * Units::GeV2eV;
        constexpr double charge          = -1.0;
        constexpr std::size_t steps      = 20000;

        RBendRep bend("RBEND");
        bend.getGeometry().setElementLength(chordLength);
        bend.getGeometry().setBendAngle(angle);
        bend.setElementLength(chordLength);
        bend.setBendAngle(angle);

        const double radius    = chordLength / (2.0 * std::sin(angle / 2.0));
        const double arcLength = radius * angle;
        const double betaGamma = betaGammaFromKineticEnergy(kineticEnergyEV, massEV);
        const double by        = dipoleFieldMagnitude(betaGamma, massEV, charge, radius);
        bend.setB(by);
        bend.setFieldAmplitude(by, 0.0);

        const double beta            = betaGamma / std::sqrt(1.0 + betaGamma * betaGamma);
        const double integrationTime = arcLength / (beta * Physics::c);
        const ReferenceState state   = integrateReferenceParticle(
                bend, massEV, charge, kineticEnergyEV, integrationTime, steps);

        const double deflection = std::atan2(state.p(0), state.p(2));
        EXPECT_NEAR(deflection, angle, 5.0e-5);
        EXPECT_NEAR(state.r(0), radius * (1.0 - std::cos(angle)), 1.0e-6);
        EXPECT_NEAR(state.r(2), radius * std::sin(angle), 1.0e-6);
    }

    TEST(SBendRep, FringeSupportProducesExpectedLinearFieldRamp) {
        constexpr double bodyLength = 1.0;
        constexpr double angle      = Physics::pi / 4.0;
        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, angle / bodyLength);
        bend.setLength(bodyLength);
        bend.setBendAngle(angle);
        bend.setFringeHalfGap(0.02);
        bend.setFringeIntegral(0.5);
        bend.setB(-1.0);

        double fieldBegin = 0.0;
        double fieldEnd   = 0.0;
        bend.getFieldExtend(fieldBegin, fieldEnd);

        Vector_t<double, 3> E(0.0);
        Vector_t<double, 3> B(0.0);
        EXPECT_NEAR(fieldBegin, -0.01, 1.0e-12);
        EXPECT_NEAR(fieldEnd, 1.01, 1.0e-12);

        EXPECT_FALSE(bend.applyToReferenceParticle(
                Vector_t<double, 3>(0.0, 0.0, -0.02), Vector_t<double, 3>(0.0), 0.0, E, B));
        EXPECT_NEAR(B(1), 0.0, 1.0e-12);

        B = Vector_t<double, 3>(0.0);
        EXPECT_FALSE(bend.applyToReferenceParticle(
                Vector_t<double, 3>(0.0, 0.0, -0.005), Vector_t<double, 3>(0.0), 0.0, E, B));
        EXPECT_NEAR(B(1), -0.5, 1.0e-12);

        B = Vector_t<double, 3>(0.0);
        EXPECT_FALSE(bend.applyToReferenceParticle(
                Vector_t<double, 3>(0.0, 0.0, 0.5), Vector_t<double, 3>(0.0), 0.0, E, B));
        EXPECT_NEAR(B(1), -1.0, 1.0e-12);

        B = Vector_t<double, 3>(0.0);
        EXPECT_FALSE(bend.applyToReferenceParticle(
                Vector_t<double, 3>(0.0, 0.0, 1.005), Vector_t<double, 3>(0.0), 0.0, E, B));
        EXPECT_NEAR(B(1), -0.5, 1.0e-12);
    }

    TEST(SBendRep, EntryFringeProducesEdgeFocusingKick) {
        constexpr double bodyLength      = 1.0;
        constexpr double bendAngle       = 0.2;
        constexpr double entryAngle      = 0.15;
        constexpr double fringeHalfGap   = 0.02;
        constexpr double fringeIntegral  = 0.5;
        constexpr double kineticEnergyEV = 590.0 * Units::MeV2eV;
        constexpr double massEV          = Physics::m_p * Units::GeV2eV;
        constexpr double charge          = 1.0;
        constexpr std::size_t steps      = 20000;
        constexpr double xOffset         = 1.0e-4;
        constexpr double yOffset         = 1.0e-4;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(entryAngle);
        bend.setExitAngle(0.0);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);

        const double effectiveLength = bend.getEffectiveFieldLength();
        const double radius          = effectiveLength / bendAngle;
        const double betaGamma       = betaGammaFromKineticEnergy(kineticEnergyEV, massEV);
        const double by              = -dipoleFieldMagnitude(betaGamma, massEV, charge, radius);
        bend.setB(by);
        bend.setFieldAmplitude(bendAngle / effectiveLength, 0.0);

        double fieldBegin = 0.0;
        double fieldEnd   = 0.0;
        bend.getFieldExtend(fieldBegin, fieldEnd);
        const double beta            = betaGamma / std::sqrt(1.0 + betaGamma * betaGamma);
        const double integrationTime = (-fieldBegin) / (beta * Physics::c);

        const ReferenceState base = integrateReferenceParticle(
                bend, massEV, charge, kineticEnergyEV, integrationTime, steps, fieldBegin);

        ReferenceState xState{};
        xState.r = Vector_t<double, 3>(xOffset, 0.0, fieldBegin);
        xState.p = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        const ReferenceState xFocused =
                integrateReferenceParticle(bend, massEV, charge, xState, integrationTime, steps);

        ReferenceState yState{};
        yState.r = Vector_t<double, 3>(0.0, yOffset, fieldBegin);
        yState.p = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        const ReferenceState yFocused =
                integrateReferenceParticle(bend, massEV, charge, yState, integrationTime, steps);

        const double signedCurvature = bendAngle / effectiveLength;
        const double psi             = signedCurvature * fringeHalfGap * fringeIntegral
                           * (1.0 + std::sin(entryAngle) * std::sin(entryAngle))
                           / std::cos(entryAngle);
        const double expectedHorizontal = signedCurvature * std::tan(entryAngle);
        const double expectedVertical   = -signedCurvature * std::tan(entryAngle - psi);
        const double deltaXP            = (xFocused.p(0) - base.p(0)) / betaGamma;
        const double deltaYP            = (yFocused.p(1) - base.p(1)) / betaGamma;

        EXPECT_NEAR(deltaXP / xOffset, expectedHorizontal, 5.0e-4);
        EXPECT_NEAR(deltaYP / yOffset, expectedVertical, 5.0e-4);
        EXPECT_LT(std::abs(expectedVertical), std::abs(expectedHorizontal));
    }

    TEST(SBendRep, ExitFringeProducesEdgeFocusingKick) {
        constexpr double bodyLength      = 1.0;
        constexpr double bendAngle       = 0.2;
        constexpr double exitAngle       = 0.14;
        constexpr double fringeHalfGap   = 0.02;
        constexpr double fringeIntegral  = 0.5;
        constexpr double kineticEnergyEV = 590.0 * Units::MeV2eV;
        constexpr double massEV          = Physics::m_p * Units::GeV2eV;
        constexpr double charge          = 1.0;
        constexpr std::size_t steps      = 20000;
        constexpr double xOffset         = 1.0e-4;
        constexpr double yOffset         = 1.0e-4;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(0.0);
        bend.setExitAngle(exitAngle);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);

        const double effectiveLength = bend.getEffectiveFieldLength();
        const double radius          = effectiveLength / bendAngle;
        const double betaGamma       = betaGammaFromKineticEnergy(kineticEnergyEV, massEV);
        const double by              = -dipoleFieldMagnitude(betaGamma, massEV, charge, radius);
        bend.setB(by);
        bend.setFieldAmplitude(bendAngle / effectiveLength, 0.0);

        double fieldBegin = 0.0;
        double fieldEnd   = 0.0;
        bend.getFieldExtend(fieldBegin, fieldEnd);
        const double beta            = betaGamma / std::sqrt(1.0 + betaGamma * betaGamma);
        const double integrationTime = (fieldEnd - bodyLength) / (beta * Physics::c);

        const ReferenceState base = integrateReferenceParticle(
                bend, massEV, charge, kineticEnergyEV, integrationTime, steps, bodyLength);

        ReferenceState xState{};
        xState.r = Vector_t<double, 3>(xOffset, 0.0, bodyLength);
        xState.p = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        const ReferenceState xFocused =
                integrateReferenceParticle(bend, massEV, charge, xState, integrationTime, steps);

        ReferenceState yState{};
        yState.r = Vector_t<double, 3>(0.0, yOffset, bodyLength);
        yState.p = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        const ReferenceState yFocused =
                integrateReferenceParticle(bend, massEV, charge, yState, integrationTime, steps);

        const double signedCurvature = bendAngle / effectiveLength;
        const double psi             = signedCurvature * fringeHalfGap * fringeIntegral
                           * (1.0 + std::sin(exitAngle) * std::sin(exitAngle))
                           / std::cos(exitAngle);
        const double expectedHorizontal = signedCurvature * std::tan(exitAngle);
        const double expectedVertical   = -signedCurvature * std::tan(exitAngle - psi);
        const double deltaXP            = (xFocused.p(0) - base.p(0)) / betaGamma;
        const double deltaYP            = (yFocused.p(1) - base.p(1)) / betaGamma;

        EXPECT_NEAR(deltaXP / xOffset, expectedHorizontal, 5.0e-4);
        EXPECT_NEAR(deltaYP / yOffset, expectedVertical, 5.0e-4);
        EXPECT_LT(std::abs(expectedVertical), std::abs(expectedHorizontal));
    }

    TEST(RBendRep, EntryFringeProducesEdgeFocusingKick) {
        constexpr double bodyLength      = 1.0;
        constexpr double bendAngle       = 0.2;
        constexpr double entryAngle      = 0.12;
        constexpr double fringeHalfGap   = 0.02;
        constexpr double fringeIntegral  = 0.5;
        constexpr double kineticEnergyEV = 1.0 * Units::GeV2eV;
        constexpr double massEV          = Physics::m_e * Units::GeV2eV;
        constexpr double charge          = -1.0;
        constexpr std::size_t steps      = 20000;
        constexpr double xOffset         = 1.0e-4;
        constexpr double yOffset         = 1.0e-4;

        RBendRep bend("RBEND");
        bend.getGeometry().setElementLength(bodyLength);
        bend.getGeometry().setBendAngle(bendAngle);
        bend.setLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(entryAngle);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);

        const double effectiveLength = bend.getEffectiveFieldLength();
        const double radius          = effectiveLength / bendAngle;
        const double betaGamma       = betaGammaFromKineticEnergy(kineticEnergyEV, massEV);
        const double by              = dipoleFieldMagnitude(betaGamma, massEV, charge, radius);
        bend.setB(by);
        bend.setFieldAmplitude(bendAngle / effectiveLength, 0.0);

        double fieldBegin = 0.0;
        double fieldEnd   = 0.0;
        bend.getFieldExtend(fieldBegin, fieldEnd);
        const double beta            = betaGamma / std::sqrt(1.0 + betaGamma * betaGamma);
        const double integrationTime = (-fieldBegin) / (beta * Physics::c);

        const ReferenceState base = integrateReferenceParticle(
                bend, massEV, charge, kineticEnergyEV, integrationTime, steps, fieldBegin);

        ReferenceState xState{};
        xState.r = Vector_t<double, 3>(xOffset, 0.0, fieldBegin);
        xState.p = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        const ReferenceState xFocused =
                integrateReferenceParticle(bend, massEV, charge, xState, integrationTime, steps);

        ReferenceState yState{};
        yState.r = Vector_t<double, 3>(0.0, yOffset, fieldBegin);
        yState.p = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        const ReferenceState yFocused =
                integrateReferenceParticle(bend, massEV, charge, yState, integrationTime, steps);

        const double signedCurvature = bendAngle / effectiveLength;
        const double psi             = signedCurvature * fringeHalfGap * fringeIntegral
                           * (1.0 + std::sin(entryAngle) * std::sin(entryAngle))
                           / std::cos(entryAngle);
        const double expectedHorizontal = signedCurvature * std::tan(entryAngle);
        const double expectedVertical   = -signedCurvature * std::tan(entryAngle - psi);
        const double deltaXP            = (xFocused.p(0) - base.p(0)) / betaGamma;
        const double deltaYP            = (yFocused.p(1) - base.p(1)) / betaGamma;

        EXPECT_NEAR(deltaXP / xOffset, expectedHorizontal, 5.0e-4);
        EXPECT_NEAR(deltaYP / yOffset, expectedVertical, 5.0e-4);
        EXPECT_LT(std::abs(expectedVertical), std::abs(expectedHorizontal));
    }

    TEST(RBendRep, ExitFringeProducesEdgeFocusingKick) {
        constexpr double bodyLength      = 1.0;
        constexpr double bendAngle       = 0.2;
        constexpr double exitAngle       = 0.11;
        constexpr double fringeHalfGap   = 0.02;
        constexpr double fringeIntegral  = 0.5;
        constexpr double kineticEnergyEV = 1.0 * Units::GeV2eV;
        constexpr double massEV          = Physics::m_e * Units::GeV2eV;
        constexpr double charge          = -1.0;
        constexpr std::size_t steps      = 20000;
        constexpr double xOffset         = 1.0e-4;
        constexpr double yOffset         = 1.0e-4;

        RBendRep bend("RBEND");
        bend.getGeometry().setElementLength(bodyLength);
        bend.getGeometry().setBendAngle(bendAngle);
        bend.setLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(bendAngle - exitAngle);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);

        const double effectiveLength = bend.getEffectiveFieldLength();
        const double radius          = effectiveLength / bendAngle;
        const double betaGamma       = betaGammaFromKineticEnergy(kineticEnergyEV, massEV);
        const double by              = dipoleFieldMagnitude(betaGamma, massEV, charge, radius);
        bend.setB(by);
        bend.setFieldAmplitude(bendAngle / effectiveLength, 0.0);

        double fieldBegin = 0.0;
        double fieldEnd   = 0.0;
        bend.getFieldExtend(fieldBegin, fieldEnd);
        const double beta            = betaGamma / std::sqrt(1.0 + betaGamma * betaGamma);
        const double integrationTime = (fieldEnd - bodyLength) / (beta * Physics::c);

        const ReferenceState base = integrateReferenceParticle(
                bend, massEV, charge, kineticEnergyEV, integrationTime, steps, bodyLength);

        ReferenceState xState{};
        xState.r = Vector_t<double, 3>(xOffset, 0.0, bodyLength);
        xState.p = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        const ReferenceState xFocused =
                integrateReferenceParticle(bend, massEV, charge, xState, integrationTime, steps);

        ReferenceState yState{};
        yState.r = Vector_t<double, 3>(0.0, yOffset, bodyLength);
        yState.p = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        const ReferenceState yFocused =
                integrateReferenceParticle(bend, massEV, charge, yState, integrationTime, steps);

        const double signedCurvature = bendAngle / effectiveLength;
        const double psi             = signedCurvature * fringeHalfGap * fringeIntegral
                           * (1.0 + std::sin(exitAngle) * std::sin(exitAngle))
                           / std::cos(exitAngle);
        const double expectedHorizontal = signedCurvature * std::tan(exitAngle);
        const double expectedVertical   = -signedCurvature * std::tan(exitAngle - psi);
        const double deltaXP            = (xFocused.p(0) - base.p(0)) / betaGamma;
        const double deltaYP            = (yFocused.p(1) - base.p(1)) / betaGamma;

        EXPECT_NEAR(deltaXP / xOffset, expectedHorizontal, 5.0e-4);
        EXPECT_NEAR(deltaYP / yOffset, expectedVertical, 5.0e-4);
        EXPECT_LT(std::abs(expectedVertical), std::abs(expectedHorizontal));
    }

}  // namespace
