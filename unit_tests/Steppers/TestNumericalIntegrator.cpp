//
// Tests of the schemes makeNumericalIntegrator() builds: Boris invariants in uniform fields, and
// the order of BORIS and BORIS-PRESTEP in a position-dependent self field. A test double of the
// space-charge algorithm prescribes the self field; the beamline holds no elements.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#include <gtest/gtest.h>

#include "Ippl.h"

#include "../BunchBeam.h"
#include "../KeywordObjects.h"
#include "AbsBeamline/ElementRegions.h"
#include "Distribution/EmittingSamplers.h"
#include "Elements/OpalBeamline.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"
#include "SpaceCharge/SpaceChargeAlgorithm.h"
#include "SpaceCharge/SpaceChargeSolver.h"
#include "Steppers/NumericalIntegrator.h"
#include "Structure/Beam.h"
#include "Utilities/OpalException.h"

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

extern Inform* gmsg;

namespace {
    using namespace opalx::spacecharge;
    using Vector = Vector_t<double, 3>;
    using State  = std::vector<std::array<Vector, 2>>;  ///< R (m) and P (beta gamma) per particle.

    constexpr double eps = std::numeric_limits<double>::epsilon();

    /// Protons with different momenta, so each has its own gamma.
    const State protons{
            {Vector(1e-3, 0.0, 0.0), Vector(0.1, 0.0, 0.0)},
            {Vector(0.0, 0.0, 0.0), Vector(0.0, 0.05, 0.02)},
            {Vector(0.0, -1e-3, 2e-3), Vector(-0.03, 0.04, 0.1)},
            {Vector(2e-3, 1e-3, -1e-3), Vector(0.2, -0.1, -0.05)}};

    /// The prescribed self field: E = e0 + k (x, y, 0) (V/m) and B = b0 (T).
    struct SelfField {
        Vector e0 = Vector(0.0), b0 = Vector(0.0);
        double k = 0.0;  ///< Transverse gradient of E (V/m^2).
    };

    /// Space-charge test double: writes the SelfField at the particle positions into E and B.
    class PrescribedField final : public SpaceChargeAlgorithm {
    public:
        PrescribedField(ParticleContainer_t& pc, const SelfField& field)
            : pc_m(pc), field_m(field) {}

        SpaceChargeSolveResult solve(const SpaceChargeSolveContext&) override {
            auto R            = pc_m.R.getView();
            auto E            = pc_m.E.getView();
            auto B            = pc_m.B.getView();
            const SelfField f = field_m;
            Kokkos::parallel_for(
                    "PrescribedField", pc_m.getLocalNum(), KOKKOS_LAMBDA(const size_t i) {
                        Vector e = f.e0;
                        e(0) += f.k * R(i)(0);
                        e(1) += f.k * R(i)(1);
                        E(i) = e;
                        B(i) = f.b0;
                    });
            Kokkos::fence();
            return {.backendSolves = 1};
        }
        SolveFrame solveFrame() const override { return SolveFrame::Reference; }

    private:
        ParticleContainer_t& pc_m;
        SelfField field_m;
    };

    class NumericalIntegratorTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
            gmsg = new Inform(nullptr, -1);
        }
        static void TearDownTestSuite() {
            delete gmsg;
            gmsg = nullptr;
            ippl::finalize();
        }

        void SetUp() override {
            const auto beam = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
            Attributes::setPredefinedString(*beam->findAttribute("PARTICLE"), "PROTON");
            bunch_m = std::make_unique<PartBunch_t>(
                    std::vector{opalx::test::bunchBeam(*beam, protons.size())}, 1.0,
                    CartesianDomainConfig3D{});
            beamline_m.buildElementMapping(defaultElementMappingGrid, defaultElementReach);
            mass_m   = pc().getReference()->getM();
            charge_m = pc().getReference()->getQ();
        }

        ParticleContainer_t& pc() { return *bunch_m->getParticleContainer(); }

        /// Put the particles at @p state and give each the time step @p dt (s).
        void load(const State& state, double dt) {
            if (pc().getLocalNum() == 0) pc().createParticles(state.size());
            auto R = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc().R.getView());
            auto P = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc().P.getView());
            for (size_t i = 0; i < state.size(); ++i) {
                R(i) = state[i][0];
                P(i) = state[i][1];
            }
            Kokkos::deep_copy(pc().R.getView(), R);
            Kokkos::deep_copy(pc().P.getView(), P);
            bunch_m->resetPcActive();
            bunch_m->setTimeStep(dt);
        }

        /// Advance @p steps steps of scheme @p name with self field @p field; the state after
        /// each step.
        std::vector<State> run(const std::string& name, const SelfField& field, int steps) {
            SpaceChargeSolver solver(std::make_unique<PrescribedField>(pc(), field), 1);
            EmittingSamplers samplers;
            const auto integrator = makeNumericalIntegrator(
                    name, {*bunch_m, beamline_m, solver, samplers, DirichletPlaneConfig{},
                           /*minBinEmitted=*/0, /*aggressiveStateSync=*/false});
            std::vector<State> states;
            for (int step = 0; step < steps; ++step) {
                integrator->advance();
                auto R = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc().R.getView());
                auto P = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc().P.getView());
                State& state = states.emplace_back(pc().getLocalNum());
                for (size_t i = 0; i < state.size(); ++i)
                    state[i] = {R(i), P(i)};
            }
            return states;
        }

        std::unique_ptr<PartBunch_t> bunch_m;
        OpalBeamline beamline_m;
        double mass_m = 0, charge_m = 0;  ///< Rest energy (eV) and charge (e) of a proton.
    };

    double norm(const Vector& v) { return std::sqrt(dot(v, v)); }

}  // namespace

// Each step rotates P about B by 2 atan(|t|), with the half-angle vector
// t = dt q c^2 B / (2 gamma m) of BorisPusher::kick (clockwise for q B_z > 0), and keeps |P|.
TEST_F(NumericalIntegratorTest, UniformMagneticFieldRotatesByTheBorisAngle) {
    const double dt       = 2e-9;  // about 0.2 rad per step in 1 T
    const SelfField field = {.b0 = Vector(0.0, 0.0, 1.0)};
    const int steps       = 50;
    for (const std::string name : {"BORIS", "BORIS-PRESTEP"}) {
        load(protons, dt);
        const std::vector<State> states = run(name, field, steps);
        for (int step = 0; step < steps; ++step) {
            for (size_t i = 0; i < protons.size(); ++i) {
                const Vector& p = (step == 0 ? protons : states[step - 1])[i][1];
                const Vector& q = states[step][i][1];
                const double t  = 0.5 * dt * charge_m * Physics::c * Physics::c
                                 / (std::sqrt(1 + dot(p, p)) * mass_m) * field.b0(2);
                const double angle =
                        std::atan2(p(0) * q(1) - p(1) * q(0), p(0) * q(0) + p(1) * q(1));
                // The rounding of the rotated P, a few eps |P| per component, seen in the plane
                // of P_perp (|P| / |P_perp| <= 2.3 here).
                EXPECT_NEAR(angle, -2 * std::atan(t), 32 * eps) << name << " " << step << " " << i;
            }
        }
        // A rotation keeps |P| exactly; each step adds a few roundings of eps |P|.
        for (size_t i = 0; i < protons.size(); ++i) {
            const double p0 = norm(protons[i][1]);
            EXPECT_NEAR(norm(states.back()[i][1]), p0, 4 * steps * eps * p0) << name << " " << i;
        }
    }
}

// With B = 0, BorisPusher::kick adds two half impulses of dt q c E / (2 m) and rotates by zero.
TEST_F(NumericalIntegratorTest, UniformElectricFieldGivesTheBorisImpulse) {
    const double dt       = 1e-11;
    const SelfField field = {.e0 = Vector(1e6, -2e6, 5e6)};
    const Vector impulse  = dt * charge_m * Physics::c / mass_m * field.e0;
    const int steps       = 20;
    for (const std::string name : {"BORIS", "BORIS-PRESTEP"}) {
        load(protons, dt);
        const std::vector<State> states = run(name, field, steps);
        for (int step = 0; step < steps; ++step) {
            for (size_t i = 0; i < protons.size(); ++i) {
                const Vector& p = (step == 0 ? protons : states[step - 1])[i][1];
                const Vector& q = states[step][i][1];
                // Two roundings of the sums, below eps (|p| + |impulse|), and the rounding of
                // the impulse coefficient, a few eps |impulse|.
                for (unsigned d = 0; d < 3; ++d)
                    EXPECT_NEAR(
                            q(d) - p(d), impulse(d),
                            8 * eps * (std::abs(p(d)) + std::abs(impulse(d))))
                            << name << " " << step << " " << i << " " << d;
            }
        }
    }
}

// A self field that does not depend on the position is the same before and after the first
// half drift, so the two schemes do the same arithmetic.
TEST_F(NumericalIntegratorTest, PositionIndependentFieldMakesTheSchemesEqual) {
    const SelfField field = {.e0 = Vector(2e6, 0.0, -1e6), .b0 = Vector(0.3, -0.2, 0.5)};
    load(protons, 1e-10);
    const State boris = run("BORIS", field, 20).back();
    load(protons, 1e-10);
    const State prestep = run("BORIS-PRESTEP", field, 20).back();
    for (size_t i = 0; i < protons.size(); ++i) {
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_EQ(boris[i][0](d), prestep[i][0](d)) << i << " " << d;
            EXPECT_EQ(boris[i][1](d), prestep[i][1](d)) << i << " " << d;
        }
    }
}

// E = k (x, y, 0) defocuses like the self field of a uniform beam; k gives the growth rate
// Omega = 1/ns, and one proton is tracked over T = 1/Omega (x grows by cosh 1). With
// Omega dt = 1/32 to 1/128, the next order shifts the measured orders by about Omega dt, and the
// errors (3e-9 m to 3e-6 m) are far above the round-off of x (1e-18 m). The reference, BORIS
// with 4096 steps, is about 1e3 times more accurate than the best run measured.
TEST_F(NumericalIntegratorTest, OrderInALinearSelfField) {
    const State start  = {{Vector(1e-3, 0.0, 0.0), Vector(0.0, 0.0, 0.1)}};
    const double omega = 1e9, T = 1 / omega;
    const double gamma    = std::sqrt(1 + dot(start[0][1], start[0][1]));
    const SelfField field = {
            .k = omega * omega * gamma * mass_m / (charge_m * Physics::c * Physics::c)};
    const auto finalX = [&](const std::string& name, int steps) {
        load(start, T / steps);
        return run(name, field, steps).back()[0][0](0);
    };
    const double reference = finalX("BORIS", 4096);
    for (const auto& [name, order] : {std::pair{"BORIS", 2.0}, std::pair{"BORIS-PRESTEP", 1.0}}) {
        std::array<double, 3> error{};
        for (int j = 0; j < 3; ++j)
            error[j] = std::abs(finalX(name, 32 << j) - reference);
        EXPECT_NEAR(std::log2(error[0] / error[1]), order, 0.2) << name;
        EXPECT_NEAR(std::log2(error[1] / error[2]), order, 0.2) << name;
    }
}

TEST_F(NumericalIntegratorTest, UnknownNameThrows) {
    EXPECT_THROW(run("LF2", SelfField{}, 1), OpalException);
}
