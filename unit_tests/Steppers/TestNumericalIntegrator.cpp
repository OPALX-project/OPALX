//
// Tests of the schemes makeNumericalIntegrator() builds: Boris invariants in uniform fields, the
// stage times of RK4 and BORIS in a time-dependent element field, the orders in a quadrupole and in
// a position-dependent self field, and RK4 in a uniform held self field. A test double of the
// space-charge algorithm prescribes the self field; a test that needs an element adds it.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#include <gtest/gtest.h>

#include "Ippl.h"

#include "../BunchBeam.h"
#include "../KeywordObjects.h"
#include "AbsBeamline/ElementRegions.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/VariableRFCavity.h"
#include "Algorithms/PolynomialTimeDependence.h"
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
#include <tuple>
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

    /// A VARIABLE_RF_CAVITY that counts the bunch calls allowed to record (StepTime::record).
    class RecordingCavity final : public VariableRFCavity {
    public:
        RecordingCavity() : VariableRFCavity("RF") {}
        void apply(const std::shared_ptr<ParticleContainer_t>& pc, const StepTime& time) override {
            records += time.record;
            VariableRFCavity::apply(pc, time);
        }
        int records = 0;
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

        /// Add @p element with its local origin at lab z = @p z0 to the beamline.
        void place(const std::shared_ptr<ElementBase>& element, double z0) {
            element->setCSTrafoGlobal2Local(
                    CoordinateSystemTrafo(Vector(0.0, 0.0, z0), Quaternion()));
            element->fixPosition();
            beamline_m.addElement(element, *bunch_m);
            beamline_m.prepareSections();
            beamline_m.buildElementMapping(defaultElementMappingGrid, defaultElementReach);
        }

        /// A cavity of 1 m^3 around all particles with E_z(t) = 10 MV/m sin(omega t + 0.3),
        /// omega = 2 pi 100 MHz, and B = 0.
        std::shared_ptr<RecordingCavity> placeCavity() {
            const auto cavity = std::make_shared<RecordingCavity>();
            cavity->setAmplitudeModel(
                    std::make_shared<PolynomialTimeDependence>(std::vector{10.0}));
            cavity->setFrequencyModel(
                    std::make_shared<PolynomialTimeDependence>(std::vector{100.0}));
            cavity->setPhaseModel(std::make_shared<PolynomialTimeDependence>(std::vector{0.3}));
            cavity->setLength(1.0);
            cavity->setWidth(1.0);
            cavity->setHeight(1.0);
            place(cavity, -0.5);
            return cavity;
        }

        /// A quadrupole of @p gradient (T/m) over lab z in [-0.2, 0.8) m with a 0.1 m bore.
        void placeQuadrupole(double gradient) {
            const auto quadrupole = std::make_shared<Multipole>("Q");
            quadrupole->getGeometry().setElementLength(1.0);
            quadrupole->setNormalComponent(1, gradient);
            quadrupole->setAperture(ApertureType::RECTANGULAR, {0.1, 0.1});
            place(quadrupole, -0.2);
        }

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

        /**
         * @brief For each scheme, one step of the particle @p start from bunch time @p t as bunch
         *        step and as reference step (advanceReference() in the beamline's fields at lab
         *        positions, which are the bunch frame here), compared within their rounding.
         *
         * The reference samples a stage at (t + dt) - dt + c_i dt instead of t + c_i dt, off by
         * less than eps (t + dt), which changes a field by eps-relative amounts; the bunch goes
         * into the element frame (z - z0) and back in each field evaluation, which rounds R_z by
         * up to eps |z0| each time, at most five times. Both get one call per step that may
         * record: the reference field callback and, if @p bunchRecords counts them, the bunch's
         * element.
         */
        void expectReferenceStepIsBunchStep(
                const std::array<Vector, 2>& start, double t, double dt, double z0,
                const int* bunchRecords) {
            for (const std::string name : {"BORIS", "BORIS-PRESTEP", "RK4"}) {
                const int bunchRecordsBefore = bunchRecords ? *bunchRecords : 0;
                load({start}, dt);
                bunch_m->setT(t);
                const auto bunchEnd = run(name, SelfField{}, 1).back()[0];
                if (bunchRecords) EXPECT_EQ(*bunchRecords - bunchRecordsBefore, 1) << name;

                SpaceChargeSolver solver(std::make_unique<PrescribedField>(pc(), SelfField{}), 1);
                EmittingSamplers samplers;
                const auto scheme = makeNumericalIntegrator(
                        name,
                        {*bunch_m, beamline_m, solver, samplers, DirichletPlaneConfig{}, 0, false});
                int calls = 0, records = 0;
                const auto end = scheme->advanceReference(
                        {start[0], start[1]}, dt, t + dt, mass_m, charge_m,
                        [&](const Vector& r, const Vector&, double time, bool record, Vector& E,
                            Vector& B) {
                            ++calls;
                            records += record;
                            return beamline_m.applyOnHost(r, time, E, B);
                        });
                EXPECT_EQ(calls, name == "RK4" ? 4 : 1) << name;
                EXPECT_EQ(records, 1) << name;
                for (unsigned d = 0; d < 3; ++d) {
                    const double gain = std::abs(bunchEnd[1](d) - start[1](d));
                    EXPECT_NEAR(
                            end.position(d), bunchEnd[0](d),
                            5 * eps * (std::abs(z0) + std::abs(bunchEnd[0](d))))
                            << name << " " << d;
                    EXPECT_NEAR(
                            end.momentum(d), bunchEnd[1](d),
                            4 * eps * (std::abs(start[1](d)) + gain))
                            << name << " " << d;
                }
            }
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
    // RK4 holds the self field of the step start through its stages, like BORIS-PRESTEP.
    for (const auto& [name, order] :
         {std::pair{"BORIS", 2.0}, std::pair{"BORIS-PRESTEP", 1.0}, std::pair{"RK4", 1.0}}) {
        std::array<double, 3> error{};
        for (int j = 0; j < 3; ++j)
            error[j] = std::abs(finalX(name, 32 << j) - reference);
        EXPECT_NEAR(std::log2(error[0] / error[1]), order, 0.2) << name;
        EXPECT_NEAR(std::log2(error[1] / error[2]), order, 0.2) << name;
    }
}

// E_z(t) = E0 sin(omega t + phi) of a VARIABLE_RF_CAVITY around all particles, B = 0: the gain of
// P_z in one step is RK4's quadrature of f(t) = q c E_z(t) / M at its stage times, Simpson's
// rule, and the midpoint rule for BORIS; P_x and P_y do not change. With omega dt = 0.63 the two
// rules differ by 1.6 % of the gain. f is sampled with the element's own formula
// (applyOnHost()) at the times StepTime::fieldTime() forms, so the residual is only the rounding
// of the weighted sum and of its addition to P_z: a few eps (|P_z| + |gain|).
TEST_F(NumericalIntegratorTest, StageTimesGiveSimpsonForRk4AndMidpointForBoris) {
    const double t = 1e-9, dt = 1e-9;
    const auto cavity = placeCavity();
    const auto f = [&](double fraction) {
        Vector E(0.0), B(0.0);
        cavity->applyOnHost(Vector(0.0, 0.0, 0.5), t + fraction * dt, E, B);
        return charge_m * Physics::c / mass_m * E(2);
    };
    const double simpson  = dt / 6 * (f(0.0) + 4 * f(0.5) + f(1.0));
    const double midpoint = dt * f(0.5);
    ASSERT_GT(std::abs(simpson - midpoint), 1e-2 * std::abs(simpson));
    for (const auto& [name, gain] : {std::pair{"RK4", simpson}, std::pair{"BORIS", midpoint}}) {
        load(protons, dt);
        bunch_m->setT(t);
        const State end = run(name, SelfField{}, 1).back();
        for (size_t i = 0; i < protons.size(); ++i) {
            const Vector& p = protons[i][1];
            EXPECT_EQ(end[i][1](0), p(0)) << name << " " << i;
            EXPECT_EQ(end[i][1](1), p(1)) << name << " " << i;
            EXPECT_NEAR(end[i][1](2) - p(2), gain, 4 * eps * (std::abs(p(2)) + std::abs(gain)))
                    << name << " " << i;
        }
    }
}

// One proton in a quadrupole of 10 T/m over T = 1/omega, the inverse betatron frequency, at
// beta = 0.1; the reference is RK4 with 1024 steps. With omega dt = 1/8 to 1/32 (RK4) and 1/32 to
// 1/128 (BORIS) the next order shifts the measured orders by less than 0.2, and the smallest RK4
// error, 8e-12 m, is far above the round-off of the position (1e-17 m per step) and of the
// reference (1e-15 m).
TEST_F(NumericalIntegratorTest, OrderInAQuadrupole) {
    const State start     = {{Vector(1e-3, 5e-4, 0.0), Vector(0.0, 0.0, 0.1)}};
    const double gradient = 10.0;  // T/m
    placeQuadrupole(gradient);
    const double gamma = std::sqrt(1 + dot(start[0][1], start[0][1]));
    const double omega = std::sqrt(
            charge_m * gradient * std::pow(Physics::c, 3) * start[0][1](2)
            / (gamma * gamma * mass_m));
    const auto finalR = [&](const std::string& name, int steps) {
        load(start, 1 / omega / steps);
        return run(name, SelfField{}, steps).back()[0][0];
    };
    const Vector reference = finalR("RK4", 1024);
    for (const auto& [name, order, steps] :
         {std::tuple{"RK4", 4.0, 8}, std::tuple{"BORIS", 2.0, 32}}) {
        std::array<double, 3> error{};
        for (int j = 0; j < 3; ++j)
            error[j] = norm(Vector(finalR(name, steps << j) - reference));
        EXPECT_NEAR(std::log2(error[0] / error[1]), order, 0.2) << name;
        EXPECT_NEAR(std::log2(error[1] / error[2]), order, 0.2) << name;
    }
}

// In a uniform B the held self field is the exact field: P rotates about B_z with
// Omega = q c^2 B / (gamma M), clockwise for q B_z > 0, and P_z stays. RK4's step factor for a
// rotation differs from exp(-i Omega dt) by (Omega dt)^5 / 120, so after n steps P is off by about
// n (Omega dt)^5 / 120 |P_perp|, 1e-7 |P_perp| here; the bound allows twice that.
TEST_F(NumericalIntegratorTest, Rk4RotatesInAUniformHeldMagneticField) {
    const double dt       = 5e-10;  // about 0.05 rad per step in 1 T
    const SelfField field = {.b0 = Vector(0.0, 0.0, 1.0)};
    const int steps       = 50;
    load(protons, dt);
    const State end = run("RK4", field, steps).back();
    for (size_t i = 0; i < protons.size(); ++i) {
        const Vector& p    = protons[i][1];
        const double angle = -charge_m * Physics::c * Physics::c * field.b0(2)
                             / (std::sqrt(1 + dot(p, p)) * mass_m) * dt * steps;
        const Vector exact(
                p(0) * std::cos(angle) - p(1) * std::sin(angle),
                p(0) * std::sin(angle) + p(1) * std::cos(angle), p(2));
        const double perp  = std::hypot(p(0), p(1));
        const double theta = std::abs(angle) / steps;
        EXPECT_LE(norm(Vector(end[i][1] - exact)), 2 * steps * std::pow(theta, 5) / 120 * perp)
                << i;
        EXPECT_EQ(end[i][1](2), p(2)) << i;
    }
}

// At t = 1.7 ns and dt = 0.9 ns, (t + dt) - dt differs from t by 2e-25 s.
TEST_F(NumericalIntegratorTest, ReferenceStepIsTheBunchStepInACavity) {
    const auto cavity = placeCavity();
    expectReferenceStepIsBunchStep(
            {Vector(1e-3, -2e-3, 0.0), Vector(0.02, -0.01, 0.1)}, 1.7e-9, 0.9e-9, -0.5,
            &cavity->records);
}

TEST_F(NumericalIntegratorTest, ReferenceStepIsTheBunchStepInAQuadrupole) {
    placeQuadrupole(10.0);
    expectReferenceStepIsBunchStep(
            {Vector(1e-3, 5e-4, 0.0), Vector(0.01, 0.0, 0.1)}, 1.7e-9, 0.9e-9, -0.2, nullptr);
}

TEST_F(NumericalIntegratorTest, UnknownNameThrows) {
    EXPECT_THROW(run("LF2", SelfField{}, 1), OpalException);
}
