/**
 * \file TestSpinTBMTPusher.cpp
 * \brief Unit tests for the Thomas-BMT spin pusher (Rodrigues-rotation kernel).
 *
 * The SpinTBMTPusher evolves the per-particle polarization vector
 * @f$\vec P@f$ over a timestep @f$\Delta t@f$ by the Thomas-BMT equation
 * @f$d\vec P / dt = \vec\Omega \times \vec P@f$, with the lab-frame angular
 * velocity
 * @f[
 *   \vec\Omega = -\frac{q}{m}\Big[\,
 *     \big(G + \frac{1}{\gamma}\big)\,\vec B
 *     - \frac{G\,\gamma}{\gamma+1}\,(\vec\beta \cdot \vec B)\,\vec\beta
 *     - \big(G + \frac{1}{\gamma+1}\big)\,\frac{\vec\beta \times \vec E}{c}
 *   \,\Big],
 * @f]
 * applied as an analytic rotation by angle @f$\phi = |\vec\Omega|\,\Delta t@f$
 * about @f$\hat n = \vec\Omega / |\vec\Omega|@f$ via Rodrigues' formula.
 * The rotation form preserves @f$|\vec P|@f$ exactly per step.
 *
 * Here @f$G = (g - 2) / 2@f$ is the magnetic moment anomaly: the dimensionless
 * deviation of the particle's gyromagnetic ratio @f$g@f$ from the Dirac value
 * @f$g = 2@f$. For the muon, @f$G = a_\mu \approx 1.166 \times 10^{-3}@f$;
 * for the proton, @f$G \approx 1.793@f$. @f$G = 0@f$ recovers pure Larmor
 * precession at @f$\omega = qB/(\gamma m)@f$, so the contribution
 * @f$G\,\vec B@f$ in @f$\vec\Omega@f$ is precisely the part the BMT equation
 * adds beyond the Dirac result.
 *
 * References:
 *  - V. Bargmann, L. Michel, V. L. Telegdi, *"Precession of the Polarization of
 *    Particles Moving in a Homogeneous Electromagnetic Field"*,
 *    Phys. Rev. Lett. 2, 435 (1959).
 *  - J. D. Jackson, *Classical Electrodynamics*, 3rd ed., §11.11 (T-BMT eq.).
 *  - PDG listings for the muon (used here for @f$a_\mu = (g-2)/2@f$):
 *    https://pdg.lbl.gov/2024/listings/rpp2024-list-muon.pdf
 *
 * The cyclotron-frequency reference used in the precession check is
 * @f$\omega_c = q\,c^2\,B / (\gamma\,m_{\text{eV}})@f$, which in OPALX
 * conventions equals the SI @f$qB/(\gamma m)@f$ because mass is stored as
 * rest energy @f$m c^2@f$ in eV.
 *
 * Tests in this file:
 *  - **SpinTBMTPusherTest::PrecessesAtCorrectFrequencyInStaticBz** — at rest in
 *    a uniform @f$B_z@f$ with @f$\vec E = 0@f$, T-BMT reduces to Larmor
 *    precession at @f$\omega_S = (1 + a_\mu)\,\omega_c@f$. Pol starts along
 *    @f$\hat x@f$; after one step the rotation angle in the @f$xy@f$ plane is
 *    compared to @f$\omega_S\,\Delta t@f$, and @f$|\vec P|@f$ is verified to
 *    remain 1.
 *  - **SpinTBMTPusherTest::NoEvolutionWhenPolParallelToBAndNoE** — with
 *    @f$\vec P \parallel \vec B@f$ and @f$\vec E = 0@f$, the cross product
 *    @f$\vec\Omega \times \vec P@f$ vanishes; Pol must be unchanged.
 *  - **SpinTBMTPusherTest::NoEvolutionWhenFieldsAreZero** — with
 *    @f$\vec E = \vec B = 0@f$, @f$\vec\Omega = 0@f$ regardless of momentum;
 *    Pol must be unchanged.
 *  - **SpinTBMTPusherTest::PreservesMagnitudeOverManySteps** — drives 1000
 *    rotations under uniform @f$B_z@f$ and checks that @f$|\vec P|@f$ stays
 *    1 to float-storage precision, validating that the Rodrigues form does
 *    not let @f$|\vec P|@f$ drift across many steps.
 */

#include <gtest/gtest.h>

#include "Physics/Physics.h"
#include "Steppers/SpinTBMTPusher.h"
#include "Utility/Inform.h"

#include <cmath>

extern Inform* gmsg;

class SpinTBMTPusherTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        if (gmsg == nullptr) {
            gmsg = new Inform(nullptr, -1);
        }
    }
    static void TearDownTestSuite() {
        if (gmsg != nullptr) {
            delete gmsg;
            gmsg = nullptr;
        }
        ippl::finalize();
    }
};

namespace {

    // Muon: rest mass in eV, charge in proton charges, anomaly G = a_mu.
    constexpr double mass_mu_eV = Physics::m_mu * 1.0e9;
    constexpr double charge_mu  = -1.0;  // mu-
    constexpr double anom_mu    = Physics::a_mu;

    // Pure cyclotron angular frequency for muon in a uniform B (omega_c, no anomaly).
    // omega_c = q * c^2 * B / (gamma * mass), with mass in eV in OPALX conventions.
    double cyclotronFreq(double B, double gamma) {
        return std::fabs(charge_mu) * Physics::c * Physics::c * B / (gamma * mass_mu_eV);
    }

}  // namespace

TEST_F(SpinTBMTPusherTest, PrecessesAtCorrectFrequencyInStaticBz) {
    // Setup: muon at rest (gamma=1, P=0). T-BMT reduces to Larmor precession
    //   omega_S = (q/m) (G + 1) B = (1 + a_mu) * omega_c
    // We compare against (1 + a_mu) * omega_c. Pol starts along x, precesses in xy plane.
    const Vector_t<double, 3> P{0.0, 0.0, 0.0};
    const Vector_t<double, 3> E{0.0, 0.0, 0.0};
    const Vector_t<double, 3> B{0.0, 0.0, 1.0};  // 1 T
    const double dt = 1.0e-12;                   // 1 ps

    ippl::Vector<float, 3> Pol{1.0f, 0.0f, 0.0f};

    SpinTBMTPusher pusher;
    pusher.evolve(Pol, P, E, B, dt, mass_mu_eV, charge_mu, anom_mu);

    const double omegaC      = cyclotronFreq(1.0, 1.0);
    const double expectedPhi = (1.0 + anom_mu) * omegaC * dt;
    // For q<0 the rotation is in the +Bf-aligned axis sense per our prefactor sign convention.
    // We just check magnitude of the rotation angle here.
    const double cosX      = static_cast<double>(Pol[0]);
    const double sinY      = static_cast<double>(Pol[1]);
    const double actualPhi = std::atan2(std::fabs(sinY), cosX);
    EXPECT_NEAR(actualPhi, expectedPhi, 1.0e-9);

    // |Pol| must remain 1.
    const double mag2 = static_cast<double>(Pol[0]) * Pol[0] + static_cast<double>(Pol[1]) * Pol[1]
                        + static_cast<double>(Pol[2]) * Pol[2];
    EXPECT_NEAR(std::sqrt(mag2), 1.0, 1.0e-6);

    // z-component must remain 0 (rotation is about z).
    EXPECT_NEAR(static_cast<double>(Pol[2]), 0.0, 1.0e-7);
}

TEST_F(SpinTBMTPusherTest, NoEvolutionWhenPolParallelToBAndNoE) {
    // B along z, Pol along z, no E -> Omega has zero cross-product effect on Pol along its axis.
    const Vector_t<double, 3> P{0.0, 0.0, 0.0};
    const Vector_t<double, 3> E{0.0, 0.0, 0.0};
    const Vector_t<double, 3> B{0.0, 0.0, 1.0};
    const double dt = 1.0e-9;

    ippl::Vector<float, 3> Pol{0.0f, 0.0f, 1.0f};

    SpinTBMTPusher pusher;
    pusher.evolve(Pol, P, E, B, dt, mass_mu_eV, charge_mu, anom_mu);

    EXPECT_NEAR(static_cast<double>(Pol[0]), 0.0, 1.0e-12);
    EXPECT_NEAR(static_cast<double>(Pol[1]), 0.0, 1.0e-12);
    EXPECT_NEAR(static_cast<double>(Pol[2]), 1.0, 1.0e-7);
}

TEST_F(SpinTBMTPusherTest, NoEvolutionWhenFieldsAreZero) {
    const Vector_t<double, 3> P{0.0, 0.0, 1.0};
    const Vector_t<double, 3> E{0.0, 0.0, 0.0};
    const Vector_t<double, 3> B{0.0, 0.0, 0.0};
    const double dt = 1.0e-9;

    ippl::Vector<float, 3> Pol{1.0f, 0.0f, 0.0f};

    SpinTBMTPusher pusher;
    pusher.evolve(Pol, P, E, B, dt, mass_mu_eV, charge_mu, anom_mu);

    EXPECT_FLOAT_EQ(Pol[0], 1.0f);
    EXPECT_FLOAT_EQ(Pol[1], 0.0f);
    EXPECT_FLOAT_EQ(Pol[2], 0.0f);
}

TEST_F(SpinTBMTPusherTest, PreservesMagnitudeOverManySteps) {
    const Vector_t<double, 3> P{0.0, 0.0, 0.0};
    const Vector_t<double, 3> E{0.0, 0.0, 0.0};
    const Vector_t<double, 3> B{0.0, 0.0, 1.0};
    const double dt = 1.0e-12;

    ippl::Vector<float, 3> Pol{1.0f, 0.0f, 0.0f};
    SpinTBMTPusher pusher;
    for (int i = 0; i < 1000; ++i) {
        pusher.evolve(Pol, P, E, B, dt, mass_mu_eV, charge_mu, anom_mu);
    }
    const double mag2 = static_cast<double>(Pol[0]) * Pol[0] + static_cast<double>(Pol[1]) * Pol[1]
                        + static_cast<double>(Pol[2]) * Pol[2];
    EXPECT_NEAR(std::sqrt(mag2), 1.0, 1.0e-5);
}
