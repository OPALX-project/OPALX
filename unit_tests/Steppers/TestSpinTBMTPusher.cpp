//
// Tests for SpinTBMTPusher: Thomas-BMT spin precession via Rodrigues rotation.
//
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
    const double dt = 1.0e-12;  // 1 ps

    ippl::Vector<float, 3> Pol{1.0f, 0.0f, 0.0f};

    SpinTBMTPusher pusher;
    pusher.evolve(Pol, P, E, B, dt, mass_mu_eV, charge_mu, anom_mu);

    const double omegaC      = cyclotronFreq(1.0, 1.0);
    const double expectedPhi = (1.0 + anom_mu) * omegaC * dt;
    // For q<0 the rotation is in the +Bf-aligned axis sense per our prefactor sign convention.
    // We just check magnitude of the rotation angle here.
    const double cosX = static_cast<double>(Pol[0]);
    const double sinY = static_cast<double>(Pol[1]);
    const double actualPhi = std::atan2(std::fabs(sinY), cosX);
    EXPECT_NEAR(actualPhi, expectedPhi, 1.0e-9);

    // |Pol| must remain 1.
    const double mag2 =
            static_cast<double>(Pol[0]) * Pol[0] + static_cast<double>(Pol[1]) * Pol[1]
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
    const double mag2 = static_cast<double>(Pol[0]) * Pol[0]
                        + static_cast<double>(Pol[1]) * Pol[1]
                        + static_cast<double>(Pol[2]) * Pol[2];
    EXPECT_NEAR(std::sqrt(mag2), 1.0, 1.0e-5);
}
