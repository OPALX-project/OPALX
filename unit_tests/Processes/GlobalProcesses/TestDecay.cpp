#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>

#include "Ippl.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Physics/MuonDecay.h"
#include "Physics/ParticleProperties.h"
#include "Physics/Physics.h"
#include "Processes/GlobalProcesses/MuonDecay.h"
#include "Processes/GlobalProcesses/PionDecay.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

namespace {

    using PC_t = ParticleContainer<double, 3>;

    class DecayTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }

        static void TearDownTestSuite() { ippl::finalize(); }

        void SetUp() override {
            oldSeed_m                = Options::seed;
            oldUseQMAttributes_m     = Options::useQMAttributes;
            Options::useQMAttributes = false;
        }

        void TearDown() override {
            Options::seed            = oldSeed_m;
            Options::useQMAttributes = oldUseQMAttributes_m;
        }

        std::shared_ptr<PC_t> makeContainer(ParticleType species = ParticleType::UNNAMED) {
            ippl::Vector<int, 3> nr        = 8;
            ippl::Vector<double, 3> rmin   = -4.0;
            ippl::Vector<double, 3> rmax   = 4.0;
            ippl::Vector<double, 3> origin = rmin;
            ippl::Vector<double, 3> hr     = (rmax - rmin) / ippl::Vector<double, 3>(nr);
            std::array<bool, 3> decomp     = {true, true, true};

            ippl::NDIndex<3> domain;
            for (unsigned i = 0; i < 3; ++i) {
                domain[i] = ippl::Index(nr[i]);
            }

            Mesh_t<3> mesh(domain, hr, origin);
            FieldLayout_t<3> fl(MPI_COMM_WORLD, domain, decomp, true);
            auto pc = std::make_shared<PC_t>(mesh, fl);
            pc->Sp  = static_cast<short>(species);
            return pc;
        }

        void createParticles(std::shared_ptr<PC_t>& pc, size_t nPart, double pz) {
            if (nPart == 0) {
                return;
            }

            pc->create(nPart);
            auto R_host = pc->R.getHostMirror();
            auto P_host = pc->P.getHostMirror();

            for (size_t i = 0; i < nPart; ++i) {
                R_host(i)[0] = 0.0;
                R_host(i)[1] = 0.0;
                R_host(i)[2] = 0.0;
                P_host(i)[0] = 0.0;
                P_host(i)[1] = 0.0;
                P_host(i)[2] = pz;
            }

            Kokkos::deep_copy(pc->R.getView(), R_host);
            Kokkos::deep_copy(pc->P.getView(), P_host);
            Kokkos::fence();
        }

        void createParticlesWithPositions(std::shared_ptr<PC_t>& pc, size_t nPart, double pz) {
            if (nPart == 0) {
                return;
            }

            pc->create(nPart);
            auto R_host = pc->R.getHostMirror();
            auto P_host = pc->P.getHostMirror();

            for (size_t i = 0; i < nPart; ++i) {
                R_host(i)[0] = 0.1 * static_cast<double>(i);
                R_host(i)[1] = 0.2 * static_cast<double>(i);
                R_host(i)[2] = 0.3 * static_cast<double>(i);
                P_host(i)[0] = 0.0;
                P_host(i)[1] = 0.0;
                P_host(i)[2] = pz;
            }

            Kokkos::deep_copy(pc->R.getView(), R_host);
            Kokkos::deep_copy(pc->P.getView(), P_host);
            Kokkos::fence();
        }

        std::pair<size_t, size_t> runDecay(
                size_t nPart, double pz, int seed, size_t containerIndex, double tau0, double dt) {
            auto pc = makeContainer();
            createParticles(pc, nPart, pz);
            const size_t before = pc->getTotalNum();

            Options::seed = seed;
            MuonDecay decay(tau0, containerIndex, Physics::m_mu);
            const size_t destroyed = decay.apply(*pc, dt, 0, containerIndex);
            const size_t after     = pc->getTotalNum();

            EXPECT_EQ(before, after + destroyed);
            return {destroyed, after};
        }

    private:
        int oldSeed_m             = -1;
        bool oldUseQMAttributes_m = false;
    };

    // =====================================================================
    // Original destroy-only tests
    // =====================================================================

    TEST_F(DecayTest, ApplyReturnsZeroForNonPositiveDt) {
        auto pc = makeContainer();
        createParticles(pc, 64, 0.0);
        const size_t before = pc->getTotalNum();

        MuonDecay decay(1.0, 0, Physics::m_mu);
        EXPECT_EQ(decay.apply(*pc, 0.0, 7, 2), 0u);
        EXPECT_EQ(decay.apply(*pc, -1.0, 7, 2), 0u);
        EXPECT_EQ(pc->getTotalNum(), before);
    }

    TEST_F(DecayTest, ApplyReturnsZeroForInvalidLifetime) {
        const std::array<double, 4> invalidTau = {
                0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::infinity()};

        for (double tau0 : invalidTau) {
            auto pc = makeContainer();
            createParticles(pc, 32, 0.0);
            const size_t before = pc->getTotalNum();

            MuonDecay decay(tau0, 1, Physics::m_mu);
            EXPECT_EQ(decay.apply(*pc, 1.0, 0, 1), 0u);
            EXPECT_EQ(pc->getTotalNum(), before);
        }
    }

    TEST_F(DecayTest, ApplyReturnsZeroForEmptyContainer) {
        auto pc = makeContainer();
        ASSERT_EQ(pc->getLocalNum(), 0u);

        MuonDecay decay(1.0, 0, Physics::m_mu);
        EXPECT_EQ(decay.apply(*pc, 1.0, 0, 0), 0u);
        EXPECT_EQ(pc->getTotalNum(), 0u);
    }

    TEST_F(DecayTest, LargeDtDecaysAllRestParticles) {
        auto [destroyed, remaining] = runDecay(256, 0.0, 1234, 0, 1.0e-12, 1.0);

        EXPECT_EQ(destroyed, 256u);
        EXPECT_EQ(remaining, 0u);
    }

    TEST_F(DecayTest, HigherMomentumDecaysFewerParticles) {
        constexpr size_t nPart                            = 2000;
        auto [lowMomentumDestroyed, lowMomentumRemaining] = runDecay(nPart, 0.0, 77, 3, 1.0, 1.0);
        auto [highMomentumDestroyed, highMomentumRemaining] =
                runDecay(nPart, 10.0, 77, 3, 1.0, 1.0);

        EXPECT_GT(lowMomentumDestroyed, highMomentumDestroyed);
        EXPECT_GT(lowMomentumDestroyed, nPart / 3);
        EXPECT_LT(highMomentumDestroyed, nPart / 3);
        EXPECT_LT(lowMomentumRemaining, highMomentumRemaining);
    }

    TEST_F(DecayTest, SameSeedAndInputsAreReproducible) {
        auto [firstDestroyed, firstRemaining]   = runDecay(512, 0.5, 98765, 4, 2.0, 0.5);
        auto [secondDestroyed, secondRemaining] = runDecay(512, 0.5, 98765, 4, 2.0, 0.5);

        EXPECT_EQ(firstDestroyed, secondDestroyed);
        EXPECT_EQ(firstRemaining, secondRemaining);
    }

    // =====================================================================
    // Destroy-only mode preserved when no daughter container is set
    // =====================================================================

    TEST_F(DecayTest, NoDaughterContainerStillDestroysOnly) {
        auto muons = makeContainer();
        muons->setM(Physics::m_mu);
        createParticles(muons, 128, 0.0);

        Options::seed = 42;
        MuonDecay decay(1.0e-12, 0, Physics::m_mu);
        // No setDaughterContainer call.
        const size_t destroyed = decay.apply(*muons, 1.0, 0, 0);
        EXPECT_EQ(destroyed, 128u);
        EXPECT_EQ(muons->getTotalNum(), 0u);
    }

    // =====================================================================
    // Multi-container daughter generation tests (muon -> electron)
    // =====================================================================

    TEST_F(DecayTest, DaughterContainerReceivesDecayedParticles) {
        auto muons     = makeContainer();
        auto electrons = makeContainer(ParticleType::ELECTRON);
        muons->setM(Physics::m_mu);
        electrons->setM(Physics::m_e);
        createParticles(muons, 256, 0.0);

        Options::seed = 1234;
        MuonDecay decay(1.0e-12, 0, Physics::m_mu);
        decay.setDaughterContainer(electrons, Physics::m_e);

        const size_t destroyed = decay.apply(*muons, 1.0, 0, 0);
        EXPECT_EQ(destroyed, 256u);
        EXPECT_EQ(muons->getTotalNum(), 0u);
        EXPECT_EQ(electrons->getTotalNum(), 256u);
    }

    TEST_F(DecayTest, DaughterPositionMatchesParent) {
        constexpr size_t nPart = 64;
        auto muons             = makeContainer();
        auto electrons         = makeContainer(ParticleType::ELECTRON);
        muons->setM(Physics::m_mu);
        electrons->setM(Physics::m_e);
        createParticlesWithPositions(muons, nPart, 0.0);

        // Copy muon positions before decay.
        auto muR_host = muons->R.getHostMirror();
        Kokkos::deep_copy(muR_host, muons->R.getView());

        Options::seed = 999;
        MuonDecay decay(1.0e-12, 0, Physics::m_mu);
        decay.setDaughterContainer(electrons, Physics::m_e);

        const size_t destroyed = decay.apply(*muons, 1.0, 0, 0);
        ASSERT_EQ(destroyed, nPart);
        ASSERT_EQ(electrons->getLocalNum(), nPart);

        auto eR_host = electrons->R.getHostMirror();
        Kokkos::deep_copy(eR_host, electrons->R.getView());

        for (size_t i = 0; i < nPart; ++i) {
            EXPECT_DOUBLE_EQ(eR_host(i)[0], muR_host(i)[0]);
            EXPECT_DOUBLE_EQ(eR_host(i)[1], muR_host(i)[1]);
            EXPECT_DOUBLE_EQ(eR_host(i)[2], muR_host(i)[2]);
        }
    }

    TEST_F(DecayTest, DaughterMomentumIsPhysicalForRestMuons) {
        constexpr size_t nPart = 1024;
        auto muons             = makeContainer();
        auto electrons         = makeContainer(ParticleType::ELECTRON);
        muons->setM(Physics::m_mu);
        electrons->setM(Physics::m_e);
        createParticles(muons, nPart, 0.0);  // muons at rest

        Options::seed = 555;
        MuonDecay decay(1.0e-12, 0, Physics::m_mu);
        decay.setDaughterContainer(electrons, Physics::m_e);

        const size_t destroyed = decay.apply(*muons, 1.0, 0, 0);
        ASSERT_EQ(destroyed, nPart);
        ASSERT_EQ(electrons->getLocalNum(), nPart);

        auto eP_host = electrons->P.getHostMirror();
        Kokkos::deep_copy(eP_host, electrons->P.getView());

        const double eMax = Physics::MuonDecay::maxElectronEnergy();  // m_mu / 2

        for (size_t i = 0; i < nPart; ++i) {
            // P is stored as beta*gamma. Physical momentum = P * m_e.
            const double bg2    = eP_host(i)[0] * eP_host(i)[0] + eP_host(i)[1] * eP_host(i)[1]
                                  + eP_host(i)[2] * eP_host(i)[2];
            const double gamma  = std::sqrt(1.0 + bg2);
            const double energy = gamma * Physics::m_e;  // [GeV]

            // Electron energy must be between m_e and m_mu/2.
            EXPECT_GE(energy, Physics::m_e - 1.0e-12);
            EXPECT_LE(energy, eMax + 1.0e-12);

            // Verify E^2 = p^2 + m_e^2 (consistency of beta*gamma storage).
            const double p2 = bg2 * Physics::m_e * Physics::m_e;
            EXPECT_NEAR(energy * energy, p2 + Physics::m_e * Physics::m_e, 1.0e-15);
        }
    }

    TEST_F(DecayTest, BoostedMuonsProduceBoostedElectrons) {
        constexpr size_t nPart = 512;
        auto muons             = makeContainer();
        auto electrons         = makeContainer(ParticleType::ELECTRON);
        muons->setM(Physics::m_mu);
        electrons->setM(Physics::m_e);
        createParticles(muons, nPart, 5.0);  // boosted muons (pz = 5 in beta*gamma)

        Options::seed = 777;
        MuonDecay decay(1.0e-6, 0, Physics::m_mu);
        decay.setDaughterContainer(electrons, Physics::m_e);

        const size_t destroyed = decay.apply(*muons, 10.0, 0, 0);
        ASSERT_GT(destroyed, 0u);
        ASSERT_EQ(electrons->getLocalNum(), destroyed);

        auto eP_host = electrons->P.getHostMirror();
        Kokkos::deep_copy(eP_host, electrons->P.getView());

        for (size_t i = 0; i < destroyed; ++i) {
            const double bg2    = eP_host(i)[0] * eP_host(i)[0] + eP_host(i)[1] * eP_host(i)[1]
                                  + eP_host(i)[2] * eP_host(i)[2];
            const double gamma  = std::sqrt(1.0 + bg2);
            const double energy = gamma * Physics::m_e;

            // E^2 = p^2 + m_e^2 must hold.
            const double p2 = bg2 * Physics::m_e * Physics::m_e;
            EXPECT_NEAR(energy * energy, p2 + Physics::m_e * Physics::m_e, 1.0e-15);

            // Electron energy should be positive and finite.
            EXPECT_GT(energy, 0.0);
            EXPECT_TRUE(std::isfinite(energy));
        }
    }

    TEST_F(DecayTest, DaughterReproducibleWithSameSeed) {
        auto runWithDaughter = [&](int seed) {
            auto muons     = makeContainer();
            auto electrons = makeContainer(ParticleType::ELECTRON);
            muons->setM(Physics::m_mu);
            electrons->setM(Physics::m_e);
            createParticles(muons, 256, 1.0);

            Options::seed = seed;
            MuonDecay decay(1.0e-6, 0, Physics::m_mu);
            decay.setDaughterContainer(electrons, Physics::m_e);
            decay.apply(*muons, 1.0, 0, 0);

            // Read back electron momenta.
            auto eP_host = electrons->P.getHostMirror();
            Kokkos::deep_copy(eP_host, electrons->P.getView());

            std::vector<std::array<double, 3>> momenta(electrons->getLocalNum());
            for (size_t i = 0; i < momenta.size(); ++i) {
                momenta[i] = {eP_host(i)[0], eP_host(i)[1], eP_host(i)[2]};
            }
            return momenta;
        };

        auto first  = runWithDaughter(12345);
        auto second = runWithDaughter(12345);

        ASSERT_EQ(first.size(), second.size());
        for (size_t i = 0; i < first.size(); ++i) {
            EXPECT_DOUBLE_EQ(first[i][0], second[i][0]);
            EXPECT_DOUBLE_EQ(first[i][1], second[i][1]);
            EXPECT_DOUBLE_EQ(first[i][2], second[i][2]);
        }
    }

    TEST_F(DecayTest, PartialDecayCreatesDaughtersOnlyForDecayed) {
        constexpr size_t nPart = 1000;
        auto muons             = makeContainer();
        auto electrons         = makeContainer(ParticleType::ELECTRON);
        muons->setM(Physics::m_mu);
        electrons->setM(Physics::m_e);
        createParticles(muons, nPart, 0.5);

        Options::seed = 303;
        MuonDecay decay(2.0, 0, Physics::m_mu);  // long lifetime, moderate dt -> partial decay
        decay.setDaughterContainer(electrons, Physics::m_e);

        const size_t destroyed = decay.apply(*muons, 0.5, 0, 0);
        EXPECT_GT(destroyed, 0u);
        EXPECT_LT(destroyed, nPart);
        EXPECT_EQ(muons->getTotalNum(), nPart - destroyed);
        EXPECT_EQ(electrons->getTotalNum(), destroyed);
    }

    // =====================================================================
    // Pion decay tests (2-body: pi -> mu + nu)
    // =====================================================================

    TEST_F(DecayTest, PionDecayDaughterCountMatchesDestroyed) {
        auto pions = makeContainer();
        auto muons = makeContainer(ParticleType::MUON);
        pions->setM(Physics::m_pi);
        muons->setM(Physics::m_mu);
        createParticles(pions, 256, 0.0);

        Options::seed = 1234;
        PionDecay decay(1.0e-12, 0, Physics::m_pi);  // very short lifetime
        decay.setDaughterContainer(muons, Physics::m_mu);

        const size_t destroyed = decay.apply(*pions, 1.0, 0, 0);
        EXPECT_EQ(destroyed, 256u);
        EXPECT_EQ(pions->getTotalNum(), 0u);
        EXPECT_EQ(muons->getTotalNum(), 256u);
    }

    TEST_F(DecayTest, PionDecayTwoBodyFixedMomentum) {
        constexpr size_t nPart = 512;
        auto pions             = makeContainer();
        auto muons             = makeContainer(ParticleType::MUON);
        pions->setM(Physics::m_pi);
        muons->setM(Physics::m_mu);
        createParticles(pions, nPart, 0.0);  // pions at rest

        Options::seed = 888;
        PionDecay decay(1.0e-12, 0, Physics::m_pi);
        decay.setDaughterContainer(muons, Physics::m_mu);

        const size_t destroyed = decay.apply(*pions, 1.0, 0, 0);
        ASSERT_EQ(destroyed, nPart);
        ASSERT_EQ(muons->getLocalNum(), nPart);

        auto muP_host = muons->P.getHostMirror();
        Kokkos::deep_copy(muP_host, muons->P.getView());

        // Expected fixed momentum in pion rest frame: p = (m_pi^2 - m_mu^2) / (2 m_pi)
        const double pExpected = (Physics::m_pi * Physics::m_pi - Physics::m_mu * Physics::m_mu)
                                 / (2.0 * Physics::m_pi);
        // In beta*gamma units: bg = p / m_mu
        const double bgExpected = pExpected / Physics::m_mu;

        for (size_t i = 0; i < nPart; ++i) {
            const double bg2 = muP_host(i)[0] * muP_host(i)[0] + muP_host(i)[1] * muP_host(i)[1]
                               + muP_host(i)[2] * muP_host(i)[2];
            const double bg  = std::sqrt(bg2);

            // All daughters must have the same |beta*gamma| (monochromatic).
            EXPECT_NEAR(bg, bgExpected, 1.0e-12);
        }

        // Verify directions are not all identical (isotropic).
        bool foundDifferent = false;
        for (size_t i = 1; i < nPart; ++i) {
            if (std::abs(muP_host(i)[0] - muP_host(0)[0]) > 1.0e-14
                || std::abs(muP_host(i)[1] - muP_host(0)[1]) > 1.0e-14) {
                foundDifferent = true;
                break;
            }
        }
        EXPECT_TRUE(foundDifferent);
    }

    TEST_F(DecayTest, PionDecayBoostedConservesEnergy) {
        constexpr size_t nPart = 512;
        auto pions             = makeContainer();
        auto muons             = makeContainer(ParticleType::MUON);
        pions->setM(Physics::m_pi);
        muons->setM(Physics::m_mu);
        createParticles(pions, nPart, 3.0);  // boosted pions

        Options::seed = 999;
        PionDecay decay(1.0e-8, 0, Physics::m_pi);
        decay.setDaughterContainer(muons, Physics::m_mu);

        const size_t destroyed = decay.apply(*pions, 10.0, 0, 0);
        ASSERT_GT(destroyed, 0u);
        ASSERT_EQ(muons->getLocalNum(), destroyed);

        auto muP_host = muons->P.getHostMirror();
        Kokkos::deep_copy(muP_host, muons->P.getView());

        for (size_t i = 0; i < destroyed; ++i) {
            const double bg2    = muP_host(i)[0] * muP_host(i)[0] + muP_host(i)[1] * muP_host(i)[1]
                                  + muP_host(i)[2] * muP_host(i)[2];
            const double gamma  = std::sqrt(1.0 + bg2);
            const double energy = gamma * Physics::m_mu;

            // E^2 = p^2 + m_mu^2 must hold.
            const double p2 = bg2 * Physics::m_mu * Physics::m_mu;
            EXPECT_NEAR(energy * energy, p2 + Physics::m_mu * Physics::m_mu, 1.0e-15);

            EXPECT_GT(energy, 0.0);
            EXPECT_TRUE(std::isfinite(energy));
        }
    }

    // =====================================================================
    // Species validation in setDaughterContainer
    // =====================================================================

    TEST_F(DecayTest, SetDaughterContainerAcceptsCorrectSpecies) {
        auto electrons = makeContainer(ParticleType::ELECTRON);
        MuonDecay muonDecay(1.0e-6, 0, Physics::m_mu);
        EXPECT_NO_THROW(muonDecay.setDaughterContainer(electrons, Physics::m_e));

        auto muons = makeContainer(ParticleType::MUON);
        PionDecay pionDecay(1.0e-8, 0, Physics::m_pi);
        EXPECT_NO_THROW(pionDecay.setDaughterContainer(muons, Physics::m_mu));
    }

    TEST_F(DecayTest, SetDaughterContainerRejectsWrongSpecies) {
        MuonDecay decay(1.0e-6, 0, Physics::m_mu);

        auto protonDaughter = makeContainer(ParticleType::PROTON);
        EXPECT_THROW(decay.setDaughterContainer(protonDaughter, Physics::m_p), OpalException);

        auto photonDaughter = makeContainer(ParticleType::PHOTON);
        EXPECT_THROW(decay.setDaughterContainer(photonDaughter, 0.0), OpalException);

        auto muonDaughter = makeContainer(ParticleType::MUON);
        EXPECT_THROW(decay.setDaughterContainer(muonDaughter, Physics::m_mu), OpalException);
    }

    TEST_F(DecayTest, SetDaughterContainerRejectsUnnamedSpecies) {
        // Default-constructed container has Sp = 0 (PHOTON) via makeContainer()
        // with UNNAMED->0 cast; use explicit UNNAMED here which is -1 (sentinel).
        auto unnamed = makeContainer(ParticleType::UNNAMED);
        MuonDecay decay(1.0e-6, 0, Physics::m_mu);
        EXPECT_THROW(decay.setDaughterContainer(unnamed, Physics::m_e), OpalException);
    }

    TEST_F(DecayTest, PionDecayRejectsElectronDaughter) {
        auto electrons = makeContainer(ParticleType::ELECTRON);
        PionDecay decay(1.0e-8, 0, Physics::m_pi);
        // PionDecay's allowed daughter is MUON, not ELECTRON.
        EXPECT_THROW(decay.setDaughterContainer(electrons, Physics::m_e), OpalException);
    }

    TEST_F(DecayTest, SetDaughterContainerAcceptsNullForDestroyOnly) {
        auto muons = makeContainer();
        muons->setM(Physics::m_mu);
        createParticles(muons, 16, 0.0);

        MuonDecay decay(1.0e-12, 0, Physics::m_mu);
        // Null daughter = destroy-only mode; must not throw regardless of species.
        EXPECT_NO_THROW(decay.setDaughterContainer(nullptr, 0.0));

        // Decay still runs in destroy-only mode.
        Options::seed          = 11;
        const size_t destroyed = decay.apply(*muons, 1.0, 0, 0);
        EXPECT_EQ(destroyed, 16u);
    }

}  // namespace
