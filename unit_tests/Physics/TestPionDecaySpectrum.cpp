/**
 * @file TestPionDecaySpectrum.cpp
 * @brief Distribution-level tests of PionDecay (deterministic two-body kinematics).
 *
 * Drives PionDecay::createDaughterParticles through the full ParticleContainer
 * + Kokkos pool pipeline. Two checks:
 *
 *  1. PionDecayRestFrameMomentumIsFixed:
 *     With pions at rest, every daughter muon must have |p| = (M_pi^2 - m_mu^2) / (2 M_pi),
 *     and the cos(theta) distribution of the decay direction must be uniform on [-1, 1].
 *
 *  2. PionDecayBoostedLabEnergyIsFlatBox:
 *     With boosted pions (pz = 5 in beta*gamma units), the lab-frame muon energy
 *     spectrum is uniform on [E_-, E_+] where
 *        E_+- = gamma * E* +- beta*gamma * p*,
 *        E*   = (M_pi^2 + m_mu^2) / (2 M_pi),
 *        p*   = (M_pi^2 - m_mu^2) / (2 M_pi).
 *
 * CSVs (consumed by tools/spectrum_plots/plot_spectrum.py):
 *   - pion_rest_costheta.csv
 *   - pion_boosted_energy_box.csv
 */

#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

#include "Ippl.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Physics/ParticleProperties.h"
#include "Physics/Physics.h"
#include "Processes/GlobalProcesses/PionDecay.h"
#include "SpectrumTestSupport.h"
#include "Utilities/Options.h"

namespace {

    using PC_t = ParticleContainer<double, 3>;

    class PionDecaySpectrumTest : public ::testing::Test {
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
            pc->setBunchStateHandler(std::make_shared<BunchStateHandler>());
            return pc;
        }

        void createParticles(std::shared_ptr<PC_t>& pc, std::size_t nPart, double pz) {
            pc->createParticles(nPart);
            auto R_host = pc->R.getHostMirror();
            auto P_host = pc->P.getHostMirror();
            for (std::size_t i = 0; i < nPart; ++i) {
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

    private:
        int oldSeed_m             = -1;
        bool oldUseQMAttributes_m = false;
    };

    TEST_F(PionDecaySpectrumTest, PionDecayRestFrameMomentumIsFixed) {
        constexpr std::size_t nPions = 100000;
        constexpr std::size_t nBins  = 50;

        auto pions = makeContainer();
        auto muons = makeContainer(ParticleType::MUON);
        pions->setM(Physics::m_pi);
        muons->setM(Physics::m_mu);
        createParticles(pions, nPions, 0.0);

        Options::seed = 20260516;
        PionDecay decay(1.0e-12, 0, Physics::m_pi);
        decay.setDaughterContainer(muons, Physics::m_mu);

        const std::size_t marked    = decay.apply(*pions, 1.0, 0, 0);
        const std::size_t destroyed = pions->deleteInvalidParticles();
        ASSERT_EQ(marked, nPions);
        ASSERT_EQ(destroyed, nPions);
        ASSERT_EQ(muons->getLocalNum(), nPions);

        auto muP_host = muons->P.getHostMirror();
        Kokkos::deep_copy(muP_host, muons->P.getView());

        const double pFixed = (Physics::m_pi * Physics::m_pi - Physics::m_mu * Physics::m_mu)
                              / (2.0 * Physics::m_pi);

        std::vector<double> cosThetaSamples;
        cosThetaSamples.reserve(nPions);
        double maxRelDev = 0.0;
        for (std::size_t i = 0; i < nPions; ++i) {
            const double bgx = muP_host(i)[0];
            const double bgy = muP_host(i)[1];
            const double bgz = muP_host(i)[2];
            const double pmag =
                    Physics::m_mu * std::sqrt(bgx * bgx + bgy * bgy + bgz * bgz);
            maxRelDev = std::max(maxRelDev, std::abs(pmag - pFixed) / pFixed);
            cosThetaSamples.push_back(Physics::m_mu * bgz / pmag);
        }
        std::cout << "[PionDecayRestFrame] |p|=" << pFixed
                  << " max rel dev=" << maxRelDev << '\n';
        EXPECT_LT(maxRelDev, 1.0e-12);

        const auto hist = opalx::test::makeHistogram(cosThetaSamples, -1.0, 1.0, nBins);
        const std::vector<double> uniform(nBins, 0.5);

        const double l1   = opalx::test::l1Distance(hist, uniform);
        const double mean = opalx::test::sampleMean(cosThetaSamples);
        const double area = opalx::test::histogramArea(hist);
        std::cout << "[PionDecayRestFrame] cos(theta) L1=" << l1 << " mean=" << mean
                  << " area=" << area << '\n';

        EXPECT_LT(l1, 0.10);
        EXPECT_NEAR(mean, 0.0, 0.02);
        EXPECT_NEAR(area, 1.0, 0.01);

        const std::string csvPath = "pion_rest_costheta.csv";
        opalx::test::writeSpectrumCsv(csvPath, hist, uniform, "cos(theta) (rest frame)");
        std::cout << "[PionDecayRestFrame] CSV written: " << csvPath << '\n';
    }

    TEST_F(PionDecaySpectrumTest, PionDecayBoostedLabEnergyIsFlatBox) {
        constexpr std::size_t nPions   = 100000;
        constexpr std::size_t nBins    = 50;
        constexpr double parentBetaGam = 5.0;  // pz of pion in beta*gamma units

        auto pions = makeContainer();
        auto muons = makeContainer(ParticleType::MUON);
        pions->setM(Physics::m_pi);
        muons->setM(Physics::m_mu);
        createParticles(pions, nPions, parentBetaGam);

        Options::seed = 20260517;
        PionDecay decay(1.0e-12, 0, Physics::m_pi);
        decay.setDaughterContainer(muons, Physics::m_mu);

        const std::size_t marked    = decay.apply(*pions, 1.0, 0, 0);
        const std::size_t destroyed = pions->deleteInvalidParticles();
        ASSERT_EQ(marked, nPions);
        ASSERT_EQ(destroyed, nPions);
        ASSERT_EQ(muons->getLocalNum(), nPions);

        auto muP_host = muons->P.getHostMirror();
        Kokkos::deep_copy(muP_host, muons->P.getView());

        // Rest-frame two-body kinematics.
        const double M     = Physics::m_pi;
        const double md    = Physics::m_mu;
        const double pStar = (M * M - md * md) / (2.0 * M);
        const double eStar = (M * M + md * md) / (2.0 * M);

        // Parent boost.
        const double gamma   = std::sqrt(1.0 + parentBetaGam * parentBetaGam);
        const double betaGam = parentBetaGam;  // = beta * gamma
        const double eMinus  = gamma * eStar - betaGam * pStar;
        const double ePlus   = gamma * eStar + betaGam * pStar;
        const double uniform = 1.0 / (ePlus - eMinus);

        std::vector<double> labEnergies;
        labEnergies.reserve(nPions);
        for (std::size_t i = 0; i < nPions; ++i) {
            const double bg2 = muP_host(i)[0] * muP_host(i)[0] + muP_host(i)[1] * muP_host(i)[1]
                               + muP_host(i)[2] * muP_host(i)[2];
            labEnergies.push_back(std::sqrt(1.0 + bg2) * Physics::m_mu);
        }

        const auto hist = opalx::test::makeHistogram(labEnergies, eMinus, ePlus, nBins);
        const std::vector<double> analyticPdf(nBins, uniform);

        const double l1      = opalx::test::l1Distance(hist, analyticPdf);
        const double mean    = opalx::test::sampleMean(labEnergies);
        const double meanRef = 0.5 * (eMinus + ePlus);
        const double area    = opalx::test::histogramArea(hist);
        std::cout << "[PionDecayBoosted] E in [" << eMinus << ", " << ePlus
                  << "] GeV, L1=" << l1 << " mean=" << mean
                  << " (analytic=" << meanRef << ") area=" << area << '\n';

        EXPECT_LT(l1, 0.10);
        EXPECT_NEAR(mean, meanRef, 0.005 * (ePlus - eMinus));
        EXPECT_NEAR(area, 1.0, 0.01);

        const std::string csvPath = "pion_boosted_energy_box.csv";
        opalx::test::writeSpectrumCsv(csvPath, hist, analyticPdf, "E_mu^lab [GeV]");
        std::cout << "[PionDecayBoosted] CSV written: " << csvPath << '\n';
    }

}  // namespace
