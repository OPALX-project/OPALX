/**
 * @file TestMuonDecaySpectrum.cpp
 * @brief Distribution-level test of MuonDecay against the analytic Michel spectrum.
 *
 * Pulls samples through the full MuonDecay::createDaughterParticles pipeline
 * (ParticleContainer + Kokkos pool seeded from Options::seed), reads the
 * daughter electron momenta back to the host, and compares the reduced-energy
 * histogram against the closed-form Michel PDF
 *   f(x) = 2 x^2 (3 - 2x) / Z,    Z = (2 x^3 - x^4) evaluated on [x_min, 1].
 *
 * Writes muon_michel_spectrum.csv into the test working directory; that CSV
 * is consumed by tools/spectrum_plots/plot_spectrum.py.
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
#include "Physics/MuonDecay.h"
#include "Physics/ParticleProperties.h"
#include "Physics/Physics.h"
#include "Processes/GlobalProcesses/MuonDecay.h"
#include "SpectrumTestSupport.h"
#include "Utilities/Options.h"

namespace {

    using PC_t = ParticleContainer<double, 3>;

    class MuonDecaySpectrumTest : public ::testing::Test {
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

        void createParticlesAtRest(std::shared_ptr<PC_t>& pc, std::size_t nPart) {
            pc->createParticles(nPart);
            auto R_host = pc->R.getHostMirror();
            auto P_host = pc->P.getHostMirror();
            for (std::size_t i = 0; i < nPart; ++i) {
                R_host(i)[0] = 0.0;
                R_host(i)[1] = 0.0;
                R_host(i)[2] = 0.0;
                P_host(i)[0] = 0.0;
                P_host(i)[1] = 0.0;
                P_host(i)[2] = 0.0;
            }
            Kokkos::deep_copy(pc->R.getView(), R_host);
            Kokkos::deep_copy(pc->P.getView(), P_host);
            Kokkos::fence();
        }

    private:
        int oldSeed_m             = -1;
        bool oldUseQMAttributes_m = false;
    };

    /// Closed-form normalization of the Michel spectrum on [a, 1]:
    ///   Z(a) = integral_a^1 2 x^2 (3 - 2 x) dx = (2 x^3 - x^4)|_a^1 = 1 - 2 a^3 + a^4.
    double michelNormalization(double a) { return 1.0 - 2.0 * a * a * a + a * a * a * a; }

    /// Closed-form mean of the Michel spectrum on [a, 1]:
    ///   E[x] = integral_a^1 x f(x) dx / Z = ((6/4) x^4 - (4/5) x^5)|_a^1 / Z
    double michelMean(double a) {
        const double a4  = a * a * a * a;
        const double a5  = a4 * a;
        const double num = (1.5 - 0.8) - (1.5 * a4 - 0.8 * a5);
        return num / michelNormalization(a);
    }

    TEST_F(MuonDecaySpectrumTest, MuonDecaySpectrumMatchesMichel) {
        constexpr std::size_t nMuons = 100000;
        constexpr std::size_t nBins  = 50;

        auto muons     = makeContainer();
        auto electrons = makeContainer(ParticleType::ELECTRON);
        muons->setM(Physics::m_mu);
        electrons->setM(Physics::m_e);
        createParticlesAtRest(muons, nMuons);

        Options::seed = 20260515;
        /**
         * Intentionally chosen timestep so all particles decay.
         * We are testing the daughter particle energy spectrum, not the decay
         * time.
         */
        MuonDecay decay(1.0e-12, 0, Physics::m_mu);
        decay.setDaughterContainer(electrons, Physics::m_e);

        const std::size_t marked    = decay.apply(*muons, 1.0, 0, 0);
        const std::size_t destroyed = muons->deleteInvalidParticles();
        ASSERT_EQ(marked, nMuons);
        ASSERT_EQ(destroyed, nMuons);
        ASSERT_EQ(electrons->getLocalNum(), nMuons);

        auto eP_host = electrons->P.getHostMirror();
        Kokkos::deep_copy(eP_host, electrons->P.getView());

        const double eMax = Physics::MuonDecay::maxElectronEnergy();  // m_mu / 2
        const double xMin = Physics::MuonDecay::minElectronX();

        std::vector<double> xSamples;
        xSamples.reserve(nMuons);
        for (std::size_t i = 0; i < nMuons; ++i) {
            const double bg2 = eP_host(i)[0] * eP_host(i)[0] + eP_host(i)[1] * eP_host(i)[1]
                               + eP_host(i)[2] * eP_host(i)[2];
            const double energy = std::sqrt(1.0 + bg2) * Physics::m_e;
            xSamples.push_back(energy / eMax);
        }

        const auto hist = opalx::test::makeHistogram(xSamples, xMin, 1.0, nBins);

        const double Z = michelNormalization(xMin);
        std::vector<double> analyticPdf(nBins, 0.0);
        for (std::size_t i = 0; i < nBins; ++i) {
            const double x = hist.center(i);
            analyticPdf[i] = Physics::MuonDecay::michelSpectrum(x) / Z;
        }

        const double l1      = opalx::test::l1Distance(hist, analyticPdf);
        const double meanX   = opalx::test::sampleMean(xSamples);
        const double area    = opalx::test::histogramArea(hist);
        const double meanRef = michelMean(xMin);

        std::cout << "[MuonDecaySpectrum] L1=" << l1 << " mean=" << meanX
                  << " (analytic=" << meanRef << ") area=" << area << '\n';

        EXPECT_LT(l1, 0.10);
        EXPECT_NEAR(meanX, meanRef, meanRef * 0.02);
        EXPECT_NEAR(area, 1.0, 0.01);

        const std::string csvPath = "muon_michel_spectrum.csv";
        opalx::test::writeSpectrumCsv(csvPath, hist, analyticPdf, "x = E_e / (m_mu/2)");
        std::cout << "[MuonDecaySpectrum] CSV written: " << csvPath << '\n';
    }

}  // namespace
