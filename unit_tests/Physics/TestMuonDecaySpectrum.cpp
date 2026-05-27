/**
 * @file TestMuonDecaySpectrum.cpp
 * @brief Distribution-level test of MuonDecay against the polarized V-A angular spectrum.
 *
 * Fires fully-polarized muons (|P| = 1, spin along +z) at rest through the full
 * MuonDecay::createDaughterParticles pipeline (ParticleContainer + Kokkos pool
 * seeded from Options::seed), reads the daughter electron momenta back to the
 * host, and compares the cos(theta) distribution (theta measured from the muon
 * spin) against the analytic angular marginal
 *   p(c) = 1/2 (1 + A c),   A = s |P| N/Z,
 *   N = \int_a^1 2 x^2 (2x - 1) dx,   Z = \int_a^1 2 x^2 (3 - 2x) dx,   a = x_min,
 * where s = asymSign is +1 for mu- and -1 for mu+. The energy (x) marginal is the
 * polarization-independent Michel spectrum and is covered elsewhere.
 *
 * Writes muon_angular_spectrum_mu_minus.csv and muon_angular_spectrum_mu_plus.csv
 * into the test working directory; those CSVs are consumed by
 * tools/spectrum_plots/plot_spectrum.py (use --overlay to compare the two).
 */

#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
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

    /// Closed-form normalization of the Michel spectrum on [a, 1] (denominator Z):
    ///   Z(a) = integral_a^1 2 x^2 (3 - 2 x) dx = 1 - 2 a^3 + a^4.
    double michelNormalization(double a) { return 1.0 - 2.0 * a * a * a + a * a * a * a; }

    /// Numerator of the integrated angular asymmetry on [a, 1]:
    ///   N(a) = integral_a^1 2 x^2 (2 x - 1) dx = 1/3 - a^4 + (2/3) a^3.
    double angularAsymmetryNumerator(double a) {
        const double a3 = a * a * a;
        const double a4 = a3 * a;
        return 1.0 / 3.0 - a4 + (2.0 / 3.0) * a3;
    }

    /// Integrated asymmetry coefficient A/(s|P|) = N(a)/Z(a) (~1/3 for a -> 0).
    double integratedAsymmetry(double a) {
        return angularAsymmetryNumerator(a) / michelNormalization(a);
    }

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

        std::shared_ptr<PC_t> makeContainer(
                ParticleType species = ParticleType::UNNAMED, bool spinEnabled = false) {
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
            auto pc = std::make_shared<PC_t>(mesh, fl, spinEnabled);
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

        /// Set the same polarization vector on every particle in the container.
        void setUniformPolarization(std::shared_ptr<PC_t>& pc, const std::array<float, 3>& pol) {
            using spin_t        = typename PC_t::spin_vector_type;
            auto Pol_host       = pc->Pol.getHostMirror();
            const std::size_t n = pc->getLocalNum();
            for (std::size_t i = 0; i < n; ++i) {
                Pol_host(i) = spin_t{pol[0], pol[1], pol[2]};
            }
            Kokkos::deep_copy(pc->Pol.getView(), Pol_host);
            Kokkos::fence();
        }

        /// Run one polarized muon-decay case and validate the cos(theta) spectrum.
        /// parentSign: -1 for mu-, +1 for mu+.
        void runAngularCase(int parentSign, const std::string& csvPath) {
            constexpr std::size_t nMuons = 100000;
            constexpr std::size_t nBins  = 50;
            constexpr double polMag      = 1.0;  // fully polarized along +z

            auto muons     = makeContainer(ParticleType::UNNAMED, /*spinEnabled=*/true);
            auto electrons = makeContainer(ParticleType::ELECTRON);
            muons->setM(Physics::m_mu);
            electrons->setM(Physics::m_e);
            createParticlesAtRest(muons, nMuons);
            setUniformPolarization(muons, {0.0f, 0.0f, static_cast<float>(polMag)});

            Options::seed = 20260515;
            // Timestep chosen so all muons decay within one step; we test the angular
            // spectrum of the daughters, not the decay time.
            MuonDecay decay(1.0e-12, 0, Physics::m_mu, parentSign);
            decay.setDaughterContainer(electrons, Physics::m_e);

            const std::size_t marked    = decay.apply(*muons, 1.0, 0, 0);
            const std::size_t destroyed = muons->deleteInvalidParticles();
            ASSERT_EQ(marked, nMuons);
            ASSERT_EQ(destroyed, nMuons);
            ASSERT_EQ(electrons->getLocalNum(), nMuons);

            auto eP_host = electrons->P.getHostMirror();
            Kokkos::deep_copy(eP_host, electrons->P.getView());

            // Muons are at rest, so the lab momentum direction equals the rest-frame
            // direction; with spin along +z, cos(theta) = P_z / |P|.
            std::vector<double> cosSamples;
            cosSamples.reserve(nMuons);
            for (std::size_t i = 0; i < nMuons; ++i) {
                const double px  = eP_host(i)[0];
                const double py  = eP_host(i)[1];
                const double pz  = eP_host(i)[2];
                const double mag = std::sqrt(px * px + py * py + pz * pz);
                if (mag > 0.0) {
                    cosSamples.push_back(pz / mag);
                }
            }

            const double xMin     = Physics::MuonDecay::minElectronX();
            const double asymSign = (parentSign < 0) ? +1.0 : -1.0;
            const double A        = asymSign * polMag * integratedAsymmetry(xMin);

            const auto hist = opalx::test::makeHistogram(cosSamples, -1.0, 1.0, nBins);
            std::vector<double> analyticPdf(nBins, 0.0);
            for (std::size_t i = 0; i < nBins; ++i) {
                analyticPdf[i] = 0.5 * (1.0 + A * hist.center(i));
            }

            const double l1      = opalx::test::l1Distance(hist, analyticPdf);
            const double meanCos = opalx::test::sampleMean(cosSamples);
            const double area    = opalx::test::histogramArea(hist);
            const double meanRef = A / 3.0;  // <cos> for p(c) = 1/2 (1 + A c) on [-1, 1]

            std::cout << "[MuonAngularSpectrum sign=" << parentSign << "] L1=" << l1
                      << " mean=" << meanCos << " (analytic=" << meanRef << ") area=" << area
                      << '\n';

            EXPECT_LT(l1, 0.10);
            EXPECT_NEAR(meanCos, meanRef, 0.01);
            EXPECT_NEAR(area, 1.0, 0.01);

            opalx::test::writeSpectrumCsv(
                    csvPath, hist, analyticPdf, "cos(theta) w.r.t. muon spin");
            std::cout << "[MuonAngularSpectrum] CSV written: " << csvPath << '\n';
        }

    private:
        int oldSeed_m             = -1;
        bool oldUseQMAttributes_m = false;
    };

    TEST_F(MuonDecaySpectrumTest, MuMinusAngularSpectrumMatchesPolarization) {
        runAngularCase(-1, "muon_angular_spectrum_mu_minus.csv");
    }

    TEST_F(MuonDecaySpectrumTest, MuPlusAngularSpectrumMatchesPolarization) {
        runAngularCase(+1, "muon_angular_spectrum_mu_plus.csv");
    }

}  // namespace
