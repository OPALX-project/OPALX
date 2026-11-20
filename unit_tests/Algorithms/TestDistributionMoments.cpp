#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "Algorithms/DistributionMoments.h"
#include "Ippl.h"
#include "PartBunch/BunchStateHandler.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"

namespace {
    using Vec3 = Vector_t<double, 3>;

    struct ExpectedMoments {
        Vec3 meanR  = 0.0;
        Vec3 meanP  = 0.0;
        Vec3 stdR   = 0.0;
        Vec3 stdP   = 0.0;
        Vec3 covRP  = 0.0;
        Vec3 epsN   = 0.0;
        Vec3 epsGeo = 0.0;
        double meanGamma{0.0};
        double meanGammaZ{0.0};
        double meanEkinMeV{0.0};
        double stdEkinMeV_likeImplementation{0.0};  // matches current code path (no /Np)
        double Dx{0.0}, DDx{0.0}, Dy{0.0}, DDy{0.0};
    };

    static ExpectedMoments computeExpected(
            const std::vector<Vec3>& R, const std::vector<Vec3>& P, double massGeV) {
        ExpectedMoments e;
        const size_t N = R.size();
        if (N == 0) {
            return e;
        }

        // Means
        for (size_t i = 0; i < N; ++i) {
            for (int d = 0; d < 3; ++d) {
                e.meanR(d) += R[i](d);
                e.meanP(d) += P[i](d);
            }
            const double p2    = P[i](0) * P[i](0) + P[i](1) * P[i](1) + P[i](2) * P[i](2);
            const double gamma = std::sqrt(p2 + 1.0);
            e.meanGamma += gamma;
            e.meanEkinMeV += (gamma - 1.0) * massGeV * Units::GeV2MeV;
            e.meanGammaZ += P[i](2);
        }
        for (int d = 0; d < 3; ++d) {
            e.meanR(d) /= static_cast<double>(N);
            e.meanP(d) /= static_cast<double>(N);
        }
        e.meanGamma /= static_cast<double>(N);
        e.meanEkinMeV /= static_cast<double>(N);
        e.meanGammaZ /= static_cast<double>(N);
        e.meanGammaZ = std::sqrt(e.meanGammaZ * e.meanGammaZ + 1.0);

        // Central covariance of [x,px,y,py,z,pz]
        double cov[6][6]   = {};
        double ncent[6][6] = {};
        for (size_t i = 0; i < N; ++i) {
            const double part[6] = {
                    R[i](0) - e.meanR(0), P[i](0) - e.meanP(0), R[i](1) - e.meanR(1),
                    P[i](1) - e.meanP(1), R[i](2) - e.meanR(2), P[i](2) - e.meanP(2),
            };
            const double part_nc[6] = {
                    R[i](0), P[i](0), R[i](1), P[i](1), R[i](2), P[i](2),
            };
            for (int a = 0; a < 6; ++a) {
                for (int b = 0; b < 6; ++b) {
                    cov[a][b] += part[a] * part[b];
                    ncent[a][b] += part_nc[a] * part_nc[b];
                }
            }

            const double p2    = P[i](0) * P[i](0) + P[i](1) * P[i](1) + P[i](2) * P[i](2);
            const double gamma = std::sqrt(p2 + 1.0);
            const double ekin  = (gamma - 1.0) * massGeV * Units::GeV2MeV;
            e.stdEkinMeV_likeImplementation += (ekin - e.meanEkinMeV) * (ekin - e.meanEkinMeV);
        }
        for (int a = 0; a < 6; ++a) {
            for (int b = 0; b < 6; ++b) {
                cov[a][b] /= static_cast<double>(N);
                ncent[a][b] /= static_cast<double>(N);
            }
        }
        e.stdEkinMeV_likeImplementation =
                std::sqrt(e.stdEkinMeV_likeImplementation / static_cast<double>(N));

        // stdR/stdP from central cov
        for (int d = 0; d < 3; ++d) {
            e.stdR(d) = std::sqrt(std::max(cov[2 * d][2 * d], 0.0));
            e.stdP(d) = std::sqrt(std::max(cov[2 * d + 1][2 * d + 1], 0.0));
        }

        // cov(R,P) using the same formula as implementation:
        // sumRP = notCentMoments(2d,2d+1) - meanR*meanP
        for (int d = 0; d < 3; ++d) {
            e.covRP(d)           = ncent[2 * d][2 * d + 1] - e.meanR(d) * e.meanP(d);
            const double squared = std::pow(e.stdR(d) * e.stdP(d), 2) - std::pow(e.covRP(d), 2);
            e.epsN(d)            = std::sqrt(std::max(squared, 0.0));
        }
        const double betaGamma = std::sqrt(std::max(e.meanGamma * e.meanGamma - 1.0, 0.0));
        e.epsGeo               = e.epsN / Vec3(betaGamma);

        // Dispersion terms as exposed by getters (central moments against pz)
        e.Dx  = cov[0][5];
        e.DDx = cov[1][5];
        e.Dy  = cov[2][5];
        e.DDy = cov[3][5];

        return e;
    }
}  // namespace

class DistributionMomentsTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() { ippl::finalize(); }

    void SetUp() override { bunchStateHandler = std::make_shared<BunchStateHandler>(); }

    std::shared_ptr<BunchStateHandler> bunchStateHandler;

    void TearDown() override { bunchStateHandler.reset(); }
};

TEST_F(DistributionMomentsTest, CompactSymmetricMatchesDenseAcrossMPIAndMassModes) {
    const size_t rank = ippl::Comm->rank();
    // Empty globally, rank-zero-only ownership, then uneven ownership and a larger reduction.
    for (int ownership = 0; ownership < 4; ++ownership) {
        const size_t local = ownership == 0   ? 0
                             : ownership == 1 ? (rank == 0 ? 17 : 0)
                                              : (ownership == 2 ? 31 : 4097) + 7 * rank;
        size_t total       = local;
        ippl::Comm->allreduce(total, 1, std::plus<size_t>());
        for (bool attributes : {false, true}) {
            for (bool referenceMass : {false, true}) {
                Kokkos::View<Vec3*> r("compact-test-r", local), p("compact-test-p", local);
                Kokkos::View<double*> mass("compact-test-m", attributes ? local : 1);
                auto rh = Kokkos::create_mirror_view(r);
                auto ph = Kokkos::create_mirror_view(p);
                auto mh = Kokkos::create_mirror_view(mass);
                for (size_t i = 0; i < local; ++i) {
                    const double t = static_cast<double>(i + 19 * rank);
                    // Off-center, nonuniform data with cross correlations in all six coordinates.
                    rh(i) =
                            Vec3(3. + .02 * std::sin(t), -2. + .03 * std::cos(.7 * t),
                                 7. + .04 * std::sin(.3 * t));
                    ph(i) =
                            Vec3(.02 + .004 * std::cos(.4 * t), -.03 + .006 * std::sin(.9 * t),
                                 .4 + .01 * std::cos(.2 * t));
                    if (attributes) mh(i) = Physics::m_e * (1. + .05 * (i % 3));
                }
                if (!attributes) mh(0) = Physics::m_e;
                Kokkos::deep_copy(r, rh);
                Kokkos::deep_copy(p, ph);
                Kokkos::deep_copy(mass, mh);
                DistributionMoments dense, compact;
                auto denseState   = bunchStateHandler->registerContainer();
                auto compactState = bunchStateHandler->registerContainer();
                dense.setContainerState(denseState);
                compact.setContainerState(compactState);
                ASSERT_FALSE(dense.usesCompactSymmetricReduction());
                compact.setCompactSymmetricReduction(true);
                dense.setEnergyReferenceMass(Physics::m_p, referenceMass);
                compact.setEnergyReferenceMass(Physics::m_p, referenceMass);
                dense.computeMoments(r, p, mass, total, local);
                compact.computeMoments(r, p, mass, total, local);
                const auto check = [](double a, double b) {
                    if (std::isnan(a))
                        EXPECT_TRUE(std::isnan(b));
                    else
                        EXPECT_NEAR(a, b, 2.e-13 + 2.e-12 * std::abs(a));
                };
                const auto cm = compact.getMoments6x6();
                const auto dm = dense.getMoments6x6();
                for (int i = 0; i < 6; ++i) {
                    EXPECT_DOUBLE_EQ(compact.getMeans()[i], dense.getMeans()[i]);
                    for (int j = 0; j < 6; ++j) {
                        check(dm(i, j), cm(i, j));
                        EXPECT_DOUBLE_EQ(cm(i, j), cm(j, i));
                    }
                }
                for (int d = 0; d < 3; ++d) {
                    check(dense.getStandardDeviationPosition()[d],
                          compact.getStandardDeviationPosition()[d]);
                    check(dense.getStandardDeviationMomentum()[d],
                          compact.getStandardDeviationMomentum()[d]);
                    check(dense.getNormalizedEmittance()[d], compact.getNormalizedEmittance()[d]);
                    check(dense.getGeometricEmittance()[d], compact.getGeometricEmittance()[d]);
                    check(dense.getStandardDeviationRP()[d], compact.getStandardDeviationRP()[d]);
                    check(dense.getHalo()[d], compact.getHalo()[d]);
                    check(dense.getRestFrameTemperature()[d], compact.getRestFrameTemperature()[d]);
                }
                EXPECT_DOUBLE_EQ(dense.getMeanKineticEnergy(), compact.getMeanKineticEnergy());
                check(dense.getStdKineticEnergy(), compact.getStdKineticEnergy());
                check(dense.getDx(), compact.getDx());
                check(dense.getDDx(), compact.getDDx());
                check(dense.getDy(), compact.getDy());
                check(dense.getDDy(), compact.getDDy());
                EXPECT_FALSE(compactState->momentsDirty);
            }
        }
    }
}

TEST_F(DistributionMomentsTest, ComputeMoments_MeanStdEmittanceEnergyDispersion) {
    // Deterministic small data set (units: R in m, P in beta*gamma, M in GeV/c^2).
    const std::vector<Vec3> R = {
            Vec3(1.0e-3, 2.0e-3, -1.0e-3), Vec3(-2.0e-3, 0.5e-3, 3.0e-3),
            Vec3(0.0e-3, -1.5e-3, 2.0e-3), Vec3(1.5e-3, -0.5e-3, -2.5e-3),
            Vec3(-0.5e-3, 1.0e-3, 0.0e-3),
    };
    const std::vector<Vec3> P = {
            Vec3(0.02, -0.01, 0.10), Vec3(-0.03, 0.02, 0.08), Vec3(0.01, 0.01, 0.11),
            Vec3(0.00, -0.02, 0.09), Vec3(0.04, 0.00, 0.105),
    };
    ASSERT_EQ(R.size(), P.size());
    const size_t N = R.size();

    // Build device views that match DistributionMoments signatures.
    Kokkos::View<Vec3*> Rview_d("Rview_d", N);
    Kokkos::View<Vec3*> Pview_d("Pview_d", N);
    Kokkos::View<double*> Mview_d("Mview_d", 1);  // scalar mass mode

    auto Rview_h = Kokkos::create_mirror_view(Rview_d);
    auto Pview_h = Kokkos::create_mirror_view(Pview_d);
    auto Mview_h = Kokkos::create_mirror_view(Mview_d);

    for (size_t i = 0; i < N; ++i) {
        Rview_h(i) = R[i];
        Pview_h(i) = P[i];
    }
    Mview_h(0) = Physics::m_e;

    Kokkos::deep_copy(Rview_d, Rview_h);
    Kokkos::deep_copy(Pview_d, Pview_h);
    Kokkos::deep_copy(Mview_d, Mview_h);

    DistributionMoments dm;
    auto containerState = bunchStateHandler->registerContainer();
    dm.setContainerState(containerState);
    dm.computeMoments(Rview_d, Pview_d, Mview_d, /*Np=*/N, /*Nlocal=*/N);

    const auto exp = computeExpected(R, P, Physics::m_e);

    const Vec3 meanR = dm.getMeanPosition();
    const Vec3 meanP = dm.getMeanMomentum();
    const Vec3 stdR  = dm.getStandardDeviationPosition();
    const Vec3 stdP  = dm.getStandardDeviationMomentum();
    const Vec3 epsN  = dm.getNormalizedEmittance();
    const Vec3 epsG  = dm.getGeometricEmittance();

    for (int d = 0; d < 3; ++d) {
        EXPECT_NEAR(meanR(d), exp.meanR(d), 1e-15);
        EXPECT_NEAR(meanP(d), exp.meanP(d), 1e-15);

        EXPECT_NEAR(stdR(d), exp.stdR(d), 1e-15);
        EXPECT_NEAR(stdP(d), exp.stdP(d), 1e-15);

        EXPECT_NEAR(epsN(d), exp.epsN(d), 1e-15);
        EXPECT_NEAR(epsG(d), exp.epsGeo(d), 1e-15);
    }

    EXPECT_NEAR(dm.getMeanGamma(), exp.meanGamma, 1e-15);
    EXPECT_NEAR(dm.getMeanGammaZ(), exp.meanGammaZ, 1e-15);
    EXPECT_NEAR(dm.getMeanKineticEnergy(), exp.meanEkinMeV, 1e-15);
    EXPECT_NEAR(dm.getStdKineticEnergy(), exp.stdEkinMeV_likeImplementation, 1e-15);

    EXPECT_NEAR(dm.getDx(), exp.Dx, 1e-15);
    EXPECT_NEAR(dm.getDDx(), exp.DDx, 1e-15);
    EXPECT_NEAR(dm.getDy(), exp.Dy, 1e-15);
    EXPECT_NEAR(dm.getDDy(), exp.DDy, 1e-15);
}

namespace {
    double dot3(const Vec3& a, const Vec3& b) { return a(0) * b(0) + a(1) * b(1) + a(2) * b(2); }

    /// Boost a rest-frame momentum u' (beta*gamma) to the frame in which the rest frame moves
    /// with velocity beta.
    Vec3 boostToLab(const Vec3& uRest, const Vec3& beta) {
        const double gamma     = 1.0 / std::sqrt(1.0 - dot3(beta, beta));
        const double gammaRest = std::sqrt(1.0 + dot3(uRest, uRest));
        const double shift = gamma * gamma / (gamma + 1.0) * dot3(beta, uRest) + gamma * gammaRest;
        Vec3 u;
        for (int d = 0; d < 3; ++d) {
            u(d) = uRest(d) + shift * beta(d);
        }
        return u;
    }
}  // namespace

// An anisotropic set with zero total momentum is boosted to the lab; computeMoments must recover
// the temperatures k_B T_i = m c^2 <u'_i^2 / gamma'> of the rest-frame set.
TEST_F(DistributionMomentsTest, RestFrameTemperature_RecoversBoostedAnisotropicSet) {
    std::vector<Vec3> uRest = {
            Vec3(1.0e-2, 0.0, 0.0), Vec3(0.0, 2.0e-2, 0.0), Vec3(0.0, 0.0, 5.0e-3),
            Vec3(5.0e-3, -7.0e-3, 2.0e-3)};
    const size_t nHalf = uRest.size();
    for (size_t i = 0; i < nHalf; ++i) {
        uRest.push_back(Vec3(-uRest[i](0), -uRest[i](1), -uRest[i](2)));
    }
    const size_t N    = uRest.size();
    const size_t rank = ippl::Comm->rank(), ranks = ippl::Comm->size();
    const size_t nLocal = rank < N ? 1 + (N - 1 - rank) / ranks : 0;

    // The macro-particle mass in M differs from the species mass, which must be used.
    const double massGeV = Physics::m_p;
    Vec3 expected(0.0);
    for (const Vec3& u : uRest) {
        const double gammaRest = std::sqrt(1.0 + dot3(u, u));
        for (int d = 0; d < 3; ++d) {
            expected(d) += massGeV * Units::GeV2eV * u(d) * u(d) / gammaRest / Physics::kB / N;
        }
    }

    const Vec3 zAxis(0.0, 0.0, 1.0);
    const Vec3 oblique(1.0 / std::sqrt(14.0), -2.0 / std::sqrt(14.0), 3.0 / std::sqrt(14.0));
    for (const Vec3& direction : {zAxis, oblique}) {
        for (const double gammaBoost : {1.0, 10.0, 100.0}) {
            const double betaBoost = std::sqrt(1.0 - 1.0 / (gammaBoost * gammaBoost));
            const Vec3 beta(
                    betaBoost * direction(0), betaBoost * direction(1), betaBoost * direction(2));

            Kokkos::View<Vec3*> Rview_d("Rview_d", nLocal);
            Kokkos::View<Vec3*> Pview_d("Pview_d", nLocal);
            Kokkos::View<double*> Mview_d("Mview_d", 1);
            auto Rview_h = Kokkos::create_mirror_view(Rview_d);
            auto Pview_h = Kokkos::create_mirror_view(Pview_d);
            auto Mview_h = Kokkos::create_mirror_view(Mview_d);
            for (size_t i = 0; i < nLocal; ++i) {
                const size_t global = rank + i * ranks;
                Rview_h(i) =
                        Vec3(1.0e-3 * global, -2.0e-3 * global * global, 3.0e-3 * (N - global));
                Pview_h(i) = boostToLab(uRest[global], beta);
            }
            Mview_h(0) = 1.0e3 * massGeV;
            Kokkos::deep_copy(Rview_d, Rview_h);
            Kokkos::deep_copy(Pview_d, Pview_h);
            Kokkos::deep_copy(Mview_d, Mview_h);

            for (bool compact : {false, true}) {
                DistributionMoments dm;
                auto containerState = bunchStateHandler->registerContainer();
                dm.setContainerState(containerState);
                dm.setCompactSymmetricReduction(compact);
                dm.setEnergyReferenceMass(massGeV, true);
                dm.computeMoments(Rview_d, Pview_d, Mview_d, /*Np=*/N, nLocal);

                const Vec3 temperature = dm.getRestFrameTemperature();
                for (int d = 0; d < 3; ++d) {
                    EXPECT_NEAR(temperature(d), expected(d), 1e-8 * expected(d))
                            << "gamma = " << gammaBoost << ", direction = " << direction
                            << ", component " << d << ", compact = " << compact;
                }
            }
        }
    }
}
