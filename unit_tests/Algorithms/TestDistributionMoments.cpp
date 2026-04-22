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
        e.stdEkinMeV_likeImplementation = std::sqrt(e.stdEkinMeV_likeImplementation / static_cast<double>(N));

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

    static void TearDownTestSuite() {
        ippl::finalize();
    }

    void SetUp() override {
        bunchStateHandler = std::make_shared<BunchStateHandler>();
    }

    std::shared_ptr<BunchStateHandler> bunchStateHandler;

    void TearDown() override {
        bunchStateHandler.reset();
    }
};

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
