#include <gtest/gtest.h>
#include <mpi.h>
#include <memory>
#include <cmath>

#include "Distribution/MultiVariateGaussian.h"
#include "Ippl.h"
#include "Utility/IpplTimings.h"


class MultiVariateGaussianTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc = 0;
        char** argv = nullptr;

        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() {
        ippl::finalize();
    }

    void SetUp() override {
        // Minimal 3D grid parameters
        nr = 64;
        ippl::Vector<double,3> rmin = -4.0;
        ippl::Vector<double,3> rmax = 4.0;
        ippl::Vector<double,3> origin = rmin;
        ippl::Vector<double,3> hr = (rmax - rmin) / ippl::Vector<double,3>(nr);
        std::array<bool,3> decomp = {true, true, true};

        ippl::NDIndex<3> domain;
        for (unsigned i = 0; i < 3; i++) {
            domain[i] = ippl::Index(this->nr[i]);
        }

        // Create FieldContainer
        auto fc = std::make_shared<FieldContainer_t>(hr, rmin, rmax, decomp, domain, origin, this->isAllPeriodic_m);

        // Create mesh and fieldlayout
        Mesh_t<3> mesh(domain, hr, origin);
        FieldLayout_t<3> fl(MPI_COMM_WORLD, domain, decomp, this->isAllPeriodic_m);

        // Create ParticleContainer
        pc = std::make_shared<ParticleContainer<double, 3>>(mesh, fl);
    }

    void TearDown() override {
        // nothing special
    }

    std::shared_ptr<ParticleContainer<double, 3>> pc;
    ippl::Vector<int,3> nr;
    bool isAllPeriodic_m = true;
};

TEST_F(MultiVariateGaussianTest, meanR_varR) {
    const Vector_t<double, 3> sigmaR_ref = 0.5;
    const Vector_t<double, 3> sigmaP_ref = 1.0;
    const Vector_t<double, 3> meanR_ref = 0.0;
    const Vector_t<double, 3> meanP_ref = 0.0;
    const Vector_t<double, 3> cutoffR_ref = 4.0;
    const Vector_t<double, 3> cutoffP_ref = 4.0;

    MultiVariateGaussian sampler(pc, meanR_ref, meanP_ref, sigmaR_ref, sigmaP_ref, cutoffR_ref, cutoffP_ref);

    size_t total_nparticles = 1000000;

    sampler.generateParticles(total_nparticles, nr);

    auto Rview = pc->R.getView();
    auto Pview = pc->P.getView();

    // compute mean of R
    double meanR[3], sumR[3];
    for(int i=0; i<3; i++){
        meanR[i] = 0.0;
        sumR[i] = 0.0;
    }

    Kokkos::parallel_reduce(
            "calc moments of particle distr.", pc->getLocalNum(),
            KOKKOS_LAMBDA(
                    const int k, double& cent0, double& cent1, double& cent2) {
                    cent0 += Rview(k)[0];
                    cent1 += Rview(k)[1];
                    cent2 += Rview(k)[2];
            },
            Kokkos::Sum<double>(sumR[0]), Kokkos::Sum<double>(sumR[1]), Kokkos::Sum<double>(sumR[2]));
    Kokkos::fence();

    MPI_Allreduce(sumR, meanR, 3, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());
    ippl::Comm->barrier();

    for(int i=0; i<3; i++){
       meanR[i] = meanR[i]/(1.0*total_nparticles);
    }

    EXPECT_NEAR(meanR[0], meanR_ref[0], 1e-15);
    EXPECT_NEAR(meanR[1], meanR_ref[1], 1e-15);
    EXPECT_NEAR(meanR[2], meanR_ref[2], 1e-15);

    // compute stddev of R
    double sumR2[3] = {0.0, 0.0, 0.0};
    double global_sumR2[3] = {0.0, 0.0, 0.0};

    Kokkos::parallel_reduce(
        "calc second moments of particle distr.", pc->getLocalNum(),
        KOKKOS_LAMBDA(const int k, double& c0, double& c1, double& c2) {
            c0 += Rview(k)[0] * Rview(k)[0];
            c1 += Rview(k)[1] * Rview(k)[1];
            c2 += Rview(k)[2] * Rview(k)[2];
        },
        Kokkos::Sum<double>(sumR2[0]),
        Kokkos::Sum<double>(sumR2[1]),
        Kokkos::Sum<double>(sumR2[2])
    );
    Kokkos::fence();

    // MPI reduce second moment
    MPI_Allreduce(sumR2, global_sumR2, 3, MPI_DOUBLE, MPI_SUM,
                ippl::Comm->getCommunicator());
    ippl::Comm->barrier();

    // compute variance = E[R^2] - (E[R])^2
    double varianceR[3];
    for (int i = 0; i < 3; i++) {
        double mean_square = global_sumR2[i] / (1.0 * total_nparticles);
        varianceR[i] = mean_square - meanR[i] * meanR[i];
    }
    double stddevR[3];
    for (int i = 0; i < 3; i++) {
        stddevR[i] = std::sqrt(varianceR[i]);
    }

    EXPECT_NEAR(stddevR[0], sigmaR_ref[0], 1e-2);
    EXPECT_NEAR(stddevR[1], sigmaR_ref[1], 1e-2);
    EXPECT_NEAR(stddevR[2], sigmaR_ref[2], 1e-2);
}


TEST_F(MultiVariateGaussianTest, cutoffR)
{
    const Vector_t<double, 3> sigmaR_ref = 0.5;
    const Vector_t<double, 3> sigmaP_ref = 1.0;
    const Vector_t<double, 3> meanR_ref = 0.0;
    const Vector_t<double, 3> meanP_ref = 0.0;
    const Vector_t<double, 3> cutoffR_ref = 4.0;
    const Vector_t<double, 3> cutoffP_ref = 4.0;

    bool fixMeanR = false;
    bool fixMeanP = false;
    MultiVariateGaussian sampler(pc, meanR_ref, meanP_ref, sigmaR_ref, sigmaP_ref, cutoffR_ref, cutoffP_ref, fixMeanR, fixMeanP);

    size_t total_nparticles = 1000000;

    sampler.generateParticles(total_nparticles, nr);

    auto Rview = pc->R.getView();

    // Track maximum absolute value in each dimension
    double local_maxAbsR[3] = {0.0, 0.0, 0.0};
    double global_maxAbsR[3] = {0.0, 0.0, 0.0};

    // Also track maximum radial magnitude if needed:
    double local_maxRadius = 0.0;
    double global_maxRadius = 0.0;

    // Reduce max along each dimension + radius
    Kokkos::parallel_reduce(
        "check cutoff R",
        pc->getLocalNum(),
        KOKKOS_LAMBDA(const int k, double& max0, double& max1, double& max2, double& maxr)
        {
            double x = Rview(k)[0];
            double y = Rview(k)[1];
            double z = Rview(k)[2];

            double ax = fabs(x);
            double ay = fabs(y);
            double az = fabs(z);
            double r  = sqrt(x*x + y*y + z*z);

            if (ax > max0) max0 = ax;
            if (ay > max1) max1 = ay;
            if (az > max2) max2 = az;
            if (r  > maxr) maxr  = r;
        },
        Kokkos::Max<double>(local_maxAbsR[0]),
        Kokkos::Max<double>(local_maxAbsR[1]),
        Kokkos::Max<double>(local_maxAbsR[2]),
        Kokkos::Max<double>(local_maxRadius)
    );
    Kokkos::fence();

    // MPI reduction across ranks
    MPI_Allreduce(local_maxAbsR, global_maxAbsR, 3,
                  MPI_DOUBLE, MPI_MAX,
                  ippl::Comm->getCommunicator());

    MPI_Allreduce(&local_maxRadius, &global_maxRadius, 1,
                  MPI_DOUBLE, MPI_MAX,
                  ippl::Comm->getCommunicator());

    ippl::Comm->barrier();

    // --------- Assertions ---------

    // Each coordinate individually should be inside cutoff
    EXPECT_LE(global_maxAbsR[0], cutoffR_ref[0]);
    EXPECT_LE(global_maxAbsR[1], cutoffR_ref[1]);
    EXPECT_LE(global_maxAbsR[2], cutoffR_ref[2]);
}

TEST_F(MultiVariateGaussianTest, meanP_and_varP)
{
    const Vector_t<double, 3> sigmaR_ref = 0.5;
    const Vector_t<double, 3> sigmaP_ref = 1.0;
    const Vector_t<double, 3> meanR_ref = 0.0;
    const Vector_t<double, 3> meanP_ref = 0.0;
    const Vector_t<double, 3> cutoffR_ref = 4.0;
    const Vector_t<double, 3> cutoffP_ref = 4.0;
    bool fixMeanR = true;
    bool fixMeanP = true;

    MultiVariateGaussian sampler(pc, meanR_ref, meanP_ref, sigmaR_ref, sigmaP_ref, cutoffR_ref, cutoffP_ref, fixMeanR, fixMeanP);

    size_t total_nparticles = 1000000;
    sampler.generateParticles(total_nparticles, nr);

    auto Pview = pc->P.getView();

    //---------------------------------------------------------------------
    // 1. Compute mean(P)
    //---------------------------------------------------------------------
    double sumP[3]      = {0.0, 0.0, 0.0};
    double global_sumP[3] = {0.0, 0.0, 0.0};

    Kokkos::parallel_reduce(
        "compute mean P", pc->getLocalNum(),
        KOKKOS_LAMBDA(const int k, double& s0, double& s1, double& s2)
        {
            s0 += Pview(k)[0];
            s1 += Pview(k)[1];
            s2 += Pview(k)[2];
        },
        Kokkos::Sum<double>(sumP[0]),
        Kokkos::Sum<double>(sumP[1]),
        Kokkos::Sum<double>(sumP[2])
    );
    Kokkos::fence();

    MPI_Allreduce(sumP, global_sumP, 3,
                  MPI_DOUBLE, MPI_SUM,
                  ippl::Comm->getCommunicator());
    ippl::Comm->barrier();

    double meanP[3];
    for (int i = 0; i < 3; i++)
        meanP[i] = global_sumP[i] / (1.0 * total_nparticles);


    //---------------------------------------------------------------------
    // 2. Compute variance(P)
    //---------------------------------------------------------------------
    double sumP2[3] = {0.0, 0.0, 0.0};
    double global_sumP2[3] = {0.0, 0.0, 0.0};

    Kokkos::parallel_reduce(
        "compute variance P", pc->getLocalNum(),
        KOKKOS_LAMBDA(const int k, double& s0, double& s1, double& s2)
        {
            double px = Pview(k)[0];
            double py = Pview(k)[1];
            double pz = Pview(k)[2];

            s0 += px * px;
            s1 += py * py;
            s2 += pz * pz;
        },
        Kokkos::Sum<double>(sumP2[0]),
        Kokkos::Sum<double>(sumP2[1]),
        Kokkos::Sum<double>(sumP2[2])
    );
    Kokkos::fence();

    MPI_Allreduce(sumP2, global_sumP2, 3,
                  MPI_DOUBLE, MPI_SUM,
                  ippl::Comm->getCommunicator());
    ippl::Comm->barrier();

    double varianceP[3];
    for (int i = 0; i < 3; i++)
    {
        double E_P2 = global_sumP2[i] / (1.0 * total_nparticles);
        varianceP[i] = E_P2 - meanP[i] * meanP[i];
    }

    //---------------------------------------------------------------------
    // 3. Expectations
    //---------------------------------------------------------------------
    EXPECT_NEAR(meanP[0], meanP_ref[0], 1e-2);
    EXPECT_NEAR(meanP[1], meanP_ref[1], 1e-2);
    EXPECT_NEAR(meanP[2], meanP_ref[2], 1e-2);

    EXPECT_NEAR(varianceP[0], sigmaP_ref[0]*sigmaP_ref[0], 1e-2);
    EXPECT_NEAR(varianceP[1], sigmaP_ref[1]*sigmaP_ref[1], 1e-2);
    EXPECT_NEAR(varianceP[2], sigmaP_ref[2]*sigmaP_ref[2], 1e-2);
}

TEST_F(MultiVariateGaussianTest, FullCovarianceTest)
{
    // Means
    Vector_t<double, 3> meanR = {0.2, -0.1, 0.4};
    Vector_t<double, 3> meanP = {-0.3, 0.5, 0.1};

    // Std dev values
    Vector_t<double, 3> sigmaR = {0.9, 2.0, 1.0};
    Vector_t<double, 3> sigmaP = {1.5, 1.0, 0.9};

    // Build L (lower-triangular)
    Matrix_t L = {0.0};

    // Diagonal entries: sqrt of variances
    for (int i = 0; i < 3; i++) {
        L[i*2][i*2]     = sigmaR[i];
        L[i*2+1][i*2+1] = sigmaP[i];
    }

    // Fill below-diagonal elements with some pattern
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < i; j++) {
            L[i][j] = (i + j + 1) * 0.1;  // arbitrary correlation
        }
    }

    // Compute covariance matrix Σ = L * L^T
    Matrix_t covExpected = {0.0};

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            double sum = 0.0;
            for (int k = 0; k < 6; k++)
                sum += L[i][k] * L[j][k];  // (L * L^T)
            covExpected[i][j] = sum;
        }
    }

    Vector_t<double, 3> cutoffR = 5.0;
    Vector_t<double, 3> cutoffP = 5.0;

    // Construct sampler
    MultiVariateGaussian sampler(pc, meanR, meanP, covExpected, cutoffR, cutoffP);

    size_t total_nparticles = 1000000;
    sampler.generateParticles(total_nparticles, nr);

    auto Rview = pc->R.getView();
    auto Pview = pc->P.getView();

    // ---------------------------
    // Compute sample means
    // ---------------------------
    double loc_sum[6] = {0.0};
    double global_sum[6] = {0.0};

    Kokkos::parallel_reduce(
        "mean", pc->getLocalNum(),
        KOKKOS_LAMBDA(const int k, double& mom0, double& mom1, double& mom2,
                    double& mom3, double& mom4, double& mom5) {
            mom0 += Rview(k)[0];
            mom1 += Pview(k)[0];
            mom2 += Rview(k)[1];
            mom3 += Pview(k)[1];
            mom4 += Rview(k)[2];
            mom5 += Pview(k)[2];
        },
        Kokkos::Sum<double>(loc_sum[0]),
        Kokkos::Sum<double>(loc_sum[1]),
        Kokkos::Sum<double>(loc_sum[2]),
        Kokkos::Sum<double>(loc_sum[3]),
        Kokkos::Sum<double>(loc_sum[4]),
        Kokkos::Sum<double>(loc_sum[5])
    );
    Kokkos::fence();

    MPI_Allreduce(loc_sum, global_sum, 6, MPI_DOUBLE, MPI_SUM,
                  ippl::Comm->getCommunicator());

    double meanSample[6];
    for (int i = 0; i < 6; i++)
        meanSample[i] = global_sum[i] / total_nparticles;

    // Check means
    EXPECT_NEAR(meanSample[0], meanR[0], 5e-3);
    EXPECT_NEAR(meanSample[1], meanP[0], 5e-3);
    EXPECT_NEAR(meanSample[2], meanR[1], 5e-3);
    EXPECT_NEAR(meanSample[3], meanP[1], 5e-3);
    EXPECT_NEAR(meanSample[4], meanR[2], 5e-3);
    EXPECT_NEAR(meanSample[5], meanP[2], 5e-3);

    // ------------------------------
    // Compute sample covariance Σ
    // ------------------------------
    double loc_moment[6][6] = {};
    double moment[6][6]  = {};

    for (unsigned i = 0; i < 6; ++i) {
            Kokkos::parallel_reduce("calc moments of particle distr.", pc->getLocalNum(),
                KOKKOS_LAMBDA(
                    const int k, double& mom0, double& mom1, double& mom2,
                    double& mom3, double& mom4, double& mom5) {
                    double part[6];
                    part[0] = Rview(k)[0]-meanSample[0];
                    part[1] = Pview(k)[0]-meanSample[1];
                    part[2] = Rview(k)[1]-meanSample[2];
                    part[3] = Pview(k)[1]-meanSample[3];
                    part[4] = Rview(k)[2]-meanSample[4];
                    part[5] = Pview(k)[2]-meanSample[5];

                    mom0 += part[i] * part[0];
                    mom1 += part[i] * part[1];
                    mom2 += part[i] * part[2];
                    mom3 += part[i] * part[3];
                    mom4 += part[i] * part[4];
                    mom5 += part[i] * part[5];
                },
                Kokkos::Sum<double>(loc_moment[i][0]),
                Kokkos::Sum<double>(loc_moment[i][1]), Kokkos::Sum<double>(loc_moment[i][2]),
                Kokkos::Sum<double>(loc_moment[i][3]), Kokkos::Sum<double>(loc_moment[i][4]),
                Kokkos::Sum<double>(loc_moment[i][5]));
            Kokkos::fence();
     }

    MPI_Allreduce(
            loc_moment, moment, 6 * 6, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());

    for (unsigned i = 0; i < 6; i++) {
            for (unsigned j = 0; j < 6; j++) {
                moment[i][j] = moment[i][j] / total_nparticles;
            }
     }

    // -----------------------------------------------------
    // Assertions: covariance comparison
    // -----------------------------------------------------
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){
            EXPECT_NEAR(moment[i][j], covExpected[i][j], 3e-2);
        }
    }
}
