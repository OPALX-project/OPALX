#include <gtest/gtest.h>
#include <mpi.h>
#include <memory>
#include <cmath>
#include <algorithm>

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
        nr = 32;
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

template <typename ViewType>
void computeMean(ViewType& view,
                  size_t nlocal,
                  size_t total_nparticles,
                  double meanR[3])
{
    double sumR[3] = {0.0, 0.0, 0.0};

    Kokkos::parallel_reduce(nlocal,
        KOKKOS_LAMBDA(const int k, double& cent0, double& cent1, double& cent2) {
                cent0 += view(k)[0];
                cent1 += view(k)[1];
                cent2 += view(k)[2];
        },
        Kokkos::Sum<double>(sumR[0]), Kokkos::Sum<double>(sumR[1]), Kokkos::Sum<double>(sumR[2])
    );
    Kokkos::fence();

    MPI_Allreduce(sumR, meanR, 3, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());
    ippl::Comm->barrier();

    for (int i=0; i<3; i++)
        meanR[i] /= total_nparticles;
}

template <typename ViewType>
void computeStdDev(ViewType &view,
                    size_t nlocal,
                    size_t total_nparticles,
                    const double meanR[3],
                    double stddevR[3])
{
    double sumR2[3] = {0.0, 0.0, 0.0};

    // Compute sum of squares locally
    Kokkos::parallel_reduce(nlocal,
        KOKKOS_LAMBDA(const int k, double& s0, double& s1, double& s2) {
            s0 += view(k)[0] * view(k)[0];
            s1 += view(k)[1] * view(k)[1];
            s2 += view(k)[2] * view(k)[2];
        },
        Kokkos::Sum<double>(sumR2[0]),
        Kokkos::Sum<double>(sumR2[1]),
        Kokkos::Sum<double>(sumR2[2])
    );
    Kokkos::fence();

    // MPI reduce to get global sums
    double global_sumR2[3] = {0.0, 0.0, 0.0};
    MPI_Allreduce(sumR2, global_sumR2, 3, MPI_DOUBLE, MPI_SUM,
                  ippl::Comm->getCommunicator());
    ippl::Comm->barrier();

    // Compute variance and stddev
    for (int i = 0; i < 3; i++) {
        double mean_square = global_sumR2[i] / static_cast<double>(total_nparticles);
        double variance = mean_square - meanR[i] * meanR[i];
        stddevR[i] = std::sqrt(variance);
    }
}

// Forward declaration so tests can call this helper before its definition.
static void preallocateParticleCapacity(
    const std::shared_ptr<ParticleContainer<double, 3>>& pc,
    size_t totalParticles);

TEST_F(MultiVariateGaussianTest, meanR_varR) {
    const Vector_t<double, 3> sigmaR_ref = 0.5;
    const Vector_t<double, 3> sigmaP_ref = 1.0;
    const Vector_t<double, 3> meanR_ref = 0.0;
    const Vector_t<double, 3> meanP_ref = 0.0;
    const Vector_t<double, 3> cutoffR_ref = 4.0;
    const Vector_t<double, 3> cutoffP_ref = 4.0;

    MultiVariateGaussian sampler(pc, meanR_ref, meanP_ref, sigmaR_ref, sigmaP_ref, cutoffR_ref, cutoffP_ref);

    size_t total_nparticles = 100000;

    preallocateParticleCapacity(pc, total_nparticles);
    sampler.generateParticles(total_nparticles, nr);

    double meanR[3];
    computeMean(pc->R.getView(), pc->getLocalNum(), total_nparticles, meanR);

    EXPECT_NEAR(meanR[0], meanR_ref[0], 1e-15);
    EXPECT_NEAR(meanR[1], meanR_ref[1], 1e-15);
    EXPECT_NEAR(meanR[2], meanR_ref[2], 1e-15);

    // compute stddev of R
    double stddevR[3];
    computeStdDev(pc->R.getView(), pc->getLocalNum(), total_nparticles, meanR, stddevR);

    EXPECT_NEAR(stddevR[0], sigmaR_ref[0], 1e-2);
    EXPECT_NEAR(stddevR[1], sigmaR_ref[1], 1e-2);
    EXPECT_NEAR(stddevR[2], sigmaR_ref[2], 1e-2);
}

template <typename ViewType>
void computeMaxAbsR(ViewType &view,
                    size_t nlocal,
                    double global_maxAbsR[3],
                    double& global_maxRadius)
{
    // Local maxima initialized
    double local_maxAbsR[3] = {0.0, 0.0, 0.0};
    double local_maxRadius = 0.0;

    // Compute local maxima using Kokkos
    Kokkos::parallel_reduce(
        nlocal,
        KOKKOS_LAMBDA(const int k,
                    double& max0, double& max1, double& max2, double& maxr) {
            double x = view(k)[0];
            double y = view(k)[1];
            double z = view(k)[2];

            double ax = Kokkos::fabs(x);
            double ay = Kokkos::fabs(y);
            double az = Kokkos::fabs(z);
            double r  = Kokkos::sqrt(x*x + y*y + z*z);

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


    // MPI reduction to get global maxima
    MPI_Allreduce(local_maxAbsR, global_maxAbsR, 3, MPI_DOUBLE, MPI_MAX,
                  ippl::Comm->getCommunicator());
    MPI_Allreduce(&local_maxRadius, &global_maxRadius, 1, MPI_DOUBLE, MPI_MAX,
                  ippl::Comm->getCommunicator());

    ippl::Comm->barrier();
}

/// Preallocate particle container capacity so SamplingBase::computeLocalEmitCount
/// can distribute particles without being constrained by zero capacity.
static void preallocateParticleCapacity(
    const std::shared_ptr<ParticleContainer<double, 3>>& pc,
    size_t totalParticles)
{
    const int nranks     = std::max(1, ippl::Comm->size());
    const size_t nranksU = static_cast<size_t>(nranks);

    const size_t maxLocalNum =
        totalParticles / nranksU + 2 * nranksU + 1;

    pc->create(maxLocalNum);
    Kokkos::View<bool*> tmp_invalid("tmp_invalid", maxLocalNum);
    pc->destroy(tmp_invalid, maxLocalNum);
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

    size_t total_nparticles = 100000;
    preallocateParticleCapacity(pc, total_nparticles);
    sampler.generateParticles(total_nparticles, nr);
    
    double global_maxAbsR[3];
    double global_maxRadius;
    computeMaxAbsR(pc->R.getView(), pc->getLocalNum(), global_maxAbsR, global_maxRadius);

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

    size_t total_nparticles = 100000;
    preallocateParticleCapacity(pc, total_nparticles);
    sampler.generateParticles(total_nparticles, nr);

    double meanP[3] = {0.0, 0.0, 0.0};
    computeMean(pc->P.getView(), pc->getLocalNum(), total_nparticles, meanP);

    double sigmaP[3] = {0.0, 0.0, 0.0};
    computeStdDev(pc->P.getView(), pc->getLocalNum(), total_nparticles, meanP,  sigmaP);

    EXPECT_NEAR(meanP[0], meanP_ref[0], 1e-2);
    EXPECT_NEAR(meanP[1], meanP_ref[1], 1e-2);
    EXPECT_NEAR(meanP[2], meanP_ref[2], 1e-2);

    EXPECT_NEAR(sigmaP[0], sigmaP_ref[0], 1e-2);
    EXPECT_NEAR(sigmaP[1], sigmaP_ref[1], 1e-2);
    EXPECT_NEAR(sigmaP[2], sigmaP_ref[2], 1e-2);
}

void computeMeanRAndP(const std::shared_ptr<ParticleContainer<double,3>>& pc,
                      size_t total_nparticles,
                      double meanSample[6])
{
    auto Rview = pc->R.getView();
    auto Pview = pc->P.getView();

    double loc_sum[6] = {0.0};

    // Kokkos parallel_reduce for both R and P

    Kokkos::parallel_reduce(
        pc->getLocalNum(),
        KOKKOS_LAMBDA(const int k,
                    double& s0, double& s1, double& s2,
                    double& s3, double& s4, double& s5) {
            s0 += Rview(k)[0];
            s1 += Pview(k)[0];
            s2 += Rview(k)[1];
            s3 += Pview(k)[1];
            s4 += Rview(k)[2];
            s5 += Pview(k)[2];
        },
        Kokkos::Sum<double>(loc_sum[0]),
        Kokkos::Sum<double>(loc_sum[1]),
        Kokkos::Sum<double>(loc_sum[2]),
        Kokkos::Sum<double>(loc_sum[3]),
        Kokkos::Sum<double>(loc_sum[4]),
        Kokkos::Sum<double>(loc_sum[5])
    );
    Kokkos::fence();

    // MPI reduction to global sum
    double global_sum[6] = {0.0};
    MPI_Allreduce(loc_sum, global_sum, 6, MPI_DOUBLE, MPI_SUM,
                  ippl::Comm->getCommunicator());
    ippl::Comm->barrier();

    // Normalize to get mean
    for (int i = 0; i < 6; i++)
        meanSample[i] = global_sum[i] / static_cast<double>(total_nparticles);
}

void computeCovariance6x6(const std::shared_ptr<ParticleContainer<double,3>>& pc,
                          const double meanSample[6],
                          size_t total_nparticles,
                          double moment[6][6])
{    
    auto Rview = pc->R.getView();
    auto Pview = pc->P.getView();

    Kokkos::Array<double, 6> mean = {
        meanSample[0],
        meanSample[1],
        meanSample[2],
        meanSample[3],
        meanSample[4],
        meanSample[5]
    };

    double loc_moment[6][6] = {0.0};

    for (unsigned i = 0; i < 6; ++i) {
        const unsigned row = i;
        Kokkos::parallel_reduce(
            pc->getLocalNum(),
            KOKKOS_LAMBDA(const int k,
                        double& s0, double& s1, double& s2,
                        double& s3, double& s4, double& s5) {
                double part[6] = {
                    Rview(k)[0] - mean[0],
                    Pview(k)[0] - mean[1],
                    Rview(k)[1] - mean[2],
                    Pview(k)[1] - mean[3],
                    Rview(k)[2] - mean[4],
                    Pview(k)[2] - mean[5]
                };
                s0 += part[row] * part[0];
                s1 += part[row] * part[1];
                s2 += part[row] * part[2];
                s3 += part[row] * part[3];
                s4 += part[row] * part[4];
                s5 += part[row] * part[5];
            },
            Kokkos::Sum<double>(loc_moment[row][0]),
            Kokkos::Sum<double>(loc_moment[row][1]),
            Kokkos::Sum<double>(loc_moment[row][2]),
            Kokkos::Sum<double>(loc_moment[row][3]),
            Kokkos::Sum<double>(loc_moment[row][4]),
            Kokkos::Sum<double>(loc_moment[row][5])
        );
        Kokkos::fence();
    }

    // MPI reduction across all ranks
    MPI_Allreduce(loc_moment, moment, 6 * 6, MPI_DOUBLE, MPI_SUM,
                  ippl::Comm->getCommunicator());
    ippl::Comm->barrier();

    // Normalize to get average covariance/moments
    for (unsigned i = 0; i < 6; ++i) {
        for (unsigned j = 0; j < 6; ++j) {
            moment[i][j] /= static_cast<double>(total_nparticles);
        }
    }
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
    using matrix_t = ippl::Vector<ippl::Vector<double, 6>, 6>;
    matrix_t L;
    matrix_t covExpected;
    
    for (int i=0;i<6;i++){
        for (int j=0;j<6;j++){
            L[i][j] = 0.0;
            covExpected[i][j] = 0.0;
        }
    }


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
    preallocateParticleCapacity(pc, total_nparticles);
    sampler.generateParticles(total_nparticles, nr);

    double meanSample[6];
    computeMeanRAndP(pc, total_nparticles, meanSample);

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
    double moment[6][6];
    computeCovariance6x6(pc, meanSample, total_nparticles, moment);

    // -----------------------------------------------------
    // Assertions: covariance comparison
    // -----------------------------------------------------
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){
            EXPECT_NEAR(moment[i][j], covExpected[i][j], 3e-2);
        }
    }
}
