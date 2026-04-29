#include <gtest/gtest.h>
#include <mpi.h>
#include <memory>
#include <cmath>
#include <algorithm>

#include "Distribution/Gaussian.h"
#include "PartBunch/BunchStateHandler.h"
#include "Ippl.h"
#include "Utility/IpplTimings.h"

class GaussianTest : public ::testing::Test {
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

        bunchStateHandler = std::make_shared<BunchStateHandler>();
        pc->setBunchStateHandler(bunchStateHandler);
    }

    void TearDown() override {
        // nothing special
    }

    std::shared_ptr<ParticleContainer<double, 3>> pc;
    std::shared_ptr<BunchStateHandler> bunchStateHandler;
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

    // Mirror TrackRun: give each rank enough headroom above its ideal share.
    const size_t maxLocalNum =
        totalParticles / nranksU + 2 * nranksU + 1;

    pc->allocateParticles(maxLocalNum);
}

TEST_F(GaussianTest, meanR_stddevR) {
    const Vector_t<double, 3> sigmaR_ref = 0.5;
    const Vector_t<double, 3> sigmaP_ref = 1.0;
    double avrgpz = 0.0;
    const Vector_t<double, 3> cutoffR = 4.0;

    Gaussian sampler(pc, sigmaR_ref, sigmaP_ref, avrgpz, cutoffR);

    size_t total_nparticles = 100000;

    preallocateParticleCapacity(pc, total_nparticles);
    sampler.generateParticles(total_nparticles, nr);

    auto Rview = pc->R.getView();
    auto Pview = pc->P.getView();

    // compute mean of R
    double meanR[3] = {0.0, 0.0, 0.0};
    computeMean(Rview, pc->getLocalNum(), total_nparticles, meanR);

    EXPECT_NEAR(meanR[0], 0.0, 1e-15);
    EXPECT_NEAR(meanR[1], 0.0, 1e-15);
    EXPECT_NEAR(meanR[2], 0.0, 1e-15);

    // compute variance of R
    double sigmaR[3] = {0.0, 0.0, 0.0};
    computeStdDev(Rview, pc->getLocalNum(), total_nparticles, meanR, sigmaR);

    EXPECT_NEAR(sigmaR[0], sigmaR_ref[0], 1e-2);
    EXPECT_NEAR(sigmaR[1], sigmaR_ref[1], 1e-2);
    EXPECT_NEAR(sigmaR[2], sigmaR_ref[2], 1e-2);
}


TEST_F(GaussianTest, cutoffR)
{
    const Vector_t<double, 3> sigmaR = 1.0;
    const Vector_t<double, 3> sigmaP = 1.0;
    double avrgpz = 0.0;
    const Vector_t<double, 3> cutoffR = 3.0;

    bool fixMeanR = false;
    Gaussian sampler(pc, sigmaR, sigmaP, avrgpz, cutoffR, fixMeanR);

    size_t total_nparticles = 100000;

    preallocateParticleCapacity(pc, total_nparticles);
    sampler.generateParticles(total_nparticles, nr);

    auto Rview = pc->R.getView();

    // compute max of abs(R)
    double global_maxAbsR[3] = {0.0, 0.0, 0.0};
    double global_maxRadius = 0.0;
    computeMaxAbsR(Rview, pc->getLocalNum(), global_maxAbsR, global_maxRadius);

    // --------- Assertions ---------

    // Each coordinate individually should be inside cutoff
    EXPECT_LE(global_maxAbsR[0], cutoffR[0]);
    EXPECT_LE(global_maxAbsR[1], cutoffR[1]);
    EXPECT_LE(global_maxAbsR[2], cutoffR[2]);
}

TEST_F(GaussianTest, meanP_and_steddevP)
{
    const Vector_t<double, 3> sigmaR_ref = 1.0;
    const Vector_t<double, 3> sigmaP_ref = 1.0;
    double avrgpz = 0.1;
    const Vector_t<double, 3> cutoffR = 4.0;

    Gaussian sampler(pc, sigmaR_ref, sigmaP_ref, avrgpz, cutoffR);

    size_t total_nparticles = 100000;
    preallocateParticleCapacity(pc, total_nparticles);
    sampler.generateParticles(total_nparticles, nr);

    auto Pview = pc->P.getView();

    //---------------------------------------------------------------------
    // 1. Compute mean(P)
    //---------------------------------------------------------------------
    double meanP[3] = {0.0, 0.0, 0.0};

    computeMean(Pview, pc->getLocalNum(), total_nparticles, meanP);

    //---------------------------------------------------------------------
    // 2. Compute stddev(P)
    //---------------------------------------------------------------------
    double sigmaP[3] = {0.0, 0.0, 0.0};
    computeStdDev(Pview, pc->getLocalNum(), total_nparticles, meanP, sigmaP);

    //---------------------------------------------------------------------
    // 3. Expectations
    //---------------------------------------------------------------------
    EXPECT_NEAR(meanP[0], 0.0, 1e-2);
    EXPECT_NEAR(meanP[1], 0.0, 1e-2);
    EXPECT_NEAR(meanP[2], avrgpz, 1e-2);

    EXPECT_NEAR(sigmaP[0], sigmaP_ref[0], 1e-2);
    EXPECT_NEAR(sigmaP[1], sigmaP_ref[1], 1e-2);
    EXPECT_NEAR(sigmaP[2], sigmaP_ref[2], 1e-2);
}
