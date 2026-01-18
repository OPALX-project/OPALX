#include <gtest/gtest.h>
#include <mpi.h>
#include <memory>
#include <cmath>

#include "Distribution/Gaussian.h"
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

TEST_F(GaussianTest, meanR_varR) {
    const Vector_t<double, 3> sigmaR = 0.5;
    const Vector_t<double, 3> sigmaP = 1.0;
    double avrgpz = 0.0;
    const Vector_t<double, 3> cutoffR = 4.0;

    Gaussian sampler(pc,
                    sigmaR, sigmaP, avrgpz, cutoffR);

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

    EXPECT_NEAR(meanR[0], 0.0, 1e-15);
    EXPECT_NEAR(meanR[1], 0.0, 1e-15);
    EXPECT_NEAR(meanR[2], 0.0, 1e-15);

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

    EXPECT_NEAR(stddevR[0], sigmaR[0], 1e-2);
    EXPECT_NEAR(stddevR[1], sigmaR[1], 1e-2);
    EXPECT_NEAR(stddevR[2], sigmaR[2], 1e-2);
}


TEST_F(GaussianTest, cutoffR)
{
    const Vector_t<double, 3> sigmaR = 1.0;
    const Vector_t<double, 3> sigmaP = 1.0;
    double avrgpz = 0.0;
    const Vector_t<double, 3> cutoffR = 3.0;

    bool fixMeanR = false;
    Gaussian sampler(pc, sigmaR, sigmaP, avrgpz, cutoffR, fixMeanR);

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
    EXPECT_LE(global_maxAbsR[0], cutoffR[0]);
    EXPECT_LE(global_maxAbsR[1], cutoffR[1]);
    EXPECT_LE(global_maxAbsR[2], cutoffR[2]);
}

TEST_F(GaussianTest, meanP_and_varP)
{
    const Vector_t<double, 3> sigmaR = 1.0;
    const Vector_t<double, 3> sigmaP = 1.0;
    double avrgpz = 0.1;
    const Vector_t<double, 3> cutoffR = 4.0;

    Gaussian sampler(pc,
                    sigmaR, sigmaP, avrgpz, cutoffR);

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
    EXPECT_NEAR(meanP[0], 0.0, 1e-2);
    EXPECT_NEAR(meanP[1], 0.0, 1e-2);
    EXPECT_NEAR(meanP[2], avrgpz, 1e-2);

    EXPECT_NEAR(varianceP[0], sigmaP[0]*sigmaP[0], 1e-2);
    EXPECT_NEAR(varianceP[1], sigmaP[1]*sigmaP[1], 1e-2);
    EXPECT_NEAR(varianceP[2], sigmaP[2]*sigmaP[2], 1e-2);
}
