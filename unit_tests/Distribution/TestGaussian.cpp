#include <gtest/gtest.h>
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

#include "Attributes/Attributes.h"
#include "Distribution/Gaussian.h"
#include "Ippl.h"
#include "PartBunch/BunchStateHandler.h"
#include "Utility/IpplTimings.h"

class GaussianTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;

        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() { ippl::finalize(); }

    void SetUp() override { resetContainer(); }

    // Each collapsed-axis case needs an empty container without reinitializing MPI/Kokkos.
    void resetContainer() {
        // Minimal 3D grid parameters
        nr                             = 32;
        ippl::Vector<double, 3> rmin   = -4.0;
        ippl::Vector<double, 3> rmax   = 4.0;
        ippl::Vector<double, 3> origin = rmin;
        ippl::Vector<double, 3> hr     = (rmax - rmin) / ippl::Vector<double, 3>(nr);
        std::array<bool, 3> decomp     = {true, true, true};

        ippl::NDIndex<3> domain;
        for (unsigned i = 0; i < 3; i++) {
            domain[i] = ippl::Index(this->nr[i]);
        }

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
    ippl::Vector<int, 3> nr;
    bool isAllPeriodic_m = true;
};

template <typename ViewType>
void computeMean(ViewType& view, size_t nlocal, size_t total_nparticles, double meanR[3]) {
    double sumR[3] = {0.0, 0.0, 0.0};

    Kokkos::parallel_reduce(
            nlocal,
            KOKKOS_LAMBDA(const int k, double& cent0, double& cent1, double& cent2) {
                cent0 += view(k)[0];
                cent1 += view(k)[1];
                cent2 += view(k)[2];
            },
            Kokkos::Sum<double>(sumR[0]), Kokkos::Sum<double>(sumR[1]),
            Kokkos::Sum<double>(sumR[2]));
    Kokkos::fence();

    MPI_Allreduce(sumR, meanR, 3, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());
    ippl::Comm->barrier();

    for (int i = 0; i < 3; i++)
        meanR[i] /= total_nparticles;
}

template <typename ViewType>
void computeStdDev(
        ViewType& view, size_t nlocal, size_t total_nparticles, const double meanR[3],
        double stddevR[3]) {
    double sumR2[3] = {0.0, 0.0, 0.0};

    // Compute sum of squares locally
    Kokkos::parallel_reduce(
            nlocal,
            KOKKOS_LAMBDA(const int k, double& s0, double& s1, double& s2) {
                s0 += view(k)[0] * view(k)[0];
                s1 += view(k)[1] * view(k)[1];
                s2 += view(k)[2] * view(k)[2];
            },
            Kokkos::Sum<double>(sumR2[0]), Kokkos::Sum<double>(sumR2[1]),
            Kokkos::Sum<double>(sumR2[2]));
    Kokkos::fence();

    // MPI reduce to get global sums
    double global_sumR2[3] = {0.0, 0.0, 0.0};
    MPI_Allreduce(sumR2, global_sumR2, 3, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());
    ippl::Comm->barrier();

    // Compute variance and stddev
    for (int i = 0; i < 3; i++) {
        double mean_square = global_sumR2[i] / static_cast<double>(total_nparticles);
        double variance    = mean_square - meanR[i] * meanR[i];
        stddevR[i]         = std::sqrt(variance);
    }
}

template <typename ViewType>
void computeMaxAbsR(
        ViewType& view, size_t nlocal, double global_maxAbsR[3], double& global_maxRadius) {
    // Local maxima initialized
    double local_maxAbsR[3] = {0.0, 0.0, 0.0};
    double local_maxRadius  = 0.0;

    // Compute local maxima using Kokkos
    Kokkos::parallel_reduce(
            nlocal,
            KOKKOS_LAMBDA(const int k, double& max0, double& max1, double& max2, double& maxr) {
                double x = view(k)[0];
                double y = view(k)[1];
                double z = view(k)[2];

                double ax = Kokkos::fabs(x);
                double ay = Kokkos::fabs(y);
                double az = Kokkos::fabs(z);
                double r  = Kokkos::sqrt(x * x + y * y + z * z);

                if (ax > max0) max0 = ax;
                if (ay > max1) max1 = ay;
                if (az > max2) max2 = az;
                if (r > maxr) maxr = r;
            },
            Kokkos::Max<double>(local_maxAbsR[0]), Kokkos::Max<double>(local_maxAbsR[1]),
            Kokkos::Max<double>(local_maxAbsR[2]), Kokkos::Max<double>(local_maxRadius));
    Kokkos::fence();

    // MPI reduction to get global maxima
    MPI_Allreduce(
            local_maxAbsR, global_maxAbsR, 3, MPI_DOUBLE, MPI_MAX, ippl::Comm->getCommunicator());
    MPI_Allreduce(
            &local_maxRadius, &global_maxRadius, 1, MPI_DOUBLE, MPI_MAX,
            ippl::Comm->getCommunicator());

    ippl::Comm->barrier();
}

/// Preallocate particle container capacity so SamplingBase::computeLocalEmitCount
/// can distribute particles without being constrained by zero capacity.
static void preallocateParticleCapacity(
        const std::shared_ptr<ParticleContainer<double, 3>>& pc, size_t totalParticles) {
    const int nranks     = std::max(1, ippl::Comm->size());
    const size_t nranksU = static_cast<size_t>(nranks);

    // Mirror TrackRun: give each rank enough headroom above its ideal share.
    const size_t maxLocalNum = totalParticles / nranksU + 2 * nranksU + 1;

    pc->allocateParticles(maxLocalNum);
}

/// Standard deviation of a unit normal truncated per axis to |x| <= c.
static double truncatedNormalSd(double c) {
    const double phi = std::exp(-0.5 * c * c) / std::sqrt(2.0 * M_PI);
    const double Z   = std::erf(c / std::sqrt(2.0));  // 2*Phi(c) - 1
    return std::sqrt(1.0 - 2.0 * c * phi / Z);
}

/// Marginal standard deviation of one coordinate of a pair of unit normals (x, y)
/// that is accepted only inside the circle x^2 + y^2 <= c^2 (OPAL's elliptical cut).
/// r^2 is exponentially distributed with mean 2, so E[x^2] = E[r^2 | r^2 <= c^2] / 2.
static double ellipticalCutMarginalSd(double c) {
    const double q = std::exp(-0.5 * c * c);
    return std::sqrt(1.0 - 0.5 * c * c * q / (1.0 - q));
}

/// Mean of pz + (px^2 + py^2) / (2 pz): the quantity OPAL's adjustPhaseSpace
/// (and the Gaussian sampler) pins to the reference momentum avrgpz.
template <typename ViewType>
double computeMeanTaylorMomentum(ViewType& view, size_t nlocal, size_t total_nparticles) {
    double loc = 0.0, glob = 0.0;
    Kokkos::parallel_reduce(
            nlocal,
            KOKKOS_LAMBDA(const int k, double& s) {
                const double px = view(k)[0], py = view(k)[1], pz = view(k)[2];
                s += (pz != 0.0) ? pz + (px * px + py * py) / (2.0 * pz) : pz;
            },
            Kokkos::Sum<double>(loc));
    Kokkos::fence();
    MPI_Allreduce(&loc, &glob, 1, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());
    return glob / total_nparticles;
}

TEST_F(GaussianTest, meanR_stddevR) {
    const Vector_t<double, 3> sigmaR_ref = 0.5;
    const Vector_t<double, 3> sigmaP_ref = 1.0;
    double avrgpz                        = 0.0;
    const Vector_t<double, 3> cutoffR    = 4.0;

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

    // (x, y) are cut on an ellipse, z per axis; both are negligible at 4 sigma
    // but we compare against the exact truncated widths anyway.
    const double tol = 4.0 * sigmaR_ref[0] / std::sqrt(2.0 * total_nparticles);
    EXPECT_NEAR(sigmaR[0], sigmaR_ref[0] * ellipticalCutMarginalSd(cutoffR[0]), tol);
    EXPECT_NEAR(sigmaR[1], sigmaR_ref[1] * ellipticalCutMarginalSd(cutoffR[1]), tol);
    EXPECT_NEAR(sigmaR[2], sigmaR_ref[2] * truncatedNormalSd(cutoffR[2]), tol);
}

TEST_F(GaussianTest, cutoffR) {
    const Vector_t<double, 3> sigmaR  = 1.0;
    const Vector_t<double, 3> sigmaP  = 1.0;
    double avrgpz                     = 0.0;
    const Vector_t<double, 3> cutoffR = 3.0;
    const Vector_t<double, 3> cutoffP(2.0, 2.5, 3.0);

    bool fixMeanR = false;
    Gaussian sampler(pc, sigmaR, sigmaP, avrgpz, cutoffR, fixMeanR, cutoffP);

    size_t total_nparticles = 100000;

    preallocateParticleCapacity(pc, total_nparticles);
    sampler.generateParticles(total_nparticles, nr);

    auto Rview = pc->R.getView();

    // compute max of abs(R)
    double global_maxAbsR[3] = {0.0, 0.0, 0.0};
    double global_maxRadius  = 0.0;
    computeMaxAbsR(Rview, pc->getLocalNum(), global_maxAbsR, global_maxRadius);

    // --------- Assertions ---------

    // Each coordinate individually should be inside cutoff
    EXPECT_LE(global_maxAbsR[0], cutoffR[0]);
    EXPECT_LE(global_maxAbsR[1], cutoffR[1]);
    EXPECT_LE(global_maxAbsR[2], cutoffR[2]);

    // OPAL semantics: (x, y) and (px, py) are cut on ellipses, z and pz per axis.
    auto r             = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
    auto p             = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
    double maxEllipseR = 0.0, maxEllipseP = 0.0, maxAbsPz = 0.0;
    for (size_t i = 0; i < pc->getLocalNum(); ++i) {
        const double ex  = r(i)[0] / (sigmaR[0] * cutoffR[0]);
        const double ey  = r(i)[1] / (sigmaR[1] * cutoffR[1]);
        const double epx = p(i)[0] / (sigmaP[0] * cutoffP[0]);
        const double epy = p(i)[1] / (sigmaP[1] * cutoffP[1]);
        maxEllipseR      = std::max(maxEllipseR, ex * ex + ey * ey);
        maxEllipseP      = std::max(maxEllipseP, epx * epx + epy * epy);
        maxAbsPz         = std::max(maxAbsPz, std::abs(p(i)[2] - avrgpz));
    }
    EXPECT_LE(maxEllipseR, 1.0);
    EXPECT_LE(maxEllipseP, 1.0);
    EXPECT_LE(maxAbsPz, sigmaP[2] * cutoffP[2]);
    // With 1e5 samples the corners of the box |x|,|y| <= 3 sigma would be populated
    // by a per-axis truncation; the ellipse excludes them.
    EXPECT_GT(maxEllipseR, 0.95);
    EXPECT_GT(maxEllipseP, 0.95);
}

TEST_F(GaussianTest, meanP_and_steddevP) {
    const Vector_t<double, 3> sigmaR_ref = 1.0;
    const Vector_t<double, 3> sigmaP_ref = 0.01;
    double avrgpz                        = 1.0;
    const Vector_t<double, 3> cutoffR    = 4.0;
    const Vector_t<double, 3> cutoffP    = 3.0;

    Gaussian sampler(pc, sigmaR_ref, sigmaP_ref, avrgpz, cutoffR, true, cutoffP);

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
    // Transverse momenta are centred exactly (OPAL adjustPhaseSpace).
    EXPECT_NEAR(meanP[0], 0.0, 1e-15);
    EXPECT_NEAR(meanP[1], 0.0, 1e-15);
    // The average energy is matched: mean(pz + (px^2+py^2)/(2 pz)) == avrgpz up to the
    // second-order residual of the single linear correction step (same as OPAL),
    // ~ |eps| * sigma_t^2 / pz^2 ~ 1e-8 here. Hence mean(pz) lies below avrgpz by the
    // mean transverse contribution.
    const double meanTaylorP =
            computeMeanTaylorMomentum(Pview, pc->getLocalNum(), total_nparticles);
    EXPECT_NEAR(meanTaylorP, avrgpz, 1e-7);
    const double expectedPzShift = (sigmaP[0] * sigmaP[0] + sigmaP[1] * sigmaP[1]) / (2.0 * avrgpz);
    EXPECT_NEAR(meanP[2], avrgpz - expectedPzShift, 1e-6);
    EXPECT_LT(meanP[2], avrgpz);

    const double tol = 4.0 * sigmaP_ref[0] / std::sqrt(2.0 * total_nparticles);
    EXPECT_NEAR(sigmaP[0], sigmaP_ref[0] * ellipticalCutMarginalSd(cutoffP[0]), tol);
    EXPECT_NEAR(sigmaP[1], sigmaP_ref[1] * ellipticalCutMarginalSd(cutoffP[1]), tol);
    EXPECT_NEAR(sigmaP[2], sigmaP_ref[2] * truncatedNormalSd(cutoffP[2]), tol);
}

// A very large cutoff (what Distribution substitutes for CUTOFF* = 0, as in OPAL)
// disables the truncation: the sample is an untruncated normal.
TEST_F(GaussianTest, largeCutoffDisablesTruncation) {
    const Vector_t<double, 3> sigmaR = 1.0;
    const Vector_t<double, 3> sigmaP = 1.0;
    const double avrgpz              = 0.0;
    const double noCut               = std::numeric_limits<double>::max();
    const Vector_t<double, 3> cutoffR(noCut, noCut, 3.0);
    const Vector_t<double, 3> cutoffP(noCut, noCut, noCut);
    Gaussian sampler(pc, sigmaR, sigmaP, avrgpz, cutoffR, false, cutoffP);

    size_t total_nparticles = 100000;
    preallocateParticleCapacity(pc, total_nparticles);
    sampler.generateParticles(total_nparticles, nr);

    auto r           = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
    auto p           = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
    double locMax[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, globMax[6];
    for (size_t i = 0; i < pc->getLocalNum(); ++i) {
        for (unsigned d = 0; d < 3; ++d) {
            locMax[d]     = std::max(locMax[d], std::abs(r(i)[d]));
            locMax[3 + d] = std::max(locMax[3 + d], std::abs(p(i)[d]));
        }
        EXPECT_TRUE(std::isfinite(r(i)[0]) && std::isfinite(r(i)[1]) && std::isfinite(p(i)[2]));
    }
    MPI_Allreduce(locMax, globMax, 6, MPI_DOUBLE, MPI_MAX, ippl::Comm->getCommunicator());

    // P(|x| > 3 sigma) = 2.7e-3, i.e. ~270 of 1e5 samples lie beyond the former cut;
    // P(|x| > 5 sigma) = 5.7e-7, so the untruncated maximum stays below 5.5 sigma.
    for (unsigned d : {0u, 1u, 3u, 4u, 5u}) {
        EXPECT_GT(globMax[d], 3.0);
        EXPECT_LT(globMax[d], 5.5);
    }
    EXPECT_LE(globMax[2], cutoffR[2]);  // z keeps its per-axis cut

    auto Rview = pc->R.getView();
    auto Pview = pc->P.getView();
    double meanR[3], meanP[3], sdR[3], sdP[3];
    computeMean(Rview, pc->getLocalNum(), total_nparticles, meanR);
    computeMean(Pview, pc->getLocalNum(), total_nparticles, meanP);
    computeStdDev(Rview, pc->getLocalNum(), total_nparticles, meanR, sdR);
    computeStdDev(Pview, pc->getLocalNum(), total_nparticles, meanP, sdP);
    const double tol = 4.0 / std::sqrt(2.0 * total_nparticles);
    EXPECT_NEAR(sdR[0], 1.0, tol);
    EXPECT_NEAR(sdR[1], 1.0, tol);
    EXPECT_NEAR(sdP[0], 1.0, tol);
    EXPECT_NEAR(sdP[1], 1.0, tol);
    EXPECT_NEAR(sdP[2], 1.0, tol);
}

// Without the phase-space adjustment the sampled means are only statistically zero
// and pz is centred on avrgpz.
TEST_F(GaussianTest, noFixMean_keepsRawSample) {
    const Vector_t<double, 3> sigmaR = 1e-3;
    const Vector_t<double, 3> sigmaP = 0.01;
    const double avrgpz              = 1.0;
    Gaussian sampler(pc, sigmaR, sigmaP, avrgpz, Vector_t<double, 3>(3.0), false);

    size_t total_nparticles = 100000;
    preallocateParticleCapacity(pc, total_nparticles);
    sampler.generateParticles(total_nparticles, nr);

    auto Rview = pc->R.getView();
    auto Pview = pc->P.getView();
    double meanR[3], meanP[3];
    computeMean(Rview, pc->getLocalNum(), total_nparticles, meanR);
    computeMean(Pview, pc->getLocalNum(), total_nparticles, meanP);
    const double seR = sigmaR[0] / std::sqrt(double(total_nparticles));
    const double seP = sigmaP[0] / std::sqrt(double(total_nparticles));
    for (unsigned d = 0; d < 3; ++d) {
        EXPECT_NEAR(meanR[d], 0.0, 5.0 * seR);
        EXPECT_NEAR(meanP[d], d == 2 ? avrgpz : 0.0, 5.0 * seP);
    }
    EXPECT_GT(std::abs(meanR[0]) + std::abs(meanR[1]) + std::abs(meanR[2]), 0.0);
}

// Exercise every collapsed-axis combination, including a point bunch and
// N < MPI ranks. Inspect device-generated values through an explicit mirror.
TEST_F(GaussianTest, ZeroWidthAxesPreserveCountOffsetsAndFiniteMomenta) {
    for (unsigned mask = 1; mask < 8; ++mask) {
        for (size_t requested : {size_t(1), size_t(10001)}) {
            SCOPED_TRACE(
                    ::testing::Message() << "collapsed-axis mask=" << mask << ", N=" << requested);
            resetContainer();
            Vector_t<double, 3> sigma(1e-3);
            for (unsigned d = 0; d < 3; ++d)
                if (mask & (1u << d)) sigma[d] = 0.0;
            Gaussian sampler(pc, sigma, Vector_t<double, 3>(0.0), 0.25, Vector_t<double, 3>(4.0));
            const Vector_t<double, 3> offset(0.01, -0.02, 0.03);
            const Vector_t<double, 3> momentumOffset(0.02, -0.01, 0.03);
            sampler.setEmissionOffsets(offset, momentumOffset, 0.0);
            preallocateParticleCapacity(pc, requested);
            sampler.generateParticles(requested, nr);
            size_t local = pc->getLocalNum(), total = 0;
            ippl::Comm->allreduce(&local, &total, 1, std::plus<size_t>());
            EXPECT_EQ(total, requested);
            const size_t nranks = static_cast<size_t>(ippl::Comm->size());
            const size_t rank   = static_cast<size_t>(ippl::Comm->rank());
            EXPECT_EQ(local, requested / nranks + (rank < requested % nranks ? 1u : 0u));
            auto r = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
            auto p = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
            for (size_t i = 0; i < local; ++i) {
                for (unsigned d = 0; d < 3; ++d) {
                    EXPECT_TRUE(std::isfinite(r(i)[d]));
                    EXPECT_TRUE(std::isfinite(p(i)[d]));
                    EXPECT_EQ(p(i)[d], (d == 2 ? 0.25 : 0.0) + momentumOffset[d]);
                    if (mask & (1u << d)) EXPECT_EQ(r(i)[d], offset[d]);
                }
            }
            double mean[3];
            auto rv = pc->R.getView();
            computeMean(rv, local, requested, mean);
            // 1e-14 m accommodates floating-point summation of the source offsets.
            for (unsigned d = 0; d < 3; ++d)
                EXPECT_NEAR(mean[d], offset[d], 1e-14);
            if (requested > 1) {
                double sd[3];
                computeStdDev(rv, local, requested, mean, sd);
                for (unsigned d = 0; d < 3; ++d)
                    // A 5% width band is wider than the ~0.7% sampling uncertainty
                    // for N=10001, but excludes a collapsed or incorrectly scaled axis.
                    if (!(mask & (1u << d))) EXPECT_NEAR(sd[d], sigma[d], 5e-5);
            }
        }
    }
}

TEST_F(GaussianTest, ZeroWidthAppendPreservesExistingParticlesAndCutoff) {
    preallocateParticleCapacity(pc, 10003);
    pc->createParticles(1);
    pc->R = 17.0;
    pc->P = 23.0;
    const Vector_t<double, 3> sigma(0.0, 2e-3, 0.0);
    Gaussian sampler(pc, sigma, Vector_t<double, 3>(0.0), 0.0, Vector_t<double, 3>(2.0), false);
    size_t requested = 10001;
    sampler.generateParticles(requested, nr);
    auto r = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
    auto p = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
    for (unsigned d = 0; d < 3; ++d) {
        EXPECT_EQ(r(0)[d], 17.0);
        EXPECT_EQ(p(0)[d], 23.0);
    }
    for (size_t i = 1; i < pc->getLocalNum(); ++i) {
        EXPECT_EQ(r(i)[0], 0.0);
        EXPECT_EQ(r(i)[2], 0.0);
        EXPECT_TRUE(std::isfinite(r(i)[1]));
        EXPECT_LE(std::abs(r(i)[1]), 4e-3);
    }
    size_t local = pc->getLocalNum(), total = 0;
    ippl::Comm->allreduce(&local, &total, 1, std::plus<size_t>());
    EXPECT_EQ(total, requested + static_cast<size_t>(ippl::Comm->size()));
}

TEST_F(GaussianTest, OmittedSpatialWidthsCollapseWithoutSuppressingMomentumSpread) {
    Distribution distribution;
    Attributes::setPredefinedString(*distribution.findAttribute("TYPE"), "GAUSS");
    const Vector_t<double, 3> sigmaP(0.02, 0.03, 0.04);
    Attributes::setReal(*distribution.findAttribute("SIGMAPX"), sigmaP[0]);
    Attributes::setReal(*distribution.findAttribute("SIGMAPY"), sigmaP[1]);
    Attributes::setReal(*distribution.findAttribute("SIGMAPZ"), sigmaP[2]);
    distribution.execute();
    for (unsigned d = 0; d < 3; ++d)
        EXPECT_EQ(distribution.getSigmaR()[d], 0.0);

    Gaussian sampler(pc, &distribution);
    sampler.setAvrgpz(0.25);
    const Vector_t<double, 3> offset(0.01, -0.02, 0.03);
    const Vector_t<double, 3> momentumOffset(0.02, -0.01, 0.03);
    sampler.setEmissionOffsets(offset, momentumOffset, 0.0);
    size_t requested = 10001;
    preallocateParticleCapacity(pc, requested);
    sampler.generateParticles(requested, nr);

    const size_t local = pc->getLocalNum();
    auto r             = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
    auto p             = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
    for (size_t i = 0; i < local; ++i) {
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_EQ(r(i)[d], offset[d]);
            EXPECT_TRUE(std::isfinite(p(i)[d]));
        }
    }
    double meanP[3], sdP[3];
    auto pv = pc->P.getView();
    computeMean(pv, local, requested, meanP);
    computeStdDev(pv, local, requested, meanP, sdP);
    for (unsigned d = 0; d < 3; ++d) {
        // Transverse means are centred exactly; the mean of pz is below avrgpz by the
        // transverse contribution removed when matching the average energy.
        double expected = momentumOffset[d];
        if (d == 2) {
            expected += 0.25 - (sdP[0] * sdP[0] + sdP[1] * sdP[1]) / (2.0 * 0.25);
        }
        // sigma_t / pz = 0.12 here, so the second-order residual of the linear
        // energy correction is ~1e-4 relative.
        EXPECT_NEAR(meanP[d], expected, d == 2 ? 1e-4 : 1e-15);
        // 5% for the sampled (3 sigma truncated) standard deviation.
        EXPECT_NEAR(sdP[d], sigmaP[d], 0.05 * sigmaP[d]);
    }
    // Energy matching holds (to second order) once the emission offsets are removed again.
    auto pShape      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
    double sumTaylor = 0.0;
    for (size_t i = 0; i < local; ++i) {
        const double px = pShape(i)[0] - momentumOffset[0];
        const double py = pShape(i)[1] - momentumOffset[1];
        const double pz = pShape(i)[2] - momentumOffset[2];
        sumTaylor += pz + (px * px + py * py) / (2.0 * pz);
    }
    double globTaylor = 0.0;
    MPI_Allreduce(&sumTaylor, &globTaylor, 1, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());
    EXPECT_NEAR(globTaylor / requested, 0.25, 1e-4);
}
