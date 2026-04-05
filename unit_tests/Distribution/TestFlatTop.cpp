#include <gtest/gtest.h>
#include <mpi.h>
#include <memory>
#include <cmath>
#include <algorithm>

#include "Distribution/FlatTop.h"
#include "PartBunch/BunchStateHandler.h"
#include "Ippl.h"
#include "Utility/IpplTimings.h"

class FlatTopTest : public ::testing::Test {
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
        std::array<bool,3> decomp = {false, false, true};

        ippl::NDIndex<3> domain;
        for (unsigned i = 0; i < 3; i++) {
            domain[i] = ippl::Index(this->nr[i]);
        }

        // Create FieldContainer
        fc = std::make_shared<FieldContainer_t>(hr, rmin, rmax, decomp, domain, origin, this->isAllPeriodic_m);

        // Create mesh and fieldlayout
        Mesh_t<3> mesh(domain, hr, origin);
        FieldLayout_t<3> fl(MPI_COMM_WORLD, domain, decomp, this->isAllPeriodic_m);

        // Create ParticleContainer
        pc = std::make_shared<ParticleContainer<double, 3>>(mesh, fl);

        bunchStateHandler = std::make_shared<BunchStateHandler>();
    }

    void TearDown() override {
        // nothing special
    }

    std::shared_ptr<ParticleContainer<double, 3>> pc;
    std::shared_ptr<FieldContainer_t> fc;
    std::shared_ptr<BunchStateHandler> bunchStateHandler;
    ippl::Vector<int,3> nr;
    bool isAllPeriodic_m = true;
};

/**
 * @brief Testing generatiioin of uniform sampples on an elliptical disk.
 *
 * This test verifies that `FlatTop::generateUniformDisk` generates particle
 * positions uniformly over the elliptical region
 *
 * @f[
 * \mathcal{D} = \left\{ (x,y) \in \mathbb{R}^2 :
 * \frac{R_x^2}{\sigma_x^2} + \frac{R_y^2}{\sigma_y^2} \le 1 \right\}.
 * @f]
 *
 * The test checks the following properties:
 *
 * - **Support**:
 *   All generated particles satisfy
 *   @f[
 *   \frac{R_x^2}{\sigma_x^2} + \frac{R_y^2}{\sigma_y^2} \le 1.
 *   @f]
 *
 * - **Zero momentum**:
 *   @f[
 *   \mathbf{P}_i = \mathbf{0}.
 *   @f]
 *
 * - **First moments**:
 *   By symmetry of the uniform disk,
 *   @f[
 *   \mathbb{E}[R_x] = \mathbb{E}[R_y] = 0.
 *   @f]
 *
 * - **Second moments**:
 *   For a uniform distribution over an elliptical disk,
 *   @f[
 *   \mathbb{E}[R_x^2] = \frac{\sigma_x^2}{4},
 *   \qquad
 *   \mathbb{E}[R_y^2] = \frac{\sigma_y^2}{4}.
 *   @f]
 *
 * The empirical sample moments are compared against these analytic values
 * using a finite-sample tolerance.
 */
TEST_F(FlatTopTest, UniformDiskStatisticsAndBounds) {
    const Vector_t<double,3> sigmaR = {0.5, 1.0, 0.0};
    const Vector_t<double,3> cutoff = 4.0;

    // Minimal FlatTop constructor usage
    FlatTop sampler(
        pc,
        /*emitting=*/false,
        /*sigmaTFall=*/1.0,
        /*sigmaTRise=*/1.0,
        cutoff,
        /*tPulseLengthFWHM=*/1.0,
        sigmaR
    );
    sampler.setBunchStateHandler(bunchStateHandler);

    const size_t nlocal = 0;
    const size_t nNew   = 100000;

    pc->create(nNew);
    sampler.generateUniformDisk(nlocal, nNew, 1.0e-12);

    auto Rview_d = pc->R.getView();
    auto Pview_d = pc->P.getView();

    auto Rview = Kokkos::create_mirror_view(Rview_d);
    auto Pview = Kokkos::create_mirror_view(Pview_d);
    Kokkos::deep_copy(Rview, Rview_d);
    Kokkos::deep_copy(Pview, Pview_d);

    double sumx = 0.0, sumy = 0.0;
    double sumx2 = 0.0, sumy2 = 0.0;

    for (size_t i = 0; i < nNew; ++i) {
        const double x = Rview(i)[0];
        const double y = Rview(i)[1];
        const double z = Rview(i)[2];

        // z-position must be zero
        EXPECT_DOUBLE_EQ(z, 0.0);

        // momentum must be zero
        EXPECT_DOUBLE_EQ(Pview(i)[0], 0.0);
        EXPECT_DOUBLE_EQ(Pview(i)[1], 0.0);
        EXPECT_DOUBLE_EQ(Pview(i)[2], 0.0);

        // inside ellipse
        const double r2 =
            (x*x)/(sigmaR[0]*sigmaR[0]) +
            (y*y)/(sigmaR[1]*sigmaR[1]);

        EXPECT_LE(r2, 1.0 + 1e-14);

        sumx  += x;
        sumy  += y;
        sumx2 += x*x;
        sumy2 += y*y;
    }

    const double meanx = sumx / nNew;
    const double meany = sumy / nNew;
    const double varx  = sumx2 / nNew;
    const double vary  = sumy2 / nNew;

    // Mean should be ~0
    EXPECT_NEAR(meanx, 0.0, 1e-2);
    EXPECT_NEAR(meany, 0.0, 1e-2);

    // Uniform disk second moments:
    // E[x^2] = σx^2 / 4, E[y^2] = σy^2 / 4
    EXPECT_NEAR(varx, sigmaR[0]*sigmaR[0] / 4.0, 1e-2);
    EXPECT_NEAR(vary, sigmaR[1]*sigmaR[1] / 4.0, 1e-2);
}


/**
 * @brief This test counts the number of particles entering the domain per rank for a FlatTop pulse.
 *
 * This test verifies that `FlatTop::countEnteringParticlesPerRank` correctly
 * counts the number of particles entering the domain within a given time window
 * when domain decomposition is disabled.
 *
 * A total of @f$N_{\mathrm{tot}}@f$ particles are emitted over the full pulse
 * duration according to a normalized FlatTop temporal profile
 * @f$f(t)@f$. The expected number of particles entering during a time window
 * @f$[t_0, t_f]@f$ is computed analytically.
 *
 * The number of particles emitted in the interval @f$[t_0, t_f]@f$ is
 *
 * @f[
 * N_{\mathrm{new}}
 * =
 * \left\lfloor
 * N_{\mathrm{tot}}
 * \frac{
 * \int_{t_0}^{t_f} f(t)\,dt
 * }{
 * \int f(t)\,dt
 * }
 * \right\rfloor.
 * @f]
 *
 * The integral over the time window is approximated using a trapezoidal rule:
 *
 * @f[
 * \int_{t_0}^{t_f} f(t)\,dt
 * \approx
 * \frac{1}{2}
 * \left( f(t_0) + f(t_f) \right)
 * (t_f - t_0).
 * @f]
 *
 * With domain decomposition disabled, particles are distributed evenly across
 * all MPI ranks. The expected number of particles per rank is therefore
 *
 * @f[
 * N_{\mathrm{new,local}}}
 * =
 * \left\lfloor
 * \frac{N_{\mathrm{new}}}{{\mathrm{\#ranks}}}
 * \right\rfloor,
 * @f]
 *
 * with the remainder
 *
 * @f[
 * remainder = N_{\mathrm{new}} - \#{\mathrm{rank}} N_{\mathrm{new,local}}
 * @f]
 *
 * assigned to rank 0, leading to
 *
 * @f[
 * N_{\mathrm{new,local}}^{(0)} = N_{\mathrm{new,local}} + remainder.
 * @f]
 *
 * and @f$N_{\mathrm{new,local}}^{(r)} = N_{\mathrm{new,local}}@f$ for @f$r > 0@f$.
 * 
 * The total number of particles entering the domain is therefore
 *
 * @f[
 * N_{\mathrm{new,global}}
 * =
 * \sum_{r=0}^{N_{\mathrm{ranks}}-1}
 * N_{\mathrm{new,local}}^{(r)}.
 * @f]
 *
 * The test verifies that
 *
 * @f[
 * N_{\mathrm{new,global}} = N_{\mathrm{expected}},
 * @f]
 *
 * The test verifies that the number of particles counted on each rank matches
 * this analytic expectation exactly.
 */
TEST_F(FlatTopTest, CountEnteringParticles_NoDomainDecomp) {
    const Vector_t<double,3> sigmaR = {1.0, 1.0, 0.0};
    const Vector_t<double,3> cutoff = {4.0, 4.0, 4.0};

    FlatTop sampler(
        pc,
        /*emitting=*/true,
        /*sigmaTFall=*/1.0,
        /*sigmaTRise=*/1.0,
        cutoff,
        /*tPulseLengthFWHM=*/10.0,
        sigmaR
    );
    sampler.setBunchStateHandler(bunchStateHandler);

    sampler.setWithDomainDecomp(false);

    // total number of particles to emit over whole pulse
    const size_t totalN = 100000;
    sampler.allocateParticles(totalN);

    // Preallocate capacity so computeLocalEmitCount can distribute totalN.
    const int nranksConst     = std::max(1, ippl::Comm->size());
    const size_t nranksU = static_cast<size_t>(nranksConst);
    const size_t maxLocalNum =
        totalN / nranksU + 2 * nranksU + 1;
    pc->create(maxLocalNum);
    Kokkos::View<bool*> tmp_invalid("tmp_invalid", maxLocalNum);
    pc->destroy(tmp_invalid, maxLocalNum);

    const double t0 = 2.0;
    const double tf = 2.1;

    size_t nlocal = sampler.countEnteringParticlesPerRank(t0, tf);

    // --- compute expected result ---
    double f0 = sampler.FlatTopProfile(t0);
    double f1 = sampler.FlatTopProfile(tf);
    double tArea = 0.5 * (f0 + f1) * (tf - t0);

    double expectedTotalNewD =
        std::floor(totalN * tArea / sampler.getDistArea());
    const size_t expectedTotalNew =
        static_cast<size_t>(expectedTotalNewD);

    int nranks;
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Mirror SamplingBase::computeLocalEmitCount for equal capacities:
    // first 'rem' ranks get one extra particle.
    size_t base = expectedTotalNew / static_cast<size_t>(nranks);
    size_t rem  = expectedTotalNew % static_cast<size_t>(nranks);
    size_t expectedLocal = base + ((static_cast<size_t>(rank) < rem) ? 1 : 0);

    EXPECT_EQ(nlocal, expectedLocal);

    size_t globalEmitted = 0;
    MPI_Allreduce(
        &nlocal,
        &globalEmitted,
        1,
        MPI_UNSIGNED_LONG,
        MPI_SUM,
        MPI_COMM_WORLD
    );

    EXPECT_EQ(globalEmitted, expectedTotalNew);
}