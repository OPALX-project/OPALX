#include <gtest/gtest.h>
#include <mpi.h>
#include <memory>
#include <cmath>

#include "Distribution/FlatTop.h"
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
    }

    void TearDown() override {
        // nothing special
    }

    std::shared_ptr<ParticleContainer<double, 3>> pc;
    std::shared_ptr<FieldContainer_t> fc;
    ippl::Vector<int,3> nr;
    bool isAllPeriodic_m = true;
};

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

    const size_t nlocal = 0;
    const size_t nNew   = 100000;

    pc->create(nNew);
    sampler.generateUniformDisk(nlocal, nNew);

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