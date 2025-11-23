#include <gtest/gtest.h>
#include <mpi.h>
#include <memory>
#include <cmath>

#include "Distribution/Gaussian.h"
#include "Distribution/Distribution.h"
#include "PartBunch/FieldContainer.hpp"
#include "PartBunch/ParticleContainer.hpp"

using ParticleContainer_t = ParticleContainer<double, 3>;
using FieldContainer_t = FieldContainer<double, 3>;
using Distribution_t = Distribution;

class GaussianTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        // Initialize MPI if not already initialized
        int mpi_initialized = 0;
        MPI_Initialized(&mpi_initialized);
        if (!mpi_initialized) {
            MPI_Init(nullptr, nullptr);
            mpi_init_here = true;
        }

        // Initialize Kokkos (required before creating any Kokkos::View)
        Kokkos::initialize();
        kokkos_initialized_here = true;
    }

    static void TearDownTestSuite() {
        // Finalize Kokkos
        if (kokkos_initialized_here) {
            Kokkos::finalize();
            kokkos_initialized_here = false;
        }

        // Finalize MPI if we initialized it
        if (mpi_init_here) {
            MPI_Finalize();
            mpi_init_here = false;
        }
    }

    void SetUp() override {
        // Minimal 3D grid parameters
        ippl::Vector<double,3> hr;     hr(0) = 1.0; hr(1) = 1.0; hr(2) = 1.0;
        ippl::Vector<double,3> rmin;   rmin(0) = -10.0; rmin(1) = -10.0; rmin(2) = -10.0;
        ippl::Vector<double,3> rmax;   rmax(0) = 10.0;  rmax(1) = 10.0;  rmax(2) = 10.0;
        ippl::Vector<double,3> origin; origin(0) = 0.0;  origin(1) = 0.0;  origin(2) = 0.0;
        std::array<bool,3> decomp = {false, false, false};

        ippl::NDIndex<3> domain;  // default constructor
        domain[0] = 4;
        domain[1] = 4;
        domain[2] = 4;

        // Create FieldContainer
        fc = std::make_shared<FieldContainer_t>(hr, rmin, rmax, decomp, domain, origin, true);

        // Create mesh and fieldlayout
        Mesh_t<3> mesh(domain, hr, origin);
        FieldLayout_t<3> fl(MPI_COMM_WORLD, domain, decomp, true);

        // Create ParticleContainer
        pc = std::make_shared<ParticleContainer_t>(mesh, fl);

        // Create Distribution
        dist = std::make_shared<Distribution_t>();
    }

    void TearDown() override {
        // nothing special
    }

    static inline bool mpi_init_here = false;
    static inline bool kokkos_initialized_here = false;

    std::shared_ptr<ParticleContainer_t> pc;
    std::shared_ptr<FieldContainer_t> fc;
    std::shared_ptr<Distribution_t> dist;
};

TEST_F(GaussianTest, GenerateParticlesBasic) {
    Gaussian sampler(pc, fc, dist);

    size_t N = 1000;
    ippl::Vector<double,3> nr; nr(0)=1; nr(1)=1; nr(2)=1;

    sampler.generateParticles(N, nr);

    auto Rview = pc->R.getView();
    auto Pview = pc->P.getView();

    auto R_h = Kokkos::create_mirror_view(Rview);
    auto P_h = Kokkos::create_mirror_view(Pview);
    Kokkos::deep_copy(R_h, Rview);
    Kokkos::deep_copy(P_h, Pview);

    size_t nlocal = R_h.extent(0);
    ASSERT_GT(nlocal, 0u);

    double meanRx=0, meanRy=0, meanRz=0, meanPz=0;
    for(size_t i=0;i<nlocal;i++){
        meanRx += R_h(i)[0];
        meanRy += R_h(i)[1];
        meanRz += R_h(i)[2];
        meanPz  += P_h(i)[2];
    }
    meanRx /= nlocal;
    meanRy /= nlocal;
    meanRz /= nlocal;
    meanPz  /= nlocal;

    EXPECT_NEAR(meanRx, 0.0, 1e-1);
    EXPECT_NEAR(meanRy, 0.0, 1e-1);
    EXPECT_NEAR(meanRz, 0.0, 1e-1);

    double avrgpz = dist->getAvrgpz();
    EXPECT_NEAR(meanPz, avrgpz, 1e-1);
}
