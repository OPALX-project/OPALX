#include <gtest/gtest.h>
#include <mpi.h>
#include <memory>
#include <cmath>

#include "Distribution/Gaussian.h"
#include "Ippl.h"
#include "Utility/IpplTimings.h"

template <typename T>
using ParticleAttrib = ippl::ParticleAttrib<T>;

using size_type = ippl::detail::size_type;

// a minimal ParticleContainer definition for testing
template <typename T, unsigned Dim>
class ParticleContainer_min : public ParticleContainer<T, Dim> {
    using Base = ParticleContainer<T, Dim>;
public:
    /// particle momenta [\beta\gamma]
    typename Base::particle_position_type P;

    ParticleContainer_min(Mesh_t<Dim>& mesh, FieldLayout_t<Dim>& FL)
        : Base(mesh, FL)  // initialize the base class
    {
        registerAttributes();
    }

    ~ParticleContainer_min() = default;

    void registerAttributes() {
        this->addAttribute(P);
    }
};

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
        pc = std::make_shared<ParticleContainer_min<double, 3>>(mesh, fl);
    }

    void TearDown() override {
        // nothing special
    }

    std::shared_ptr<ParticleContainer_min<double, 3>> pc;
    ippl::Vector<int,3> nr;
    bool isAllPeriodic_m = true;
};

TEST_F(GaussianTest, meanR) {
    const Vector_t<double, 3> sigmaR = 1.0;
    const Vector_t<double, 3> sigmaP = 1.0;
    double avrgpz = 0.0;
    const Vector_t<double, 3> cutoffR = 3.0;

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
}
