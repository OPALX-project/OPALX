#include <gtest/gtest.h>

#include "Distribution/Uniform.h"
#include "PartBunch/BunchStateHandler.h"

namespace {

    void checkColdUniformEllipsoid(const std::shared_ptr<ParticleContainer_t>& pc) {
        constexpr size_t total = 100000;
        const Vector_t<double, 3> axes(1.0, 2.0, 3.0);
        constexpr double pz = 0.002;

        pc->allocateParticles(total + 1);
        Uniform sampler(pc, axes, pz);
        size_t n = total;
        sampler.generateParticles(n, Vector_t<double, 3>(16));

        const auto R = pc->R.getView();
        const auto P = pc->P.getView();
        Vector_t<double, 3> sum(0.0);
        Vector_t<double, 3> sum2(0.0);
        double maxEllipsoidRadius = 0.0;
        double maxMomentumError   = 0.0;
        Kokkos::parallel_reduce(
                pc->getLocalNum(),
                KOKKOS_LAMBDA(
                        const size_t i, Vector_t<double, 3>& localSum,
                        Vector_t<double, 3>& localSum2, double& maxRadius, double& maxPError) {
                    localSum += R(i);
                    localSum2 += R(i) * R(i);
                    const double radius = R(i)[0] * R(i)[0] / (axes[0] * axes[0])
                                          + R(i)[1] * R(i)[1] / (axes[1] * axes[1])
                                          + R(i)[2] * R(i)[2] / (axes[2] * axes[2]);
                    maxRadius = Kokkos::max(maxRadius, radius);
                    maxPError = Kokkos::max(
                            maxPError, Kokkos::abs(P(i)[0]) + Kokkos::abs(P(i)[1])
                                               + Kokkos::abs(P(i)[2] - pz));
                },
                Kokkos::Sum<Vector_t<double, 3>>(sum), Kokkos::Sum<Vector_t<double, 3>>(sum2),
                Kokkos::Max<double>(maxEllipsoidRadius), Kokkos::Max<double>(maxMomentumError));

        ippl::Comm->allreduce(sum[0], 3, std::plus<double>());
        ippl::Comm->allreduce(sum2[0], 3, std::plus<double>());
        ippl::Comm->allreduce(maxEllipsoidRadius, 1, std::greater<double>());
        ippl::Comm->allreduce(maxMomentumError, 1, std::greater<double>());

        sum /= static_cast<double>(total);
        sum2 /= static_cast<double>(total);
        EXPECT_LT(maxEllipsoidRadius, 1.0 + 1.0e-12);
        EXPECT_DOUBLE_EQ(maxMomentumError, 0.0);
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_NEAR(sum[d], 0.0, 1.5e-2 * axes[d]);
            EXPECT_NEAR(sum2[d], axes[d] * axes[d] / 5.0, 1.0e-2 * axes[d] * axes[d]);
        }
    }

}  // namespace

class UniformTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() { ippl::finalize(); }

    void SetUp() override {
        const Vector_t<int, 3> nr(16);
        const Vector_t<double, 3> rmin(-1.0);
        const Vector_t<double, 3> hr(2.0 / 16.0);
        ippl::NDIndex<3> domain;
        for (unsigned d = 0; d < 3; ++d) {
            domain[d] = ippl::Index(nr[d]);
        }
        std::array<bool, 3> decomp{true, true, true};
        Mesh_t<3> mesh(domain, hr, rmin);
        FieldLayout_t<3> layout(MPI_COMM_WORLD, domain, decomp, false);
        pc = std::make_shared<ParticleContainer_t>(mesh, layout);
        pc->setBunchStateHandler(std::make_shared<BunchStateHandler>());
    }

    std::shared_ptr<ParticleContainer_t> pc;
};

TEST_F(UniformTest, SamplesColdUniformEllipsoid) { checkColdUniformEllipsoid(pc); }
