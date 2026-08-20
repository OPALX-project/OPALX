#include <gtest/gtest.h>

#include "BeamlineCore/ConstantFocusingRep.h"
#include "PartBunch/BunchStateHandler.h"
#include "PartBunch/ParticleContainer.hpp"

class ConstantFocusingTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() { ippl::finalize(); }

    void SetUp() override {
        const Vector_t<int, 3> nr(8);
        const Vector_t<double, 3> rmin(-1.0);
        const Vector_t<double, 3> hr(0.25);
        ippl::NDIndex<3> domain;
        for (unsigned d = 0; d < 3; ++d) {
            domain[d] = ippl::Index(nr[d]);
        }
        std::array<bool, 3> decomp{true, true, true};
        ippl::UniformCartesian<double, 3> mesh(domain, hr, rmin);
        ippl::FieldLayout<3> layout(MPI_COMM_WORLD, domain, decomp, false);
        pc = std::make_shared<ParticleContainer_t>(mesh, layout);
        pc->setBunchStateHandler(std::make_shared<BunchStateHandler>());
        pc->allocateParticles(3);
        pc->createParticles(2);
    }

    std::shared_ptr<ParticleContainer_t> pc;
};

TEST_F(ConstantFocusingTest, UsesInitialSpaceChargeScaleAboutCentroid) {
    auto R     = pc->R.getView();
    auto E     = pc->E.getView();
    auto hostR = Kokkos::create_mirror_view(R);
    auto hostE = Kokkos::create_mirror_view(E);
    hostR(0)   = Vector_t<double, 3>(1.0, 2.0, 3.0);
    hostR(1)   = Vector_t<double, 3>(3.0, 4.0, 5.0);
    hostE(0)   = Vector_t<double, 3>(1.0, -2.0, 2.0);
    hostE(1)   = Vector_t<double, 3>(-1.0, 2.0, -2.0);
    Kokkos::deep_copy(R, hostR);
    Kokkos::deep_copy(E, hostE);

    ConstantFocusingRep focusing("focus");
    focusing.getGeometry().setElementLength(6.0);
    focusing.setRadius(2.0);
    focusing.setStrength(3.0);
    focusing.apply(pc);

    Kokkos::deep_copy(hostE, E);
    EXPECT_DOUBLE_EQ(focusing.getGradient(), 4.5);
    const Vector_t<double, 3> initialE0(1.0, -2.0, 2.0);
    const Vector_t<double, 3> initialE1(-1.0, 2.0, -2.0);
    for (unsigned d = 0; d < 3; ++d) {
        EXPECT_DOUBLE_EQ(hostE(0)[d], initialE0[d] - 4.5);
        EXPECT_DOUBLE_EQ(hostE(1)[d], initialE1[d] + 4.5);
    }
}

TEST_F(ConstantFocusingTest, AppliesOnlyInsideLongitudinalElementExtent) {
    pc->createParticles(1);

    auto R     = pc->R.getView();
    auto E     = pc->E.getView();
    auto hostR = Kokkos::create_mirror_view(R);
    auto hostE = Kokkos::create_mirror_view(E);

    // Only the middle particle has local z in [0, 1); the upper boundary is excluded.
    hostR(0) = Vector_t<double, 3>(-1.0, 0.0, -0.25);
    hostR(1) = Vector_t<double, 3>(2.0, 0.0, 0.25);
    hostR(2) = Vector_t<double, 3>(-1.0, 0.0, 1.0);
    hostE(0) = Vector_t<double, 3>(1.0, 0.0, 0.0);
    hostE(1) = Vector_t<double, 3>(1.0, 0.0, 0.0);
    hostE(2) = Vector_t<double, 3>(1.0, 0.0, 0.0);
    Kokkos::deep_copy(R, hostR);
    Kokkos::deep_copy(E, hostE);

    ConstantFocusingRep focusing("focus");
    focusing.getGeometry().setElementLength(1.0);
    focusing.setRadius(1.0);
    focusing.setStrength(1.0);
    focusing.apply(pc);

    Kokkos::deep_copy(hostE, E);
    EXPECT_DOUBLE_EQ(focusing.getGradient(), 1.0);
    const Vector_t<double, 3> expectedOutside(1.0, 0.0, 0.0);
    const Vector_t<double, 3> expectedInside(3.0, 0.0, -1.0 / 12.0);
    for (unsigned d = 0; d < 3; ++d) {
        EXPECT_DOUBLE_EQ(hostE(0)[d], expectedOutside[d]);
        EXPECT_DOUBLE_EQ(hostE(1)[d], expectedInside[d]);
        EXPECT_DOUBLE_EQ(hostE(2)[d], expectedOutside[d]);
    }
}
