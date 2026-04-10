#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>

#include "Ippl.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Processes/GlobalProcesses/Decay.h"
#include "Utilities/Options.h"

namespace {

using PC_t = ParticleContainer<double, 3>;

class DecayTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() { ippl::finalize(); }

    void SetUp() override {
        oldSeed_m            = Options::seed;
        oldUseQMAttributes_m = Options::useQMAttributes;
        Options::useQMAttributes = false;
    }

    void TearDown() override {
        Options::seed            = oldSeed_m;
        Options::useQMAttributes = oldUseQMAttributes_m;
    }

    std::shared_ptr<PC_t> makeContainer() {
        ippl::Vector<int, 3> nr        = 8;
        ippl::Vector<double, 3> rmin   = -4.0;
        ippl::Vector<double, 3> rmax   = 4.0;
        ippl::Vector<double, 3> origin = rmin;
        ippl::Vector<double, 3> hr     = (rmax - rmin) / ippl::Vector<double, 3>(nr);
        std::array<bool, 3> decomp     = {true, true, true};

        ippl::NDIndex<3> domain;
        for (unsigned i = 0; i < 3; ++i) {
            domain[i] = ippl::Index(nr[i]);
        }

        Mesh_t<3> mesh(domain, hr, origin);
        FieldLayout_t<3> fl(MPI_COMM_WORLD, domain, decomp, true);
        return std::make_shared<PC_t>(mesh, fl);
    }

    void createParticles(std::shared_ptr<PC_t>& pc, size_t nPart, double pz) {
        if (nPart == 0) {
            return;
        }

        pc->create(nPart);
        auto R_host = pc->R.getHostMirror();
        auto P_host = pc->P.getHostMirror();

        for (size_t i = 0; i < nPart; ++i) {
            R_host(i)[0] = 0.0;
            R_host(i)[1] = 0.0;
            R_host(i)[2] = 0.0;
            P_host(i)[0] = 0.0;
            P_host(i)[1] = 0.0;
            P_host(i)[2] = pz;
        }

        Kokkos::deep_copy(pc->R.getView(), R_host);
        Kokkos::deep_copy(pc->P.getView(), P_host);
        Kokkos::fence();
    }

    std::pair<size_t, size_t> runDecay(size_t nPart,
                                       double pz,
                                       int seed,
                                       size_t containerIndex,
                                       double tau0,
                                       double dt) {
        auto pc = makeContainer();
        createParticles(pc, nPart, pz);
        const size_t before = pc->getTotalNum();

        Options::seed = seed;
        Decay decay(tau0, containerIndex);
        const size_t destroyed = decay.apply(*pc, dt, 0, containerIndex);
        const size_t after     = pc->getTotalNum();

        EXPECT_EQ(before, after + destroyed);
        return {destroyed, after};
    }

private:
    int oldSeed_m = -1;
    bool oldUseQMAttributes_m = false;
};

TEST_F(DecayTest, ApplyReturnsZeroForNonPositiveDt) {
    auto pc = makeContainer();
    createParticles(pc, 64, 0.0);
    const size_t before = pc->getTotalNum();

    Decay decay(1.0, 0);
    EXPECT_EQ(decay.apply(*pc, 0.0, 7, 2), 0u);
    EXPECT_EQ(decay.apply(*pc, -1.0, 7, 2), 0u);
    EXPECT_EQ(pc->getTotalNum(), before);
}

TEST_F(DecayTest, ApplyReturnsZeroForInvalidLifetime) {
    const std::array<double, 4> invalidTau = {
        0.0,
        -1.0,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity()};

    for (double tau0 : invalidTau) {
        auto pc = makeContainer();
        createParticles(pc, 32, 0.0);
        const size_t before = pc->getTotalNum();

        Decay decay(tau0, 1);
        EXPECT_EQ(decay.apply(*pc, 1.0, 0, 1), 0u);
        EXPECT_EQ(pc->getTotalNum(), before);
    }
}

TEST_F(DecayTest, ApplyReturnsZeroForEmptyContainer) {
    auto pc = makeContainer();
    ASSERT_EQ(pc->getLocalNum(), 0u);

    Decay decay(1.0, 0);
    EXPECT_EQ(decay.apply(*pc, 1.0, 0, 0), 0u);
    EXPECT_EQ(pc->getTotalNum(), 0u);
}

TEST_F(DecayTest, LargeDtDecaysAllRestParticles) {
    auto [destroyed, remaining] = runDecay(
        256, 0.0, 1234, 0, 1.0e-12, 1.0);

    EXPECT_EQ(destroyed, 256u);
    EXPECT_EQ(remaining, 0u);
}

TEST_F(DecayTest, HigherMomentumDecaysFewerParticles) {
    constexpr size_t nPart = 2000;
    auto [lowMomentumDestroyed, lowMomentumRemaining] = runDecay(
        nPart, 0.0, 77, 3, 1.0, 1.0);
    auto [highMomentumDestroyed, highMomentumRemaining] = runDecay(
        nPart, 10.0, 77, 3, 1.0, 1.0);

    EXPECT_GT(lowMomentumDestroyed, highMomentumDestroyed);
    EXPECT_GT(lowMomentumDestroyed, nPart / 3);
    EXPECT_LT(highMomentumDestroyed, nPart / 3);
    EXPECT_LT(lowMomentumRemaining, highMomentumRemaining);
}

TEST_F(DecayTest, SameSeedAndInputsAreReproducible) {
    auto [firstDestroyed, firstRemaining] = runDecay(
        512, 0.5, 98765, 4, 2.0, 0.5);
    auto [secondDestroyed, secondRemaining] = runDecay(
        512, 0.5, 98765, 4, 2.0, 0.5);

    EXPECT_EQ(firstDestroyed, secondDestroyed);
    EXPECT_EQ(firstRemaining, secondRemaining);
}

}  // namespace
