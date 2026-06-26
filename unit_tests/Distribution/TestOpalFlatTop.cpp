#include <gtest/gtest.h>
#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

#include "Distribution/OpalFlatTop.h"
#include "Ippl.h"
#include "PartBunch/BunchStateHandler.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

class OpalFlatTopTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() { ippl::finalize(); }

    void SetUp() override {
        nr                             = 32;
        ippl::Vector<double, 3> rmin   = -4.0;
        ippl::Vector<double, 3> rmax   = 4.0;
        ippl::Vector<double, 3> origin = rmin;
        ippl::Vector<double, 3> hr     = (rmax - rmin) / ippl::Vector<double, 3>(nr);
        std::array<bool, 3> decomp     = {false, false, true};

        ippl::NDIndex<3> domain;
        for (unsigned i = 0; i < 3; ++i) {
            domain[i] = ippl::Index(this->nr[i]);
        }

        fc = std::make_shared<FieldContainer_t>(
                hr, rmin, rmax, decomp, domain, origin, this->isAllPeriodic_m);

        Mesh_t<3> mesh(domain, hr, origin);
        FieldLayout_t<3> fl(MPI_COMM_WORLD, domain, decomp, this->isAllPeriodic_m);
        pc = std::make_shared<ParticleContainer<double, 3>>(mesh, fl);

        bunchStateHandler = std::make_shared<BunchStateHandler>();
        pc->setBunchStateHandler(bunchStateHandler);
    }

    size_t globalLocalNum() const {
        size_t local  = pc->getLocalNum();
        size_t global = 0;
        MPI_Allreduce(&local, &global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        return global;
    }

    std::shared_ptr<ParticleContainer<double, 3>> pc;
    std::shared_ptr<FieldContainer_t> fc;
    std::shared_ptr<BunchStateHandler> bunchStateHandler;
    ippl::Vector<int, 3> nr;
    bool isAllPeriodic_m = true;
};

TEST_F(OpalFlatTopTest, BuildsSortedInventoryWithExactSizeAndEmissionBounds) {
    const Vector_t<double, 3> sigmaR = {0.5, 0.5, 0.0};
    const Vector_t<double, 3> cutoff = {4.0, 4.0, 4.0};
    OpalFlatTop sampler(
            pc, fc,
            /*emitting=*/true,
            /*sigmaTFall=*/0.15,
            /*sigmaTRise=*/0.25, cutoff,
            /*tPulseLengthFWHM=*/1.2, sigmaR);

    size_t totalParticles = 1001;
    const size_t nranks   = static_cast<size_t>(std::max(1, ippl::Comm->size()));
    pc->allocateParticles(totalParticles / nranks + 2 * nranks + 16);

    sampler.generateParticles(totalParticles, nr);

    const auto& birthTimes = sampler.getBirthTimes();
    ASSERT_EQ(birthTimes.size(), totalParticles);
    EXPECT_TRUE(std::is_sorted(birthTimes.begin(), birthTimes.end()));
    EXPECT_GE(birthTimes.front(), -0.5 * sampler.getEmissionTime() - 1.0e-15);
    EXPECT_LE(birthTimes.back(), 0.5 * sampler.getEmissionTime() + 1.0e-15);
    EXPECT_DOUBLE_EQ(sampler.getGlobalTimeShift(), 0.5 * sampler.getEmissionTime());
    EXPECT_DOUBLE_EQ(sampler.getEmissionTimeStep(), sampler.getEmissionTime() / 100.0);

    EXPECT_EQ(sampler.getNextGlobalIndex(), 0u);
    EXPECT_EQ(globalLocalNum(), 0u);
}

TEST_F(OpalFlatTopTest, EmitsAllInventoryWithVariableDtAndNoRemainderLoss) {
    const Vector_t<double, 3> sigmaR = {0.5, 0.5, 0.0};
    const Vector_t<double, 3> cutoff = {4.0, 4.0, 4.0};
    OpalFlatTop sampler(
            pc, fc,
            /*emitting=*/true,
            /*sigmaTFall=*/0.1,
            /*sigmaTRise=*/0.1, cutoff,
            /*tPulseLengthFWHM=*/1.0, sigmaR);
    sampler.setEmissionOffsets(0.0, Vector_t<double, 3>({0.0, 0.0, 0.1}), 0.0, "NONE");

    const size_t totalParticles = 1007;
    const size_t nranks         = static_cast<size_t>(std::max(1, ippl::Comm->size()));
    pc->allocateParticles(totalParticles / nranks + 2 * nranks + 16);

    size_t mutableTotal = totalParticles;
    sampler.generateParticles(mutableTotal, nr);

    double t           = -sampler.getGlobalTimeShift();
    const double dts[] = {0.07, 0.11, 0.05, 0.13};
    for (int step = 0; step < 1000 && !sampler.isEmissionDone(t); ++step) {
        const double dt = dts[step % 4];
        sampler.emitParticles(t, dt);
        t += dt;
    }

    EXPECT_TRUE(sampler.isEmissionDone(t));
    EXPECT_EQ(sampler.getNextGlobalIndex(), totalParticles);
    EXPECT_EQ(globalLocalNum(), totalParticles);
}

TEST_F(OpalFlatTopTest, RejectsOverdueBirthTimesInsteadOfPreAgingParticles) {
    const Vector_t<double, 3> sigmaR = {0.5, 0.5, 0.0};
    const Vector_t<double, 3> cutoff = {4.0, 4.0, 4.0};
    OpalFlatTop sampler(
            pc, fc,
            /*emitting=*/true,
            /*sigmaTFall=*/1.0e-12,
            /*sigmaTRise=*/1.0e-12, cutoff,
            /*tPulseLengthFWHM=*/4.0e-12, sigmaR);

    const Vector_t<double, 3> P0 = {0.0, 0.0, 0.1};
    sampler.setEmissionOffsets(0.0, P0, 0.0, "NONE");
    sampler.setBirthTimesForTest({-4.0e-13, 1.0e-13, 4.0e-13, 9.0e-13});

    const size_t nranks = static_cast<size_t>(std::max(1, ippl::Comm->size()));
    pc->allocateParticles(4 / nranks + 2 * nranks + 8);

    EXPECT_THROW(sampler.emitParticles(0.0, 5.0e-13), OpalException);
    EXPECT_EQ(sampler.getNextGlobalIndex(), 0u);
    EXPECT_EQ(globalLocalNum(), 0u);
}

TEST_F(OpalFlatTopTest, UsesStoredBirthTimesForParticleDtAndZCorrection) {
    const Vector_t<double, 3> sigmaR = {0.5, 0.5, 0.0};
    const Vector_t<double, 3> cutoff = {4.0, 4.0, 4.0};
    OpalFlatTop sampler(
            pc, fc,
            /*emitting=*/true,
            /*sigmaTFall=*/1.0e-12,
            /*sigmaTRise=*/1.0e-12, cutoff,
            /*tPulseLengthFWHM=*/4.0e-12, sigmaR);

    const Vector_t<double, 3> P0 = {0.0, 0.0, 0.1};
    sampler.setEmissionOffsets(0.0, P0, 0.0, "NONE");
    sampler.setBirthTimesForTest({1.0e-13, 4.0e-13, 9.0e-13});

    const size_t nranks = static_cast<size_t>(std::max(1, ippl::Comm->size()));
    pc->allocateParticles(3 / nranks + 2 * nranks + 8);

    sampler.emitParticles(0.0, 5.0e-13);

    const size_t totalNew       = 2;
    const size_t rank           = static_cast<size_t>(ippl::Comm->rank());
    const size_t base           = totalNew / nranks;
    const size_t rem            = totalNew % nranks;
    const size_t nlocalExpected = base + (rank < rem ? 1 : 0);
    const size_t localOffset    = rank * base + std::min(rank, rem);

    ASSERT_EQ(pc->getLocalNum(), nlocalExpected);

    auto RviewDevice  = pc->R.getView();
    auto PviewDevice  = pc->P.getView();
    auto dtviewDevice = pc->dt.getView();
    auto Rview        = Kokkos::create_mirror_view(RviewDevice);
    auto Pview        = Kokkos::create_mirror_view(PviewDevice);
    auto dtview       = Kokkos::create_mirror_view(dtviewDevice);
    Kokkos::deep_copy(Rview, RviewDevice);
    Kokkos::deep_copy(Pview, PviewDevice);
    Kokkos::deep_copy(dtview, dtviewDevice);

    const std::vector<double> birthTimes = {1.0e-13, 4.0e-13};
    const double gamma                   = std::sqrt(1.0 + P0[2] * P0[2]);
    const double betaZ                   = P0[2] / gamma;

    for (size_t i = 0; i < nlocalExpected; ++i) {
        const double expectedDt = 5.0e-13 - birthTimes[localOffset + i];
        EXPECT_NEAR(dtview(i), expectedDt, 1.0e-18);
        EXPECT_DOUBLE_EQ(Pview(i)[0], P0[0]);
        EXPECT_DOUBLE_EQ(Pview(i)[1], P0[1]);
        EXPECT_DOUBLE_EQ(Pview(i)[2], P0[2]);
        EXPECT_NEAR(Rview(i)[2], 0.5 * betaZ * Physics::c * expectedDt, 1.0e-12);
    }

    EXPECT_EQ(sampler.getNextGlobalIndex(), totalNew);
    EXPECT_EQ(globalLocalNum(), totalNew);
}

TEST_F(OpalFlatTopTest, ProvidesOldOpalInitialReferenceMomentum) {
    const Vector_t<double, 3> sigmaR = {0.5, 0.5, 0.0};
    const Vector_t<double, 3> cutoff = {4.0, 4.0, 4.0};
    OpalFlatTop sampler(
            pc, fc,
            /*emitting=*/true,
            /*sigmaTFall=*/1.0e-12,
            /*sigmaTRise=*/1.0e-12, cutoff,
            /*tPulseLengthFWHM=*/4.0e-12, sigmaR);

    const Vector_t<double, 3> P0 = {0.01, 0.02, 0.1};
    sampler.setEmissionOffsets(0.0, P0, 0.0, "NONE");
    EXPECT_TRUE(sampler.hasInitialReferenceMomentum());
    EXPECT_DOUBLE_EQ(sampler.getInitialReferenceMomentum()[0], P0[0]);
    EXPECT_DOUBLE_EQ(sampler.getInitialReferenceMomentum()[1], P0[1]);
    EXPECT_DOUBLE_EQ(sampler.getInitialReferenceMomentum()[2], P0[2]);

    sampler.setEmissionOffsets(0.0, P0, 0.0, "ASTRA");
    const Vector_t<double, 3> refP = sampler.getInitialReferenceMomentum();
    EXPECT_DOUBLE_EQ(refP[0], 0.0);
    EXPECT_DOUBLE_EQ(refP[1], 0.0);
    EXPECT_DOUBLE_EQ(refP[2], 0.5 * std::sqrt(dot(P0, P0)));
}
