#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>

#include "Distribution/EmittedFromFile.h"
#include "Ippl.h"
#include "PartBunch/BunchStateHandler.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "Utility/IpplTimings.h"

class EmittedFromFileTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() { ippl::finalize(); }

    void SetUp() override {
        nr                             = 8;
        ippl::Vector<double, 3> rmin   = -1.0;
        ippl::Vector<double, 3> rmax   = 1.0;
        ippl::Vector<double, 3> origin = rmin;
        ippl::Vector<double, 3> hr     = (rmax - rmin) / ippl::Vector<double, 3>(nr);
        std::array<bool, 3> decomp     = {false, false, true};

        ippl::NDIndex<3> domain;
        for (unsigned i = 0; i < 3; ++i) {
            domain[i] = ippl::Index(nr[i]);
        }

        auto fc = std::make_shared<FieldContainer_t>(
                hr, rmin, rmax, decomp, domain, origin, isAllPeriodic_m);

        Mesh_t<3> mesh(domain, hr, origin);
        FieldLayout_t<3> fl(MPI_COMM_WORLD, domain, decomp, isAllPeriodic_m);

        pc                = std::make_shared<ParticleContainer<double, 3>>(mesh, fl);
        bunchStateHandler = std::make_shared<BunchStateHandler>();
        pc->setBunchStateHandler(bunchStateHandler);
        tempFilename = "emittedfromfile_test_input.dat";
    }

    void TearDown() override {
        if (!tempFilename.empty()) {
            std::remove(tempFilename.c_str());
        }
    }

    void allocate(size_t globalCapacity) {
        const int nranks           = std::max(1, ippl::Comm->size());
        const size_t nranksU       = static_cast<size_t>(nranks);
        const size_t localCapacity = globalCapacity / nranksU + 2 * nranksU + 1;
        pc->allocateParticles(localCapacity);
    }

    size_t globalParticleCount() const {
        unsigned long local  = static_cast<unsigned long>(pc->getLocalNum());
        unsigned long global = 0;
        MPI_Allreduce(&local, &global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        return static_cast<size_t>(global);
    }

    void writeOldOpalDump() {
        std::ofstream out(tempFilename);
        ASSERT_TRUE(out.is_open());
        out << "# x [m]            px [betax gamma] y [m]            py [betay gamma] "
               "t [s]            pz [betaz gamma] Bin Number\n";
        out << "1.0e-3 1.0e-1 2.0e-3 2.0e-1 -2.0e-12 3.0e-1 1\n";
        out << "3.0e-3 4.0e-1 4.0e-3 5.0e-1 -1.0e-12 6.0e-1 1\n";
    }

    void writeCountHeaderDump() {
        std::ofstream out(tempFilename);
        ASSERT_TRUE(out.is_open());
        out << "3\n";
        out << "x px y py t pz\n";
        out << "1.0e-3 0.0 2.0e-3 0.0 -2.0e-12 2.0e-1\n";
        out << "3.0e-3 0.0 4.0e-3 0.0 -5.0e-13 4.0e-1\n";
        out << "5.0e-3 0.0 6.0e-3 0.0 -1.0e-12 6.0e-1\n";
    }

    std::shared_ptr<ParticleContainer<double, 3>> pc;
    std::shared_ptr<BunchStateHandler> bunchStateHandler;
    ippl::Vector<int, 3> nr;
    bool isAllPeriodic_m = true;
    std::string tempFilename;
};

TEST_F(EmittedFromFileTest, ParsesOldOpalDumpAndProvidesReferenceMomentum) {
    writeOldOpalDump();
    allocate(4);

    auto fc = std::shared_ptr<FieldContainer_t>();
    EmittedFromFile sampler(pc, fc, tempFilename);
    sampler.setEmissionOffsets(
            Vector_t<double, 3>(0.0, 0.0, 0.0), Vector_t<double, 3>(0.01, 0.02, 0.03), 0.0, "NONE");

    size_t requested = 0;
    sampler.generateParticles(requested, nr);

    EXPECT_EQ(requested, static_cast<size_t>(2));
    EXPECT_EQ(globalParticleCount(), static_cast<size_t>(0));
    ASSERT_TRUE(sampler.hasInitialReferenceMomentum());

    const Vector_t<double, 3> refP = sampler.getInitialReferenceMomentum();
    EXPECT_DOUBLE_EQ(refP[0], 0.25 + 0.01);
    EXPECT_DOUBLE_EQ(refP[1], 0.35 + 0.02);
    EXPECT_DOUBLE_EQ(refP[2], 0.45 + 0.03);
    EXPECT_DOUBLE_EQ(sampler.getEmissionTime(), 1.0e-12);
    EXPECT_DOUBLE_EQ(sampler.getGlobalTimeShift(), 0.5e-12);
    EXPECT_DOUBLE_EQ(sampler.getEmissionTimeStep(), 1.0e-14);
}

TEST_F(EmittedFromFileTest, EmitsSortedRecordsWithFractionalDtAndHalfStepDrift) {
    writeCountHeaderDump();
    allocate(4);

    auto fc = std::shared_ptr<FieldContainer_t>();
    EmittedFromFile sampler(pc, fc, tempFilename);
    sampler.setEmissionOffsets(
            Vector_t<double, 3>(0.1, 0.2, 0.3), Vector_t<double, 3>(0.0), 0.0, "NONE");

    size_t requested = 0;
    sampler.generateParticles(requested, nr);
    ASSERT_EQ(requested, static_cast<size_t>(3));

    const double tStart = -sampler.getGlobalTimeShift();
    sampler.emitParticles(tStart, 1.0e-12);
    EXPECT_EQ(globalParticleCount(), static_cast<size_t>(2));

    if (ippl::Comm->size() == 1) {
        auto Rview_d  = pc->R.getView();
        auto Pview_d  = pc->P.getView();
        auto dtView_d = pc->dt.getView();
        auto Rview    = Kokkos::create_mirror_view(Rview_d);
        auto Pview    = Kokkos::create_mirror_view(Pview_d);
        auto dtView   = Kokkos::create_mirror_view(dtView_d);
        Kokkos::deep_copy(Rview, Rview_d);
        Kokkos::deep_copy(Pview, Pview_d);
        Kokkos::deep_copy(dtView, dtView_d);

        const double firstDt     = 1.0e-12;
        const double secondDt    = 0.5e-12;
        const double firstGamma  = std::sqrt(1.0 + 0.4 * 0.4);
        const double secondGamma = std::sqrt(1.0 + 0.6 * 0.6);

        EXPECT_NEAR(dtView(0), firstDt, 1.0e-18);
        EXPECT_NEAR(dtView(1), secondDt, 1.0e-18);

        EXPECT_DOUBLE_EQ(Pview(0)[2], 0.4);
        EXPECT_DOUBLE_EQ(Pview(1)[2], 0.6);
        EXPECT_NEAR(Rview(0)[0], 0.103, 1.0e-15);
        EXPECT_NEAR(Rview(0)[1], 0.204, 1.0e-15);
        EXPECT_NEAR(Rview(0)[2], 0.3 + 0.5 * Physics::c * firstDt * 0.4 / firstGamma, 1.0e-15);
        EXPECT_NEAR(Rview(1)[2], 0.3 + 0.5 * Physics::c * secondDt * 0.6 / secondGamma, 1.0e-15);
    }

    sampler.emitParticles(tStart + 1.0e-12, 1.0e-12);
    EXPECT_EQ(globalParticleCount(), static_cast<size_t>(3));
    EXPECT_TRUE(sampler.isEmissionDone(tStart + 2.0e-12));

    sampler.emitParticles(tStart + 2.0e-12, 1.0e-12);
    EXPECT_EQ(globalParticleCount(), static_cast<size_t>(3));
}

TEST_F(EmittedFromFileTest, HonorsRequestedParticleLimitBeforeSorting) {
    writeCountHeaderDump();
    allocate(4);

    auto fc = std::shared_ptr<FieldContainer_t>();
    EmittedFromFile sampler(pc, fc, tempFilename);
    sampler.setEmissionOffsets(Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0), 0.0, "NONE");

    size_t requested = 2;
    sampler.generateParticles(requested, nr);
    EXPECT_EQ(requested, static_cast<size_t>(2));

    const Vector_t<double, 3> refP = sampler.getInitialReferenceMomentum();
    EXPECT_DOUBLE_EQ(refP[0], 0.0);
    EXPECT_DOUBLE_EQ(refP[1], 0.0);
    EXPECT_DOUBLE_EQ(refP[2], 0.3);

    sampler.emitParticles(-sampler.getGlobalTimeShift(), 3.0e-12);
    EXPECT_EQ(globalParticleCount(), static_cast<size_t>(2));
}

TEST_F(EmittedFromFileTest, RejectsOverdueBirthTimesInsteadOfPreAgingParticles) {
    writeCountHeaderDump();
    allocate(4);

    auto fc = std::shared_ptr<FieldContainer_t>();
    EmittedFromFile sampler(pc, fc, tempFilename);
    sampler.setEmissionOffsets(Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0), 0.0, "NONE");

    size_t requested = 0;
    sampler.generateParticles(requested, nr);

    EXPECT_THROW(sampler.emitParticles(1.0e-12, 1.0e-12), OpalException);
    EXPECT_EQ(globalParticleCount(), static_cast<size_t>(0));
}

TEST_F(EmittedFromFileTest, RejectsUnsupportedEmissionModel) {
    writeOldOpalDump();
    allocate(4);

    auto fc = std::shared_ptr<FieldContainer_t>();
    EmittedFromFile sampler(pc, fc, tempFilename);
    sampler.setEmissionOffsets(Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0), 0.0, "ASTRA");

    size_t requested = 0;
    EXPECT_THROW(sampler.generateParticles(requested, nr), OpalException);
}
