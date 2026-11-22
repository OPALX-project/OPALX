#include <gtest/gtest.h>
#include "../TemporaryDirectory.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
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
        pc           = makeContainer(false);
        tempFilename = (directory_m.create("opalx-emitted-fromfile") / "input.dat").string();
    }

    /// Container on the test mesh; @p spin turns on the per-particle polarization Pol.
    std::shared_ptr<ParticleContainer<double, 3>> makeContainer(bool spin) {
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

        // Particle layouts borrow these objects; migration tests require them
        // to outlive the container rather than referring to expired stack locals.
        meshes.emplace_back(std::make_unique<Mesh_t<3>>(domain, hr, origin));
        layouts.emplace_back(
                std::make_unique<FieldLayout_t<3>>(
                        MPI_COMM_WORLD, domain, decomp, isAllPeriodic_m));
        auto& mesh = *meshes.back();
        auto& fl   = *layouts.back();

        auto container    = std::make_shared<ParticleContainer<double, 3>>(mesh, fl, spin);
        bunchStateHandler = std::make_shared<BunchStateHandler>();
        container->setBunchStateHandler(bunchStateHandler);
        return container;
    }

    void TearDown() override { directory_m.reset(); }

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

    void writeMultiBinOldOpalDump() {
        std::ofstream out(tempFilename);
        ASSERT_TRUE(out.is_open());
        out << "# x px y py t pz Bin Number\n";
        out << "0.0 0.0 0.0 0.0 -0.5e-12 0.1 1\n";
        out << "0.0 0.0 0.0 0.0 -1.0e-12 0.2 2\n";
        out << "0.0 0.0 0.0 0.0 -3.0e-12 0.3 3\n";
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

    void writeExplicitBirthDump() {
        std::ofstream out(tempFilename);
        ASSERT_TRUE(out.is_open());
        out << "2\n";
        out << "x y z px py pz birth_time\n";
        out << "1.0e-3 2.0e-3 3.0e-3 0.0 0.0 4.0e-1 -2.5e-13\n";
        out << "5.0e-3 6.0e-3 -4.0e-3 0.0 0.0 6.0e-1 2.5e-13\n";
    }

    std::vector<std::unique_ptr<Mesh_t<3>>> meshes;
    std::vector<std::unique_ptr<FieldLayout_t<3>>> layouts;
    std::shared_ptr<ParticleContainer<double, 3>> pc;
    std::shared_ptr<BunchStateHandler> bunchStateHandler;
    ippl::Vector<int, 3> nr;
    bool isAllPeriodic_m = true;
    opalx::test::TemporaryDirectory directory_m;
    std::string tempFilename;
};

TEST_F(EmittedFromFileTest, ParsesOldOpalDumpAndProvidesReferenceMomentum) {
    writeOldOpalDump();
    allocate(4);

    EmittedFromFile sampler(pc, tempFilename);
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
    EXPECT_DOUBLE_EQ(sampler.getEmissionTime(), 2.0e-12);
    EXPECT_DOUBLE_EQ(sampler.getGlobalTimeShift(), 1.0e-12);
    EXPECT_DOUBLE_EQ(sampler.getEmissionTimeStep(), 2.0e-14);
}

TEST_F(EmittedFromFileTest, ReconstructsOldOpalEmissionWindowFromBins) {
    writeMultiBinOldOpalDump();
    allocate(4);

    EmittedFromFile sampler(pc, tempFilename);
    sampler.setEmissionOffsets(Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0), 0.0, "NONE");

    size_t requested = 0;
    sampler.generateParticles(requested, nr);

    EXPECT_EQ(requested, static_cast<size_t>(3));
    EXPECT_DOUBLE_EQ(sampler.getEmissionTime(), 3.0e-12);
    EXPECT_DOUBLE_EQ(sampler.getGlobalTimeShift(), 1.5e-12);
    EXPECT_DOUBLE_EQ(sampler.getEmissionTimeStep(), 3.0e-14);
}

// The external AWA regressions use a count followed by an uncommented header
// with repeated unit labels and the two-word name "Bin Number". Test all
// count/comment combinations, preserving the legacy bin-derived time window.
TEST_F(EmittedFromFileTest, AnnotatedLegacyHeadersPreserveBinsAndPhaseSpace) {
    for (bool count : {false, true}) {
        for (bool comment : {false, true}) {
            SCOPED_TRACE(::testing::Message() << "count=" << count << " comment=" << comment);
            {
                std::ofstream out(tempFilename);
                if (count) out << "3\n";
                if (comment) out << "# ";
                out << "X [m] PX [betax gamma] Y [m] PY [betay gamma] "
                       "T [s] PZ [betaz gamma] Bin Number\n";
                out << "0.001 0 0.002 0 -0.5e-12 0.1 1\n";
                out << "0.003 0 0.004 0 -1.0e-12 0.2 2\n";
                out << "0.005 0 0.006 0 -3.0e-12 0.3 3\n";
            }
            pc = makeContainer(false);
            allocate(4);
            EmittedFromFile sampler(pc, tempFilename);
            sampler.setEmissionOffsets(
                    Vector_t<double, 3>(0.1, 0.2, 0.3), Vector_t<double, 3>(0.0), 0.0, "NONE");
            size_t requested = 0;
            sampler.generateParticles(requested, nr);
            EXPECT_EQ(requested, 3u);
            EXPECT_DOUBLE_EQ(sampler.getEmissionTime(), 3e-12);
            EXPECT_DOUBLE_EQ(sampler.getGlobalTimeShift(), 1.5e-12);
            EXPECT_NEAR(sampler.getInitialReferenceMomentum()[2], 0.2, 1e-16);
            EXPECT_DOUBLE_EQ(sampler.getInitialReferencePosition()[2], 0.3);
            sampler.emitParticles(-1.5e-12, 3.25e-12);
            EXPECT_EQ(globalParticleCount(), 3u);
            EXPECT_TRUE(sampler.isEmissionDone(1.75e-12));
            const auto r =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
            const auto p =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
            const auto dt =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->dt.getView());
            const std::array<double, 3> fractional{2.75e-12, 2.25e-12, 0.25e-12};
            std::array<int, 3> localSeen{}, globalSeen{};
            for (size_t i = 0; i < pc->getLocalNum(); ++i) {
                const int row = static_cast<int>(std::lround((r(i)[0] - 0.101) / 0.002));
                ASSERT_GE(row, 0);
                ASSERT_LT(row, 3);
                ++localSeen[row];
                const double pz = 0.1 * (row + 1);
                EXPECT_NEAR(r(i)[0], 0.101 + 0.002 * row, 1e-15);
                EXPECT_NEAR(r(i)[1], 0.202 + 0.002 * row, 1e-15);
                EXPECT_DOUBLE_EQ(p(i)[0], 0.0);
                EXPECT_DOUBLE_EQ(p(i)[1], 0.0);
                EXPECT_NEAR(p(i)[2], pz, 1e-16);
                EXPECT_NEAR(dt(i), fractional[row], 1e-27);
                EXPECT_NEAR(
                        r(i)[2],
                        0.3 + 0.5 * Physics::c * fractional[row] * pz / std::sqrt(1.0 + pz * pz),
                        1e-15);
            }
            MPI_Allreduce(localSeen.data(), globalSeen.data(), 3, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            for (int seen : globalSeen)
                EXPECT_EQ(seen, 1);
        }
    }
}

TEST_F(EmittedFromFileTest, AnnotatedExplicitHeadersPreserveNamedColumnsAndBirthTimes) {
    for (bool count : {false, true}) {
        for (bool comment : {false, true}) {
            SCOPED_TRACE(::testing::Message() << "count=" << count << " comment=" << comment);
            {
                std::ofstream out(tempFilename);
                if (count) out << "1\n";
                if (comment) out << "# ";
                out << "PZ [beta gamma] T [s] Z [m] PY [beta gamma] "
                       "X [m] PX [beta gamma] Y [m]\n";
                out << "0.6 2.5e-13 0.003 0.2 0.001 0.1 0.002\n";
            }
            pc = makeContainer(false);
            allocate(2);
            EmittedFromFile sampler(pc, tempFilename);
            sampler.setEmissionOffsets(
                    Vector_t<double, 3>(0.1, 0.2, 0.3), Vector_t<double, 3>(0.0), 0.5e-12, "NONE");
            size_t requested = 0;
            sampler.generateParticles(requested, nr);
            EXPECT_EQ(requested, 1u);
            EXPECT_DOUBLE_EQ(sampler.getGlobalTimeShift(), 0.0);
            EXPECT_DOUBLE_EQ(sampler.getInitialReferencePosition()[2], 0.303);
            sampler.emitParticles(0.0, 0.5e-12);
            EXPECT_EQ(globalParticleCount(), 0u);
            sampler.emitParticles(0.5e-12, 0.5e-12);
            EXPECT_EQ(globalParticleCount(), 1u);
            const auto r =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
            const auto p =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
            const auto dt =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->dt.getView());
            const Vector_t<double, 3> initial(0.101, 0.202, 0.303), momentum(0.1, 0.2, 0.6);
            for (size_t i = 0; i < pc->getLocalNum(); ++i) {
                EXPECT_NEAR(dt(i), 0.25e-12, 1e-27);
                for (unsigned d = 0; d < 3; ++d) {
                    EXPECT_DOUBLE_EQ(p(i)[d], momentum[d]);
                    EXPECT_NEAR(
                            r(i)[d],
                            initial[d]
                                    + 0.5 * Physics::c * 0.25e-12 * momentum[d] / std::sqrt(1.41),
                            1e-15);
                }
            }
        }
    }
}

TEST_F(EmittedFromFileTest, RejectsDuplicateColumnsAndMalformedAnnotatedHeaders) {
    const std::array<const char*, 5> headers{
            "x [m] y [m] z [m] px py pz birth_time x", "x y z px py pz birth_time X",
            "x [m y z px py pz birth_time", "x m] y z px py pz birth_time", "x y z px py pz"};
    for (const char* header : headers) {
        SCOPED_TRACE(header);
        {
            std::ofstream out(tempFilename);
            out << "1\n" << header << "\n0 0 0 0 0 0 0\n";
        }
        EXPECT_THROW(EmittedFromFile sampler(pc, tempFilename), OpalException);
    }
}

TEST_F(EmittedFromFileTest, RejectsDuplicateCommentHeaderInsteadOfFallingBackToLegacy) {
    {
        std::ofstream out(tempFilename);
        out << "# x [m] y [m] z [m] px py pz birth_time X\n";
        // A positive integer last value would look like a valid legacy bin
        // if the malformed explicit comment were ignored instead of rejected.
        out << "0 0 0 0 0 0 1\n";
    }
    EXPECT_THROW(EmittedFromFile sampler(pc, tempFilename), OpalException);
}

TEST_F(EmittedFromFileTest, EmitsSortedRecordsWithFractionalDtAndHalfStepDrift) {
    writeCountHeaderDump();
    allocate(4);

    EmittedFromFile sampler(pc, tempFilename);
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
        EXPECT_NEAR(Rview(0)[0], 0.103, 1.0e-12);
        EXPECT_NEAR(Rview(0)[1], 0.204, 1.0e-12);
        EXPECT_NEAR(Rview(0)[2], 0.3 + 0.5 * Physics::c * firstDt * 0.4 / firstGamma, 1.0e-12);
        EXPECT_NEAR(Rview(1)[2], 0.3 + 0.5 * Physics::c * secondDt * 0.6 / secondGamma, 1.0e-12);
    }

    sampler.emitParticles(tStart + 1.0e-12, 1.0e-12);
    EXPECT_EQ(globalParticleCount(), static_cast<size_t>(3));
    EXPECT_TRUE(sampler.isEmissionDone(tStart + 2.0e-12));

    sampler.emitParticles(tStart + 2.0e-12, 1.0e-12);
    EXPECT_EQ(globalParticleCount(), static_cast<size_t>(3));
}

TEST_F(EmittedFromFileTest, PreservesExplicitBirthSpacetimeAndFractionalFirstStep) {
    writeExplicitBirthDump();
    allocate(4);

    EmittedFromFile sampler(pc, tempFilename);
    sampler.setEmissionOffsets(
            Vector_t<double, 3>(0.1, 0.2, 0.3), Vector_t<double, 3>(0.0), 1.0e-12, "NONE");

    size_t requested = 0;
    sampler.generateParticles(requested, nr);
    ASSERT_EQ(requested, static_cast<size_t>(2));
    EXPECT_DOUBLE_EQ(sampler.getEmissionTime(), 5.0e-13);
    EXPECT_DOUBLE_EQ(sampler.getGlobalTimeShift(), 0.0);

    const Vector_t<double, 3> refR = sampler.getInitialReferencePosition();
    EXPECT_NEAR(refR[0], 0.103, 1.0e-15);
    EXPECT_NEAR(refR[1], 0.204, 1.0e-15);
    EXPECT_NEAR(refR[2], 0.2995, 1.0e-15);

    sampler.emitParticles(0.5e-12, 0.6e-12);
    EXPECT_EQ(globalParticleCount(), static_cast<size_t>(1));

    if (ippl::Comm->size() == 1) {
        auto Rview_d  = pc->R.getView();
        auto dtView_d = pc->dt.getView();
        auto Rview    = Kokkos::create_mirror_view(Rview_d);
        auto dtView   = Kokkos::create_mirror_view(dtView_d);
        Kokkos::deep_copy(Rview, Rview_d);
        Kokkos::deep_copy(dtView, dtView_d);

        const double fractionalDt = 0.35e-12;
        const double gamma        = std::sqrt(1.0 + 0.4 * 0.4);
        EXPECT_NEAR(dtView(0), fractionalDt, 1.0e-18);
        EXPECT_NEAR(Rview(0)[0], 0.101, 1.0e-15);
        EXPECT_NEAR(Rview(0)[1], 0.202, 1.0e-15);
        EXPECT_NEAR(Rview(0)[2], 0.303 + 0.5 * Physics::c * fractionalDt * 0.4 / gamma, 1.0e-15);
    }

    sampler.emitParticles(1.1e-12, 0.2e-12);
    EXPECT_EQ(globalParticleCount(), static_cast<size_t>(2));
    EXPECT_TRUE(sampler.isEmissionDone(1.3e-12));
}

// Each identity must be born exactly once, including births exactly at a step
// boundary (zero fractional dt). Validate on every rank, not only rank one.
TEST_F(EmittedFromFileTest, ReorderedCommentHeaderSelectsExplicitTimesAndCoordinates) {
    {
        std::ofstream out(tempFilename);
        out << "1\n# PZ t Z PY X PX Y\n";
        out << "0.6 2e-13 0.003 0.2 0.001 0.1 0.002\n";
    }
    allocate(2);
    EmittedFromFile sampler(pc, tempFilename);
    sampler.setEmissionOffsets(Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0), 0.0, "NONE");
    size_t requested = 0;
    sampler.generateParticles(requested, nr);
    ASSERT_EQ(requested, 1u);
    const auto r = sampler.getInitialReferencePosition();
    const auto p = sampler.getInitialReferenceMomentum();
    EXPECT_DOUBLE_EQ(r[0], 0.001);
    EXPECT_DOUBLE_EQ(r[1], 0.002);
    EXPECT_DOUBLE_EQ(r[2], 0.003);
    EXPECT_DOUBLE_EQ(p[0], 0.1);
    EXPECT_DOUBLE_EQ(p[1], 0.2);
    EXPECT_DOUBLE_EQ(p[2], 0.6);
    sampler.emitParticles(0.0, 1e-13);
    EXPECT_EQ(globalParticleCount(), 0u);
    sampler.emitParticles(1e-13, 1e-13);
    EXPECT_EQ(globalParticleCount(), 1u);
}

TEST_F(EmittedFromFileTest, TwelveTimedBirthsPreserveIdentityAndBoundarySteps) {
    constexpr double step = 1e-12;
    const std::array<double, 12> births{0.0,  0.25, 1.0, 1.25, 1.5, 2.0,
                                        2.25, 2.5,  3.0, 3.25, 3.5, 4.0};
    {
        std::ofstream out(tempFilename);
        ASSERT_TRUE(out.is_open());
        out << "12\nx y z px py pz birth_time\n" << std::setprecision(17);
        for (std::size_t i = 0; i < births.size(); ++i)
            out << (i + 1) * 1e-3 << " 0 0.002 0 0 0.4 " << births[i] * step << '\n';
    }
    allocate(12);
    EmittedFromFile sampler(pc, tempFilename);
    sampler.setEmissionOffsets(Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0), 0.0, "NONE");
    size_t requested = 12;
    sampler.generateParticles(requested, nr);
    ASSERT_EQ(requested, 12u);
    std::array<int, 12> localIdentities{}, globalIdentities{};
    std::array<long long, 12> acceptedIds;
    acceptedIds.fill(-1);
    for (unsigned interval = 0; interval < 4; ++interval) {
        const auto before = pc->getLocalNum();
        sampler.emitParticles(interval * step, step);
        EXPECT_EQ(globalParticleCount(), 3u * (interval + 1));
        const auto r  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
        const auto p  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
        const auto dt = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->dt.getView());
        for (std::size_t i = before; i < pc->getLocalNum(); ++i) {
            const int identity = static_cast<int>(std::lround(r(i)[0] / 1e-3)) - 1;
            EXPECT_GE(identity, 0);
            EXPECT_LT(identity, 12);
            if (identity < 0 || identity >= 12) continue;
            ++localIdentities[identity];
            const double fractional = (interval + 1) * step - births[identity] * step;
            EXPECT_NEAR(dt(i), fractional, 1e-27);
            EXPECT_DOUBLE_EQ(p(i)[2], 0.4);
            EXPECT_NEAR(
                    r(i)[2], 0.002 + 0.5 * Physics::c * fractional * 0.4 / std::sqrt(1.16), 1e-15);
        }
        // IDs need not match between MPI sizes, but each file row must retain
        // its birth ID across later emission batches and spatial migration.
        pc->update();
        const auto movedR =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
        const auto ids = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->ID.getView());
        std::array<long long, 12> localIds, globalIds;
        localIds.fill(-1);
        for (size_t i = 0; i < pc->getLocalNum(); ++i) {
            const int row = static_cast<int>(std::lround(movedR(i)[0] / 1e-3)) - 1;
            ASSERT_GE(row, 0);
            ASSERT_LT(row, 12);
            localIds[row] = ids(i);
        }
        MPI_Allreduce(
                localIds.data(), globalIds.data(), 12, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        for (size_t row = 0; row < 12; ++row) {
            if (acceptedIds[row] >= 0) {
                EXPECT_EQ(globalIds[row], acceptedIds[row]);
            }
            if (globalIds[row] >= 0) acceptedIds[row] = globalIds[row];
        }
    }
    EXPECT_TRUE(sampler.isEmissionDone(4 * step));
    sampler.emitParticles(4 * step, step);
    EXPECT_EQ(globalParticleCount(), 12u);
    MPI_Allreduce(
            localIdentities.data(), globalIdentities.data(), 12, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    for (const auto count : globalIdentities)
        EXPECT_EQ(count, 1);
}

TEST_F(EmittedFromFileTest, EmittedParticlesGetBeamPolarization) {
    pc = makeContainer(/*spin=*/true);
    writeCountHeaderDump();
    allocate(4);

    EmittedFromFile sampler(pc, tempFilename);
    sampler.setEmissionOffsets(Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0), 0.0, "NONE");
    sampler.setInitialPolarization({0.0, 0.0, 1.0});

    size_t requested = 0;
    sampler.generateParticles(requested, nr);
    const double tStart = -sampler.getGlobalTimeShift();
    sampler.emitParticles(tStart, 1.0e-12);
    sampler.emitParticles(tStart + 1.0e-12, 1.0e-12);
    ASSERT_EQ(globalParticleCount(), static_cast<size_t>(3));

    auto pol = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->Pol.getView());
    for (size_t i = 0; i < pc->getLocalNum(); ++i) {
        EXPECT_DOUBLE_EQ(pol(i)[0], 0.0) << "particle " << i;
        EXPECT_DOUBLE_EQ(pol(i)[1], 0.0) << "particle " << i;
        EXPECT_DOUBLE_EQ(pol(i)[2], 1.0) << "particle " << i;
    }
}

TEST_F(EmittedFromFileTest, HonorsRequestedParticleLimitBeforeSorting) {
    writeCountHeaderDump();
    allocate(4);

    EmittedFromFile sampler(pc, tempFilename);
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

    EmittedFromFile sampler(pc, tempFilename);
    sampler.setEmissionOffsets(Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0), 0.0, "NONE");

    size_t requested = 0;
    sampler.generateParticles(requested, nr);

    EXPECT_THROW(sampler.emitParticles(1.0e-12, 1.0e-12), OpalException);
    EXPECT_EQ(globalParticleCount(), static_cast<size_t>(0));
}

TEST_F(EmittedFromFileTest, RejectsUnsupportedEmissionModel) {
    writeOldOpalDump();
    allocate(4);

    EmittedFromFile sampler(pc, tempFilename);
    sampler.setEmissionOffsets(Vector_t<double, 3>(0.0), Vector_t<double, 3>(0.0), 0.0, "ASTRA");

    size_t requested = 0;
    EXPECT_THROW(sampler.generateParticles(requested, nr), OpalException);
}
