#include <gtest/gtest.h>
#include <mpi.h>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <memory>

#include "Distribution/Distribution.h"
#include "Distribution/FromFile.h"
#include "Ippl.h"
#include "PartBunch/BunchStateHandler.h"
#include "Utility/IpplTimings.h"

class FromFileTest : public ::testing::Test {
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
        for (unsigned i = 0; i < 3; i++) {
            domain[i] = ippl::Index(this->nr[i]);
        }

        auto fc = std::make_shared<FieldContainer_t>(
                hr, rmin, rmax, decomp, domain, origin, this->isAllPeriodic_m);

        Mesh_t<3> mesh(domain, hr, origin);
        FieldLayout_t<3> fl(MPI_COMM_WORLD, domain, decomp, this->isAllPeriodic_m);

        pc = std::make_shared<ParticleContainer<double, 3>>(mesh, fl);

        bunchStateHandler = std::make_shared<BunchStateHandler>();
        pc->setBunchStateHandler(bunchStateHandler);
        tempFilename = "fromfile_test_input.dat";
    }

    void TearDown() override {
        if (!tempFilename.empty()) {
            std::remove(tempFilename.c_str());
        }
    }

    void writeSampleFile() {
        // Strict format: N, then header, then data (comments/blanks skipped in data section)
        std::ofstream out(tempFilename);
        ASSERT_TRUE(out.is_open());

        out << "5\n";
        out << "x y z px py pz\n";

        // Comment and blank lines in data section are skipped
        out << "# This is a comment and must be skipped\n";
        out << "\n";

        // Five data lines with 6 columns each: x y z px py pz
        out << "7.3797510200e-04 2.8811575176e-02 3.3682581175e-04 9.7591095318e-03 "
               "3.3682581175e-04 -2.2083905829e-02\n";
        out << "1.5399764838e-04 1.0245287933e-02 2.7921666726e-04 1.4009747125e-02 "
               "2.7921666726e-04 4.9442454400e-02\n";
        out << "6.1570316063e-05 6.5451409548e-04 2.0566658643e-03 4.6015050865e-02 "
               "2.0566658643e-03 2.4489916771e-02\n";
        out << "-8.5549718574e-04 -3.0244620065e-02 2.5710406506e-04 1.1943422591e-02 "
               "2.5710406506e-04 -3.1888283063e-02\n";
        out << "-4.0851931221e-04 -1.3879415332e-02 2.2037752464e-03 4.9513250233e-02 "
               "2.2037752464e-03 1.4595333192e-02\n";

        out.close();
    }

    void writeFileWithHeader(
            const std::string& headerLine, const std::string& dataLine, size_t numParticles = 1) {
        std::ofstream out(tempFilename);
        ASSERT_TRUE(out.is_open());
        out << numParticles << "\n";
        out << headerLine << "\n";
        out << dataLine << "\n";
        out.close();
    }

    std::shared_ptr<ParticleContainer<double, 3>> pc;
    std::shared_ptr<BunchStateHandler> bunchStateHandler;
    ippl::Vector<int, 3> nr;
    bool isAllPeriodic_m = true;
    std::string tempFilename;
};

TEST_F(FromFileTest, GeneratesParticlesFromAsciiFile) {
    writeSampleFile();

    auto fc = std::shared_ptr<FieldContainer_t>();
    FromFile sampler(pc, fc, tempFilename);

    size_t requested = 10;  // larger than available in file to test clamping

    // Preallocate capacity so SamplingBase::computeLocalEmitCount
    // can distribute up to the requested number of particles.
    const int nranks         = std::max(1, ippl::Comm->size());
    const size_t nranksU     = static_cast<size_t>(nranks);
    const size_t maxLocalNum = requested / nranksU + 2 * nranksU + 1;
    pc->allocateParticles(maxLocalNum);

    sampler.generateParticles(requested, nr);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Sum local particles over all ranks
    size_t localN  = pc->getLocalNum();
    size_t globalN = 0;
    MPI_Allreduce(&localN, &globalN, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

    EXPECT_EQ(globalN, 5u);

    // Copy views to host space for GPU compatibility
    auto Rview_d = pc->R.getView();
    auto Pview_d = pc->P.getView();

    auto Rview = Kokkos::create_mirror_view(Rview_d);
    auto Pview = Kokkos::create_mirror_view(Pview_d);
    Kokkos::deep_copy(Rview, Rview_d);
    Kokkos::deep_copy(Pview, Pview_d);

    if (localN > 0) {
        // Just check that first local particle has finite values
        EXPECT_TRUE(std::isfinite(Rview(0)[0]));
        EXPECT_TRUE(std::isfinite(Rview(0)[1]));
        EXPECT_TRUE(std::isfinite(Rview(0)[2]));
        EXPECT_TRUE(std::isfinite(Pview(0)[0]));
        EXPECT_TRUE(std::isfinite(Pview(0)[1]));
        EXPECT_TRUE(std::isfinite(Pview(0)[2]));
    }
}

TEST_F(FromFileTest, ParseHeader_ReorderedColumns) {
    // Test parseHeader: header "py px pz y x z" => col0=py, col1=px, col2=pz, col3=y, col4=x,
    // col5=z. Data "10 20 30 40 50 60" => x=50, y=40, z=60, px=20, py=10, pz=30.
    writeFileWithHeader("py px pz y x z", "10 20 30 40 50 60");

    auto fc = std::shared_ptr<FieldContainer_t>();
    FromFile sampler(pc, fc, tempFilename);

    size_t requested = 1;
    sampler.generateParticles(requested, nr);

    size_t localN = pc->getLocalNum();
    if (localN > 0) {
        auto Rview_d = pc->R.getView();
        auto Pview_d = pc->P.getView();
        auto Rview   = Kokkos::create_mirror_view(Rview_d);
        auto Pview   = Kokkos::create_mirror_view(Pview_d);
        Kokkos::deep_copy(Rview, Rview_d);
        Kokkos::deep_copy(Pview, Pview_d);

        EXPECT_DOUBLE_EQ(Rview(0)[0], 50.0);  // x
        EXPECT_DOUBLE_EQ(Rview(0)[1], 40.0);  // y
        EXPECT_DOUBLE_EQ(Rview(0)[2], 60.0);  // z
        EXPECT_DOUBLE_EQ(Pview(0)[0], 20.0);  // px
        EXPECT_DOUBLE_EQ(Pview(0)[1], 10.0);  // py
        EXPECT_DOUBLE_EQ(Pview(0)[2], 30.0);  // pz
    }
}

TEST_F(FromFileTest, ParseHeader_StandardOrderAndAlternateNames) {
    writeFileWithHeader("x y z px/momentumx py/momentumy pz/momentumz", "1.0 2.0 3.0 4.0 5.0 6.0");

    auto fc = std::shared_ptr<FieldContainer_t>();
    FromFile sampler(pc, fc, tempFilename);

    size_t requested = 1;
    sampler.generateParticles(requested, nr);

    size_t localN = pc->getLocalNum();
    if (localN > 0) {
        auto Rview_d = pc->R.getView();
        auto Pview_d = pc->P.getView();
        auto Rview   = Kokkos::create_mirror_view(Rview_d);
        auto Pview   = Kokkos::create_mirror_view(Pview_d);
        Kokkos::deep_copy(Rview, Rview_d);
        Kokkos::deep_copy(Pview, Pview_d);

        EXPECT_DOUBLE_EQ(Rview(0)[0], 1.0);
        EXPECT_DOUBLE_EQ(Rview(0)[1], 2.0);
        EXPECT_DOUBLE_EQ(Rview(0)[2], 3.0);
        EXPECT_DOUBLE_EQ(Pview(0)[0], 4.0);
        EXPECT_DOUBLE_EQ(Pview(0)[1], 5.0);
        EXPECT_DOUBLE_EQ(Pview(0)[2], 6.0);
    }
}
