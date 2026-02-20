#include <gtest/gtest.h>
#include <mpi.h>
#include <memory>
#include <fstream>
#include <cstdio>

#include "Distribution/FromFile.h"
#include "Distribution/Distribution.h"
#include "Ippl.h"
#include "Utility/IpplTimings.h"

class FromFileTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc = 0;
        char **argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() {
        ippl::finalize();
    }

    void SetUp() override {
        nr = 8;
        ippl::Vector<double, 3> rmin = -1.0;
        ippl::Vector<double, 3> rmax = 1.0;
        ippl::Vector<double, 3> origin = rmin;
        ippl::Vector<double, 3> hr = (rmax - rmin) / ippl::Vector<double, 3>(nr);
        std::array<bool, 3> decomp = {false, false, true};

        ippl::NDIndex<3> domain;
        for (unsigned i = 0; i < 3; i++) {
            domain[i] = ippl::Index(this->nr[i]);
        }

        auto fc = std::make_shared<FieldContainer_t>(hr, rmin, rmax, decomp, domain, origin, this->isAllPeriodic_m);

        Mesh_t<3> mesh(domain, hr, origin);
        FieldLayout_t<3> fl(MPI_COMM_WORLD, domain, decomp, this->isAllPeriodic_m);

        pc = std::make_shared<ParticleContainer<double, 3>>(mesh, fl);

        tempFilename = "fromfile_test_input.dat";
    }

    void TearDown() override {
        if (!tempFilename.empty()) {
            std::remove(tempFilename.c_str());
        }
    }

    void writeSampleFile() {
        // Small sample copied from tests/Drift-1-fromfile/inputdistr.dat
        std::ofstream out(tempFilename);
        ASSERT_TRUE(out.is_open());

        // First line: number of particles
        out << "5\n";

        // Comment and blank lines that should be ignored
        out << "# This is a comment and must be skipped\n";
        out << "\n";

        // Five data lines with 6 columns each: x y z px py pz
        out << "7.3797510200e-04 2.8811575176e-02 3.3682581175e-04 9.7591095318e-03 3.3682581175e-04 -2.2083905829e-02\n";
        out << "1.5399764838e-04 1.0245287933e-02 2.7921666726e-04 1.4009747125e-02 2.7921666726e-04 4.9442454400e-02\n";
        out << "6.1570316063e-05 6.5451409548e-04 2.0566658643e-03 4.6015050865e-02 2.0566658643e-03 2.4489916771e-02\n";
        out << "-8.5549718574e-04 -3.0244620065e-02 2.5710406506e-04 1.1943422591e-02 2.5710406506e-04 -3.1888283063e-02\n";
        out << "-4.0851931221e-04 -1.3879415332e-02 2.2037752464e-03 4.9513250233e-02 2.2037752464e-03 1.4595333192e-02\n";

        out.close();
    }

    std::shared_ptr<ParticleContainer<double, 3>> pc;
    ippl::Vector<int, 3> nr;
    bool isAllPeriodic_m = true;
    std::string tempFilename;
};

TEST_F(FromFileTest, GeneratesParticlesFromAsciiFile) {
    writeSampleFile();

    auto fc = std::shared_ptr<FieldContainer_t>();
    FromFile sampler(pc, fc, tempFilename);

    size_t requested = 10; // larger than available in file to test clamping
    sampler.generateParticles(requested, nr);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Sum local particles over all ranks
    size_t localN = pc->getLocalNum();
    size_t globalN = 0;
    MPI_Allreduce(&localN, &globalN, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

    EXPECT_EQ(globalN, 5u);

    auto Rview = pc->R.getView();
    auto Pview = pc->P.getView();

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

