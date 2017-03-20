#include "gtest/gtest.h"

#include "mpi.h"

#include "Utility/Inform.h" // ippl

extern Inform* gmsg;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    gmsg = new Inform("UnitTests: ", std::cerr);
    if (!gmsg) {
        return 1;
    }
    MPI_Init(&argc, &argv);
    int test_out = RUN_ALL_TESTS();
    MPI_Finalize();

    return test_out;
}
