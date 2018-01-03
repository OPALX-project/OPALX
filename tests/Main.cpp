#include "gtest/gtest.h"

#include "mpi.h"

#include "Utility/IpplInfo.h" // ippl

Ippl *ippl;
Inform* gmsg;

class NewLineAdder: public ::testing::EmptyTestEventListener {
    virtual void OnTestPartResult(const ::testing::TestPartResult &test_part_result) {
        if (test_part_result.failed())
            printf("\n");
    }
};

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    gmsg = new Inform("UnitTests: ", std::cerr);
    if (!gmsg) {
        return 1;
    }
    MPI_Init(&argc, &argv);

    ::testing::TestEventListeners &listeners =
          ::testing::UnitTest::GetInstance()->listeners();
    listeners.Append(new NewLineAdder);

    int test_out = RUN_ALL_TESTS();
    MPI_Finalize();

    return test_out;
}