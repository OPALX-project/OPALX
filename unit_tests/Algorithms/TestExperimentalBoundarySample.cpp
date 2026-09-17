#include "Algorithms/ExperimentalBoundarySample.h"
#include "Ippl.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

class ExperimentalBoundarySampleTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }
    static void TearDownTestSuite() { ippl::finalize(); }
};

TEST_F(ExperimentalBoundarySampleTest, ParsesOnlyUnsignedPowerOfTwoStrides) {
    using experimental_boundary::parseStride;
    EXPECT_EQ(parseStride(nullptr), 1u);
    for (const char* text : {"0", "1", "00", "0001"})
        EXPECT_EQ(parseStride(text), 1u);
    for (unsigned bit = 0; bit < 64; ++bit) {
        const auto stride = UINT64_C(1) << bit;
        EXPECT_EQ(parseStride(std::to_string(stride).c_str()), stride);
    }
    for (const char* text :
         {"", " ", " 2", "2 ", "\t4", "8\n", "-1", "+2", "1.0", "0x10", "2e3", "3", "6",
          "18446744073709551615", "18446744073709551616", "9223372036854775809",
          "999999999999999999999999999999999999"})
        EXPECT_EQ(parseStride(text), 0u) << text;
}

TEST_F(ExperimentalBoundarySampleTest, FixedHashAnchorsPreserveAllIdBits) {
    const std::array<std::array<std::uint64_t, 2>, 8> anchors{
            {{{UINT64_C(0), UINT64_C(0xe220a8397b1dcdaf)}},
             {{UINT64_C(1), UINT64_C(0x910a2dec89025cc1)}},
             {{UINT64_C(2), UINT64_C(0x975835de1c9756ce)}},
             {{UINT64_C(3), UINT64_C(0x1d0b14e4db018fed)}},
             {{UINT64_C(17), UINT64_C(0x808475f02ee37363)}},
             {{UINT64_C(0x8000000000000000), UINT64_C(0x481ec0a212a9f3db)}},
             {{UINT64_C(0x7fffffffffffffff), UINT64_C(0x2a67d7552e039ea7)}},
             {{UINT64_C(0xffffffffffffffff), UINT64_C(0xe4d971771b652c20)}}}};
    Kokkos::View<std::uint64_t*> ids("hash anchor IDs", anchors.size());
    Kokkos::View<std::uint64_t*> hashes("hash anchor values", anchors.size());
    auto host = Kokkos::create_mirror_view(ids);
    for (std::size_t i = 0; i < anchors.size(); ++i)
        host(i) = anchors[i][0];
    Kokkos::deep_copy(ids, host);
    Kokkos::parallel_for(
            "boundary sample hash anchors", anchors.size(),
            KOKKOS_LAMBDA(std::size_t i) { hashes(i) = experimental_boundary::hash(ids(i)); });
    const auto result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), hashes);
    for (std::size_t i = 0; i < anchors.size(); ++i) {
        EXPECT_EQ(experimental_boundary::hash(anchors[i][0]), anchors[i][1]);
        EXPECT_EQ(result(i), anchors[i][1]);
    }
}

TEST_F(ExperimentalBoundarySampleTest, FullPopulationAndNestedMembershipAreDeterministic) {
    using experimental_boundary::selected;
    for (std::uint64_t id = 0; id < 4096; ++id) {
        EXPECT_TRUE(selected(id, 1));
        for (unsigned bit = 1; bit < 64; ++bit) {
            const auto stride = UINT64_C(1) << bit;
            if (selected(id, stride)) EXPECT_TRUE(selected(id, stride / 2));
        }
    }
    EXPECT_TRUE(selected(std::numeric_limits<std::uint64_t>::max(), 1));
    // Membership is not an ID prefix, and ID zero is not a substitute for the
    // separate always-tracked reference. Small samples may select no particles.
    EXPECT_FALSE(selected(0, 2));
    EXPECT_TRUE(selected(2, 2));
    EXPECT_FALSE(selected(3, 2));
}

TEST_F(ExperimentalBoundarySampleTest, DeviceSelectionSurvivesReorderingAndRankOwnershipChanges) {
    // Explicit stable IDs, rather than a FROMFILE row-order assumption: native
    // IDs are assigned rank-strided, whereas FROMFILE distributes row blocks.
    constexpr int population = 1031;
    const int rank = ippl::Comm->rank(), ranks = ippl::Comm->size();
    const auto communicator = ippl::Comm->getCommunicator();
    for (const std::uint64_t stride : {UINT64_C(1), UINT64_C(8), UINT64_C(32)}) {
        std::vector<int> expected(population);
        for (int index = 0; index < population; ++index)
            expected[index] = experimental_boundary::selected(
                    UINT64_C(0x100000000) + 17 * static_cast<std::uint64_t>(index), stride);
        for (int arrangement = 0; arrangement < 3; ++arrangement) {
            std::vector<int> owned;
            for (int index = 0; index < population; ++index) {
                const int owner =
                        arrangement == 0 ? index % ranks
                        : arrangement == 1
                                ? static_cast<int>(
                                          static_cast<long long>(index) * ranks / population)
                                : (index + 1) % ranks;
                if (owner == rank) owned.push_back(index);
            }
            if (arrangement != 0) std::reverse(owned.begin(), owned.end());
            Kokkos::View<std::int64_t*> ids("migrated stable IDs", owned.size());
            Kokkos::View<int*> indices("migrated logical particle indices", owned.size());
            Kokkos::View<int*> flags("local boundary membership", population);
            auto idHost    = Kokkos::create_mirror_view(ids);
            auto indexHost = Kokkos::create_mirror_view(indices);
            for (std::size_t i = 0; i < owned.size(); ++i) {
                idHost(i)    = INT64_C(0x100000000) + 17 * owned[i];
                indexHost(i) = owned[i];
            }
            Kokkos::deep_copy(ids, idHost);
            Kokkos::deep_copy(indices, indexHost);
            Kokkos::deep_copy(flags, 0);
            Kokkos::parallel_for(
                    "partitioned boundary sample", owned.size(), KOKKOS_LAMBDA(std::size_t i) {
                        flags(indices(i)) = experimental_boundary::selected(ids(i), stride) ? 1 : 0;
                    });
            const auto flagHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), flags);
            std::vector<int> global(population);
            MPI_Allreduce(
                    flagHost.data(), global.data(), population, MPI_INT, MPI_SUM, communicator);
            EXPECT_EQ(global, expected) << "stride=" << stride << ", arrangement=" << arrangement;
        }
    }
}
