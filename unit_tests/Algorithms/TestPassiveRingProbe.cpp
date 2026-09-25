#include "Algorithms/PassiveRingProbe.h"
#include "Ippl.h"
#include "PartBunch/ParticleContainer.hpp"
#include "gtest/gtest.h"

#include <unistd.h>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

namespace {
    using Vector = passive_probe::Vector;

    class RestoreEnvironment {
    public:
        RestoreEnvironment() {
            const char* value = std::getenv("OPALX_TEST_PASSIVE_RING_PROBES");
            existed           = value != nullptr;
            if (value) previous = value;
        }
        ~RestoreEnvironment() {
            if (existed)
                setenv("OPALX_TEST_PASSIVE_RING_PROBES", previous.c_str(), 1);
            else
                unsetenv("OPALX_TEST_PASSIVE_RING_PROBES");
        }

    private:
        bool existed;
        std::string previous;
    };

    std::vector<std::vector<std::string>> readCsv(const std::string& path) {
        std::ifstream file(path);
        std::vector<std::vector<std::string>> rows;
        std::string line;
        std::getline(file, line);  // Header.
        while (std::getline(file, line)) {
            std::vector<std::string> row;
            std::istringstream fields(line);
            std::string field;
            while (std::getline(fields, field, ','))
                row.push_back(field);
            rows.push_back(row);
        }
        return rows;
    }
}  // namespace

class PassiveRingProbeTest : public ::testing::Test {
protected:
    inline static std::string directory;
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        if (ippl::Comm->rank() == 0) {
            directory = (std::filesystem::temp_directory_path()
                         / ("opalx-passive-probe-test-" + std::to_string(getpid())))
                                .string();
            std::filesystem::create_directory(directory);
        }
        int length = static_cast<int>(directory.size());
        MPI_Bcast(&length, 1, MPI_INT, 0, ippl::Comm->getCommunicator());
        directory.resize(length);
        MPI_Bcast(directory.data(), length, MPI_CHAR, 0, ippl::Comm->getCommunicator());
    }
    static void TearDownTestSuite() {
        MPI_Barrier(ippl::Comm->getCommunicator());
        if (ippl::Comm->rank() == 0) std::filesystem::remove_all(directory);
        ippl::finalize();
    }
    static PassiveRingProbe::Config config(
            const std::string& name, const std::string& ids = "17 29") {
        return PassiveRingProbe::parseConfig(
                "OUTPUT \"" + directory + "/" + name + ".csv\"\nIDS " + ids
                + "\nPLANE 8 0 0 0 0 0 1\n");
    }
    static std::shared_ptr<ParticleContainer_t> container() {
        ippl::NDIndex<3> domain;
        for (unsigned d = 0; d < 3; ++d)
            domain[d] = ippl::Index(8);
        Mesh_t<3> mesh(domain, Vector(1), Vector(-4));
        std::array<bool, 3> decomposition{true, true, true};
        FieldLayout_t<3> layout(ippl::Comm->getCommunicator(), domain, decomposition, true);
        auto pc = std::make_shared<ParticleContainer_t>(mesh, layout);
        pc->setBunchStateHandler(std::make_shared<BunchStateHandler>());
        pc->createParticles(3);
        return pc;
    }
    static void populate(
            const std::shared_ptr<ParticleContainer_t>& pc, int arrangement, double labZ,
            bool missing = false, bool duplicate = false) {
        // Different frames on successive calls, same laboratory trajectory. This
        // catches accidental interpolation between moving reference-frame values.
        const double angle = 0.13 * (arrangement + 1);
        const CoordinateSystemTrafo toLab(
                Vector(0.2 * arrangement, -0.1, 0.3),
                Quaternion(std::cos(angle / 2), 0, std::sin(angle / 2), 0));
        pc->setToLabTrafo(toLab);
        auto r         = Kokkos::create_mirror_view(pc->R.getView());
        auto p         = Kokkos::create_mirror_view(pc->P.getView());
        auto ids       = Kokkos::create_mirror_view(pc->ID.getView());
        auto dt        = Kokkos::create_mirror_view(pc->dt.getView());
        const int rank = ippl::Comm->rank(), ranks = ippl::Comm->size();
        for (unsigned i = 0; i < 3; ++i) {
            ids(i) = 1000 + 3 * rank + i;
            r(i)   = Vector(i, rank, -3);
            p(i)   = Vector(1, 2, 3);
            dt(i)  = 1e-12 * (i + 1);
        }
        unsigned local = 0;
        if (!missing)
            for (unsigned index = 0; index < 2; ++index) {
                if (int((index + arrangement) % ranks) != rank) continue;
                // Reverse the local order in alternating arrangements.
                const unsigned slot = arrangement % 2 ? 2 - local++ : local++;
                ids(slot)           = index == 0 ? 17 : 29;
                r(slot) = toLab.transformFrom(Vector(0.01 * index, -0.02 * index, labZ));
                p(slot) = toLab.rotateFrom(Vector(0, 0, 0.3));
            }
        if (duplicate && rank == 0) {
            // Existing ID 17 is on rank arrangement % ranks, so this also tests
            // cross-rank duplicates when arrangement=1 and MPI has two ranks.
            ids(0) = 17;
            ids(1) = 17;
        }
        Kokkos::deep_copy(pc->R.getView(), r);
        Kokkos::deep_copy(pc->P.getView(), p);
        Kokkos::deep_copy(pc->ID.getView(), ids);
        Kokkos::deep_copy(pc->dt.getView(), dt);
    }
};

TEST_F(PassiveRingProbeTest, StrictConfigurationValidatesIdsPlanesAndSyntax) {
    const auto parsed = PassiveRingProbe::parseConfig(
            "# internal selected diagnostics\n\nOUTPUT \"path with space.csv\"\n"
            "IDS 29 0 17\nPLANE 8 1 2 3 0 0 2\nPLANE 9 0 0 0 1 0 0 1e-8\n");
    EXPECT_EQ(parsed.output, "path with space.csv");
    EXPECT_EQ(parsed.ids, (std::vector<std::int64_t>{0, 17, 29}));
    ASSERT_EQ(parsed.planes.size(), 2u);
    EXPECT_EQ(parsed.planes[0].id, 8);
    EXPECT_DOUBLE_EQ(parsed.planes[0].plane.origin(2), 3);
    EXPECT_DOUBLE_EQ(parsed.planes[0].plane.normal(2), 2);
    EXPECT_DOUBLE_EQ(parsed.planes[0].plane.tolerance, 1e-9);
    EXPECT_DOUBLE_EQ(parsed.planes[1].plane.tolerance, 1e-8);
    const std::string prefix = "OUTPUT result.csv\n";
    const std::string plane  = "PLANE 0 0 0 0 0 0 1\n";
    for (const std::string ids :
         {"", "-1", "+1", "1.5", "1e2", "0x10", "1 1", "9223372036854775808", "xyz"})
        EXPECT_THROW(
                PassiveRingProbe::parseConfig(prefix + "IDS " + ids + "\n" + plane),
                std::invalid_argument)
                << ids;
    for (const std::string& text : std::vector<std::string>{
                 "", prefix + plane, "IDS 0\n" + plane, prefix + "IDS 0\n",
                 prefix + "OUTPUT again.csv\nIDS 0\n" + plane, prefix + "IDS 0\nIDS 1\n" + plane,
                 prefix + "IDS 0\n" + plane + plane, "OUTPUT x y\nIDS 0\n" + plane,
                 "OUTPUT \"\"\nIDS 0\n" + plane, prefix + "IDS 0\nPLANE -1 0 0 0 0 0 1\n",
                 prefix + "IDS 0\nPLANE 0 0 0 0 0 0 0\n",
                 prefix + "IDS 0\nPLANE 0 0 0 0 0 0 1 -1\n",
                 prefix + "IDS 0\nPLANE 0 0 0 0 0 0 1 0\n",
                 prefix + "IDS 0\nPLANE 0 0 0 0 0 0 1 1e-9 trailing\n",
                 prefix + "IDS 0\nPLANE 0 0 0 0 0 0\n", prefix + "IDS 0\nPLANE 0 0 0 0 0 0 nan\n",
                 prefix + "IDS 0\nOTHER 0\n" + plane})
        EXPECT_THROW(PassiveRingProbe::parseConfig(text), std::invalid_argument) << text;
}

TEST_F(PassiveRingProbeTest, StableIdsSurviveReorderingOwnershipAndMovingFramesWithoutMutation) {
    const auto settings = config("migration");
    PassiveRingProbe recorder(settings);
    const auto pc = container();
    for (int stage = 0; stage < 2; ++stage) {
        populate(pc, stage, stage == 0 ? -0.25 : 0.75);
        // Independent allocations: HostSpace builds can alias mirror-and-copy
        // views, which would otherwise hide an accidental particle mutation.
        auto rBefore   = Kokkos::create_mirror(pc->R.getView());
        auto pBefore   = Kokkos::create_mirror(pc->P.getView());
        auto dtBefore  = Kokkos::create_mirror(pc->dt.getView());
        auto idsBefore = Kokkos::create_mirror(pc->ID.getView());
        Kokkos::deep_copy(rBefore, pc->R.getView());
        Kokkos::deep_copy(pBefore, pc->P.getView());
        Kokkos::deep_copy(dtBefore, pc->dt.getView());
        Kokkos::deep_copy(idsBefore, pc->ID.getView());
        const auto beforeFrame = pc->getToLabTrafo();
        recorder.observe(pc, double(stage));
        const auto rAfter =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
        const auto pAfter =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
        const auto dtAfter =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->dt.getView());
        const auto idsAfter =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->ID.getView());
        for (unsigned i = 0; i < 3; ++i) {
            EXPECT_EQ(idsAfter(i), idsBefore(i));
            EXPECT_DOUBLE_EQ(dtAfter(i), dtBefore(i));
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_DOUBLE_EQ(rAfter(i)(d), rBefore(i)(d));
                EXPECT_DOUBLE_EQ(pAfter(i)(d), pBefore(i)(d));
            }
        }
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_DOUBLE_EQ(pc->getToLabTrafo().getOrigin()(d), beforeFrame.getOrigin()(d));
            EXPECT_DOUBLE_EQ(
                    pc->getToLabTrafo().transformTo(Vector(1))(d),
                    beforeFrame.transformTo(Vector(1))(d));
        }
    }
    if (ippl::Comm->rank() == 0) {
        const auto rows = readCsv(settings.output);
        ASSERT_EQ(rows.size(), 2u);
        for (unsigned i = 0; i < 2; ++i) {
            ASSERT_EQ(rows[i].size(), 14u);
            EXPECT_EQ(std::stoll(rows[i][0]), i == 0 ? 17 : 29);
            EXPECT_EQ(rows[i][1], "8");
            EXPECT_EQ(rows[i][2], "1");
            EXPECT_NEAR(std::stod(rows[i][3]), 0.25, 3e-16);
            EXPECT_NEAR(std::stod(rows[i][4]), 0.01 * i, 2e-16);
            EXPECT_NEAR(std::stod(rows[i][5]), -0.02 * i, 2e-16);
            EXPECT_NEAR(std::stod(rows[i][6]), 0, 2e-16);
            EXPECT_NEAR(std::stod(rows[i][9]), 0.3, 2e-16);
            EXPECT_EQ(rows[i][10], "crossing");
            EXPECT_DOUBLE_EQ(std::stod(rows[i][11]), 0);
            EXPECT_DOUBLE_EQ(std::stod(rows[i][12]), 1);
            EXPECT_NEAR(std::stod(rows[i][13]), 0.25, 3e-16);
        }
    }
}

TEST_F(PassiveRingProbeTest, MissingIdsAreExplicitAndCannotBridgeUnobservedIntervals) {
    const auto settings = config("missing");
    PassiveRingProbe recorder(settings);
    const auto pc = container();
    populate(pc, 0, -1);
    recorder.observe(pc, 0);
    populate(pc, 1, 0, true);
    recorder.observe(pc, 1);
    recorder.observe(pc, 2);  // Still missing: no repeated records.
    populate(pc, 3, 1);
    recorder.observe(pc, 3);  // Reappears beyond the plane: must not count a turn.
    populate(pc, 4, -1);
    recorder.observe(pc, 4);
    populate(pc, 5, 1);
    recorder.observe(pc, 5);
    if (ippl::Comm->rank() == 0) {
        const auto rows = readCsv(settings.output);
        ASSERT_EQ(rows.size(), 4u);
        EXPECT_EQ(rows[0][10], "missing");
        EXPECT_EQ(rows[1][10], "missing");
        EXPECT_TRUE(std::isnan(std::stod(rows[0][4])));
        EXPECT_EQ(rows[2][10], "crossing");
        EXPECT_EQ(rows[3][10], "crossing");
        EXPECT_EQ(rows[2][2], "1");
        EXPECT_NEAR(std::stod(rows[2][3]), 4.5, 1e-15);
    }
}

TEST_F(PassiveRingProbeTest, CollectiveErrorsRejectDuplicatesInvalidTimesAndExistingOutput) {
    const auto settings = config("errors");
    PassiveRingProbe recorder(settings);
    EXPECT_THROW({ PassiveRingProbe duplicate(settings); }, std::runtime_error);
    const auto pc = container();
    populate(pc, 0, -1);
    recorder.observe(pc, 1);
    EXPECT_NO_THROW(recorder.observe(pc, 1));
    try {
        recorder.observe(pc, 0.5);
        FAIL() << "Decreasing observation time must be rejected";
    } catch (const std::runtime_error& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("particle ID 17"), std::string::npos);
        EXPECT_NE(message.find("plane ID 8"), std::string::npos);
        EXPECT_NE(message.find("status=DecreasingTime"), std::string::npos);
        EXPECT_NE(message.find("previous_time_s=1"), std::string::npos);
        EXPECT_NE(message.find("current_time_s=0.5"), std::string::npos);
    }
    EXPECT_THROW(
            recorder.observe(pc, std::numeric_limits<double>::quiet_NaN()), std::runtime_error);
    EXPECT_THROW(recorder.observe(nullptr, 2), std::runtime_error);
    if (ippl::Comm->size() > 1)
        EXPECT_THROW(recorder.observe(pc, 2 + ippl::Comm->rank()), std::runtime_error);
    populate(pc, 1, 1, false, true);
    EXPECT_THROW(recorder.observe(pc, 2), std::runtime_error);
}

TEST_F(PassiveRingProbeTest, EnvironmentIsOptInAndEligibilityIsCheckedCollectively) {
    RestoreEnvironment restore;
    unsetenv("OPALX_TEST_PASSIVE_RING_PROBES");
    EXPECT_EQ(PassiveRingProbe::fromEnvironment(false), nullptr);
    setenv("OPALX_TEST_PASSIVE_RING_PROBES", "", 1);
    EXPECT_THROW(PassiveRingProbe::fromEnvironment(true), std::runtime_error);
    const std::string path = directory + "/config.txt";
    setenv("OPALX_TEST_PASSIVE_RING_PROBES", path.c_str(), 1);
    EXPECT_THROW(PassiveRingProbe::fromEnvironment(true), std::runtime_error);
    if (ippl::Comm->rank() == 0) {
        std::ofstream file(path);
        file << "OUTPUT \"" << directory
             << "/environment.csv\"\n"
                "IDS 17\nPLANE 0 0 0 0 0 0 1\n";
    }
    MPI_Barrier(ippl::Comm->getCommunicator());
    EXPECT_THROW(PassiveRingProbe::fromEnvironment(ippl::Comm->rank() != 0), std::runtime_error);
    const auto recorder = PassiveRingProbe::fromEnvironment(true);
    EXPECT_NE(recorder, nullptr);
}
