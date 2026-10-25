#include <unistd.h>
#include <array>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include "Algorithms/CyclotronProbeRecorder.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/ProbeRep.h"
#include "Elements/OpalProbe.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"
#include "ValueDefinitions/RealVariable.h"
#include "gtest/gtest.h"

extern Inform* gmsg;

namespace {
    using namespace cyclotron_probe;
    Endpoint point(double x, double z, double px, double pz, double time) {
        return {Vector(x, .2, z), Vector(px, 0, pz), time};
    }
    struct CrossingKernel {
        Kokkos::View<State*> states;
        Kokkos::View<Sample*> samples;
        Kokkos::View<int*> statuses;
        KOKKOS_INLINE_FUNCTION void operator()(int i) const {
            Segment segment{Vector(0, 0, -1), Vector(0, 0, 1)};
            if (i == 1) segment = {Vector(0, 0, 1), Vector(0, 0, -1)};
            if (i == 2) segment = {Vector(-1, 0, -1), Vector(1, 0, 1)};
            if (i == 3) segment = {Vector(-1, 0, 0), Vector(1, 0, 0)};
            Endpoint before{Vector(-1, .2, 0), Vector(1, 0, 0), 0};
            Endpoint after{Vector(1, .2, 0), Vector(1, 0, 0), 2e-9};
            if (i == 3) {
                before.position = Vector(1, .2, -1);
                after.position  = Vector(1, .2, 1);
                before.momentum = after.momentum = Vector(0, 0, 1);
            }
            if (i == 4) {
                before.position(2) = 2;
                after.position(2)  = 2;
            }
            if (i == 5) {
                after.position(0) = -1;
            }
            if (i == 6) {
                before.momentum = after.momentum = Vector(0);
            }
            if (i == 7) {
                before.position = Vector(1, .2, 0);
                after.position  = Vector(-1, .2, 0);
                before.momentum = after.momentum = Vector(-1, 0, 0);
            }
            update(segment, states(i), i, before, samples(i));
            statuses(i) = static_cast<int>(update(segment, states(i), i, after, samples(i)));
        }
    };
    void populatePosition(
            const std::shared_ptr<ParticleContainer_t>& pc, double x, double px, double z) {
        const auto r = pc->R.getView(), p = pc->P.getView();
        Kokkos::parallel_for(
                "test lab endpoints", pc->getLocalNum(), KOKKOS_LAMBDA(std::size_t i) {
                    r(i) = Vector(x, .2, z);
                    p(i) = Vector(px, 0, 0);
                });
    }
    std::vector<std::vector<double>> rows(const std::string& path) {
        std::ifstream file(path);
        std::string line;
        std::vector<std::vector<double>> result;
        while (std::getline(file, line)) {
            if (line.empty() || line.front() == '#') continue;
            std::istringstream input(line);
            std::vector<double> row;
            double value;
            while (input >> value)
                row.push_back(value);
            result.push_back(row);
        }
        return result;
    }
}  // namespace

class CyclotronProbeTest : public ::testing::Test {
protected:
    inline static std::string directory;
    inline static std::unique_ptr<Inform> message;
    inline static Inform* savedMessage = nullptr;
    bool ascii, percentiles;
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        savedMessage = gmsg;
        message      = std::make_unique<Inform>("CyclotronProbeTest");
        gmsg         = message.get();
        if (ippl::Comm->rank() == 0) {
            directory = (std::filesystem::temp_directory_path()
                         / ("opalx-cyclotron-probe-" + std::to_string(getpid())))
                                .string();
            std::filesystem::create_directory(directory);
        }
        int size = directory.size();
        MPI_Bcast(&size, 1, MPI_INT, 0, ippl::Comm->getCommunicator());
        directory.resize(size);
        MPI_Bcast(directory.data(), size, MPI_CHAR, 0, ippl::Comm->getCommunicator());
    }
    static void TearDownTestSuite() {
        ippl::Comm->barrier();
        if (ippl::Comm->rank() == 0) std::filesystem::remove_all(directory);
        gmsg = savedMessage;
        message.reset();
        ippl::finalize();
    }
    void SetUp() override {
        ascii                       = Options::asciidump;
        percentiles                 = Options::computePercentiles;
        Options::computePercentiles = false;
    }
    void TearDown() override {
        Options::asciidump          = ascii;
        Options::computePercentiles = percentiles;
    }
    static CyclotronProbeRecorder::Config config(const std::string& name) {
        return {{Vector(0, 0, -2), Vector(0, 0, 2)}, .1, directory + '/' + name};
    }
    static std::shared_ptr<ParticleContainer_t> container() {
        ippl::NDIndex<3> domain;
        for (unsigned d = 0; d < 3; ++d)
            domain[d] = ippl::Index(8);
        auto mesh = std::make_shared<Mesh_t<3>>(domain, Vector(1), Vector(-4));
        std::array<bool, 3> decomposition{true, false, false};
        auto layout = std::make_shared<FieldLayout_t<3>>(
                ippl::Comm->getCommunicator(), domain, decomposition, true);
        // IPPL's spatial layout borrows the mesh/layout. Retain both through
        // every real migration update, not just container construction.
        auto pc = std::shared_ptr<ParticleContainer_t>(
                new ParticleContainer_t(*mesh, *layout),
                [mesh, layout](ParticleContainer_t* value) {
                    delete value;
                });
        pc->setBunchStateHandler(std::make_shared<BunchStateHandler>());
        pc->createParticles(ippl::Comm->rank() == 0 ? 1 : 0);
        pc->setQ(1.6e-19);
        pc->setM(.938);
        return pc;
    }
    static void position(
            const std::shared_ptr<ParticleContainer_t>& pc, double x, double px, double z = .25) {
        populatePosition(pc, x, px, z);
    }
};

TEST_F(CyclotronProbeTest, DeviceCrossingsClippingDirectionsAndIds) {
    // Actual device path: orthogonal and oblique segments, reversed endpoints,
    // boundary inclusion, finite-span miss, parallel motion and zero momentum.
    Kokkos::View<State*> states("test histories", 8);
    Kokkos::View<Sample*> samples("test samples", 8);
    Kokkos::View<int*> statuses("test statuses", 8);
    Kokkos::parallel_for("test segment update", 8, CrossingKernel{states, samples, statuses});
    auto status = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), statuses);
    auto hits   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), samples);
    for (int i : {0, 1, 2, 3, 7}) {
        EXPECT_EQ(status(i), static_cast<int>(Status::Crossing));
        EXPECT_NEAR(hits(i).crossing.time, 1e-9, 1e-24);
        EXPECT_DOUBLE_EQ(hits(i).crossing.position(1), .2);
    }
    for (int i : {4, 5, 6})
        EXPECT_EQ(status(i), static_cast<int>(Status::Advanced));
}

TEST_F(CyclotronProbeTest, EndpointDuplicatesInitialPlaneJitterAndInvalidInput) {
    Segment segment{Vector(0, 0, -1), Vector(0, 0, 1)};
    State state;
    Sample sample;
    EXPECT_EQ(update(segment, state, 17, point(0, 0, 1, 0, 0), sample), Status::Initialized);
    EXPECT_EQ(update(segment, state, 17, point(1, 0, -1, 0, 1), sample), Status::Advanced);
    EXPECT_EQ(update(segment, state, 17, point(-1, 0, -1, 0, 2), sample), Status::Crossing);
    EXPECT_EQ(update(segment, state, 17, point(-1, 0, 1, 0, 2.5), sample), Status::Advanced);
    EXPECT_EQ(update(segment, state, 17, point(0, 0, 1, 0, 3), sample), Status::Crossing);
    EXPECT_EQ(sample.turn, 2u);
    EXPECT_EQ(update(segment, state, 17, point(1, 0, 1, 0, 4), sample), Status::Advanced);
    EXPECT_EQ(update(segment, state, 17, point(-1e-12, 0, -1, 0, 5), sample), Status::Crossing);
    EXPECT_EQ(update(segment, state, 17, point(1e-12, 0, 1, 0, 6), sample), Status::Advanced);
    EXPECT_EQ(update(segment, state, 18, point(-1, 0, 1, 0, 7), sample), Status::Initialized);
    EXPECT_EQ(state.passages, 0);
    auto bad = point(1, 0, 1, 0, 6);
    EXPECT_EQ(update(segment, state, 18, bad, sample), Status::DecreasingTime);
    EXPECT_DOUBLE_EQ(state.forward.previous.time, 7);
    bad.time        = 8;
    bad.position(0) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(update(segment, state, 18, bad, sample), Status::InvalidInput);
    state.passages = INT_MAX;
    EXPECT_EQ(update(segment, state, 18, point(1, 0, 1, 0, 8), sample), Status::TurnOverflow);
    EXPECT_DOUBLE_EQ(state.forward.previous.time, 7);
    segment.end = segment.start;
    EXPECT_EQ(update(segment, state, 18, point(1, 0, 1, 0, 8), sample), Status::InvalidInput);
}

TEST_F(CyclotronProbeTest, StationaryRefreshPreservesArmingAndResumedCrossingTime) {
    Segment segment{Vector(0, 0, -1), Vector(0, 0, 1)};
    State state;
    Sample sample;
    EXPECT_EQ(update(segment, state, 17, point(-1, 0, 1, 0, 0), sample), Status::Initialized);
    EXPECT_EQ(update(segment, state, 17, point(-1e-12, 0, 1, 0, 1e-9), sample), Status::Advanced);
    // Another container advanced the global clock while this one was stopped.
    EXPECT_EQ(update(segment, state, 17, point(-1e-12, 0, 1, 0, 10e-9), sample), Status::Advanced);
    EXPECT_EQ(update(segment, state, 17, point(1e-12, 0, 1, 0, 12e-9), sample), Status::Crossing);
    EXPECT_NEAR(sample.crossing.time, 11e-9, 1e-23);
}

TEST_F(CyclotronProbeTest, ConfigurationValidationAndLegacyMapping) {
    ProbeRep probe("mapping");
    probe.setDimensions(.1, .2, .3, .4);
    probe.setStep(.001);
    probe.setOutputFN(directory + "/mapping.h5");
    auto c = CyclotronProbeRecorder::configuration(probe);
    EXPECT_DOUBLE_EQ(c.segment.start(0), .1);
    EXPECT_DOUBLE_EQ(c.segment.start(2), .3);
    EXPECT_DOUBLE_EQ(c.segment.start(1), 0);
    EXPECT_DOUBLE_EQ(c.step, .001);
    EXPECT_EQ(c.output, directory + "/mapping");
    c      = config("invalid");
    c.step = 0;
    EXPECT_THROW(CyclotronProbeRecorder::validate(c), std::invalid_argument);
    c.step = -1;
    EXPECT_THROW(CyclotronProbeRecorder::validate(c), std::invalid_argument);
    c.step = 1e-30;
    EXPECT_THROW(CyclotronProbeRecorder::validate(c), std::invalid_argument);
    c             = config("invalid");
    c.segment.end = c.segment.start;
    EXPECT_THROW(CyclotronProbeRecorder::validate(c), std::invalid_argument);
    c                = config("invalid");
    c.segment.end(0) = std::numeric_limits<double>::infinity();
    EXPECT_THROW(CyclotronProbeRecorder::validate(c), std::invalid_argument);
    c = config("duplicate");
    EXPECT_THROW(CyclotronProbeRecorder({c, c}), std::invalid_argument);
}

TEST_F(CyclotronProbeTest, InputConvertsMillimetresAndIgnoresPathLength) {
    [[maybe_unused]] static RealVariable* p0Owner = new RealVariable();
    OpalProbe prototype;
    EXPECT_NO_THROW(prototype.update());
    std::unique_ptr<OpalProbe> probe(prototype.clone("PARSER_PROBE"));
    Attributes::setReal(*probe->findAttribute("XSTART"), 100);
    Attributes::setReal(*probe->findAttribute("XEND"), 200);
    Attributes::setReal(*probe->findAttribute("YSTART"), 300);
    Attributes::setReal(*probe->findAttribute("YEND"), 400);
    Attributes::setReal(*probe->findAttribute("STEP"), 2);
    Attributes::setReal(*probe->findAttribute("L"), 1);
    EXPECT_NO_THROW(probe->update());
    const auto* representation = dynamic_cast<const Probe*>(probe->getElement());
    ASSERT_NE(representation, nullptr);
    auto c = CyclotronProbeRecorder::configuration(*representation);
    EXPECT_DOUBLE_EQ(c.segment.start(0), .1);
    EXPECT_DOUBLE_EQ(c.segment.start(2), .3);
    EXPECT_DOUBLE_EQ(c.step, .002);
    EXPECT_DOUBLE_EQ(probe->getLength(), 0);
    EXPECT_DOUBLE_EQ(representation->getGeometry().getElementLength(), 0);
    Attributes::setReal(*probe->findAttribute("STEP"), 0);
    EXPECT_THROW(probe->update(), std::invalid_argument);
}

TEST_F(CyclotronProbeTest, AsciiEmptyRanksMigrationFramesAndTurnGroups) {
    Options::asciidump    = true;
    const auto c          = config("ascii");
    auto miss             = config("miss");
    miss.segment.start(0) = miss.segment.end(0) = 3;
    auto pc                                     = container();
    position(pc, -1, 1);
    CyclotronProbeRecorder recorder({c, miss});
    recorder.initialize({pc});
    recorder.observe(pc, 0, 2, 3);
    // Force IPPL ownership updates across the X decomposition with registered
    // histories. The owning rank changes; empty ranks still participate in save.
    pc->update();
    position(pc, 1, 1);
    pc->update();
    recorder.observe(pc, 2e-9, 2, 3);
    // Reverse the momentum at the accepted endpoint before the return drift.
    // A single interval with +p and -p would interpolate zero crossing momentum.
    position(pc, 1, -1);
    recorder.observe(pc, 2e-9, 2, 3);
    // Change tracking frame while preserving the same laboratory trajectory.
    CoordinateSystemTrafo frame(
            Vector(.3, -.4, .5), Quaternion(std::cos(.15), Vector(0, std::sin(.15), 0)));
    position(pc, -1, -1, .85);
    auto r = pc->R.getView(), p = pc->P.getView();
    frame.inverted().transformBunchTo(r, pc->getLocalNum());
    frame.inverted().rotateBunchTo(p, pc->getLocalNum());
    pc->setToLabTrafo(frame);
    recorder.observe(pc, 4e-9, 2, 7);
    const auto before = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), r);
    recorder.save();
    const auto after = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), r);
    for (std::size_t i = 0; i < pc->getLocalNum(); ++i)
        for (int d = 0; d < 3; ++d)
            EXPECT_DOUBLE_EQ(before(i)(d), after(i)(d));
    if (ippl::Comm->rank() == 0) {
        auto data = rows(c.output + ".loss");
        ASSERT_EQ(data.size(), 2u);
        for (std::size_t i = 0; i < 2; ++i) {
            ASSERT_EQ(data[i].size(), 10u);
            EXPECT_NEAR(data[i][0], 0, 1e-14);
            EXPECT_NEAR(data[i][1], .2, 1e-14);
            EXPECT_NEAR(data[i][2], i == 0 ? .25 : .55, 1e-14);
            EXPECT_EQ(data[i][7], i == 0 ? 3 : 7);
            EXPECT_EQ(data[i][8], 2);
            EXPECT_NEAR(data[i][9], (1 + 2 * i) * 1e-9, 1e-15);
        }
        auto hist = rows(c.output + ".hist");
        ASSERT_EQ(hist.size(), 21u);
        EXPECT_EQ(hist[2][0], 1);
        EXPECT_EQ(hist[5][0], 1);
        double total = 0;
        for (const auto& row : hist)
            total += row[0];
        EXPECT_EQ(total, 2);
        auto peaks = rows(c.output + ".peaks");
        ASSERT_EQ(peaks.size(), 2u);
        EXPECT_NEAR(peaks[0][0], 250, 1e-10);
        EXPECT_NEAR(peaks[1][0], 550, 1e-10);
        EXPECT_FALSE(std::filesystem::exists(miss.output + ".loss"));
        EXPECT_TRUE(rows(miss.output + ".peaks").empty());
        for (const auto& row : rows(miss.output + ".hist"))
            EXPECT_EQ(row[0], 0);
    }
}

TEST_F(CyclotronProbeTest, CapacityGrowthPreservesArmedHistoryAndNewParticleIds) {
    Options::asciidump = false;
    const auto c       = config("growth");
    auto pc            = container();
    position(pc, -1, 1);
    const auto ids = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->ID.getView());
    std::int64_t expectedId = ippl::Comm->rank() == 0 ? ids(0) : 0;
    MPI_Bcast(&expectedId, 1, MPI_INT64_T, 0, ippl::Comm->getCommunicator());
    CyclotronProbeRecorder recorder({c});
    recorder.initialize({pc});
    recorder.observe(pc, 0, 0, 0);

    // Grow registered history after arming a particle. Reallocation must keep
    // its previous endpoint; new IDs must initialize without inventing a hit.
    const auto oldCapacity = pc->R.size();
    pc->createParticles(ippl::Comm->rank() == 0 ? oldCapacity + 1 : 0);
    if (ippl::Comm->rank() == 0) EXPECT_GT(pc->R.size(), oldCapacity);
    position(pc, 1, 1);
    recorder.observe(pc, 2e-9, 0, 0);
    recorder.save();

    h5_prop_t props = H5CreateFileProp();
    auto comm       = ippl::Comm->getCommunicator();
    H5SetPropFileMPIOCollective(props, &comm);
    auto file = H5OpenFile((c.output + ".h5").c_str(), H5_O_RDONLY, props);
    H5CloseProp(props);
    ASSERT_NE(file, static_cast<h5_file_t>(H5_ERR));
    EXPECT_EQ(H5GetNumSteps(file), 1);
    ASSERT_EQ(H5SetStep(file, 0), H5_SUCCESS);
    // Read the global count before setting a one-particle view on every rank.
    EXPECT_EQ(H5PartGetNumParticles(file), 1);
    ASSERT_EQ(H5PartSetView(file, 0, 0), H5_SUCCESS);
    h5_int64_t readId = -1;
    EXPECT_EQ(H5PartReadDataInt64(file, "id", &readId), H5_SUCCESS);
    EXPECT_EQ(readId, expectedId);
    double time = 0;
    EXPECT_EQ(H5PartReadDataFloat64(file, "time", &time), H5_SUCCESS);
    EXPECT_NEAR(time, 1e-9, 1e-24);
    H5CloseFile(file);
}

TEST_F(CyclotronProbeTest, CollectiveErrorsAndFreshOutputGuard) {
    Options::asciidump = true;
    auto c             = config("errors");
    auto pc            = container();
    position(pc, -1, 1);
    CyclotronProbeRecorder recorder({c});
    recorder.initialize({pc});
    recorder.observe(pc, 0, 0, 0);
    position(pc, std::numeric_limits<double>::quiet_NaN(), 1);
    // Only one rank owns the invalid particle; all ranks must receive the error.
    EXPECT_THROW(recorder.observe(pc, 1e-9, 0, 0), std::runtime_error);
    auto unavailable = config("nonexistent/subdir/output");
    CyclotronProbeRecorder noOutput({unavailable});
    noOutput.initialize({pc});
    EXPECT_THROW(noOutput.save(), std::runtime_error);
    auto stale = config("stale");
    if (ippl::Comm->rank() == 0) {
        std::ofstream old(stale.output + ".h5");
        old << "old output";
    }
    ippl::Comm->barrier();
    EXPECT_THROW(CyclotronProbeRecorder({stale}), std::runtime_error);
    auto alias   = c;
    alias.output = directory + "/./errors";
    EXPECT_THROW(CyclotronProbeRecorder({c, alias}), std::invalid_argument);
}

TEST_F(CyclotronProbeTest, Hdf5ContainsCrossingUnitsChargeMassAndStableId) {
    Options::asciidump = false;
    auto c             = config("hdf5");
    auto pc            = container();
    position(pc, -1, 1);
    auto id = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->ID.getView());
    std::int64_t expectedId = ippl::Comm->rank() == 0 ? id(0) : 0;
    MPI_Bcast(&expectedId, 1, MPI_INT64_T, 0, ippl::Comm->getCommunicator());
    CyclotronProbeRecorder recorder({c});
    recorder.initialize({pc});
    recorder.observe(pc, 0, 0, 0);
    position(pc, 1, 1);
    recorder.observe(pc, 2e-9, 0, 0);
    recorder.save();
    recorder.save();  // No new hits: no duplicate HDF5 step.
    position(pc, 1, -1);
    recorder.observe(pc, 3e-9, 0, 1);  // Arm reverse direction away from the plane.
    position(pc, -1, -1);
    recorder.observe(pc, 4e-9, 0, 1);
    recorder.save();
    // Read back collectively, each rank reading its local output slice.
    h5_prop_t props = H5CreateFileProp();
    auto comm       = ippl::Comm->getCommunicator();
    H5SetPropFileMPIOCollective(props, &comm);
    auto file = H5OpenFile((c.output + ".h5").c_str(), H5_O_RDONLY, props);
    H5CloseProp(props);
    ASSERT_NE(file, static_cast<h5_file_t>(H5_ERR));
    ASSERT_EQ(H5SetStep(file, 0), H5_SUCCESS);
    H5PartSetView(file, 0, 0);
    double value;
    for (const auto& field : std::vector<std::pair<std::string, double>>{
                 {"x", 0},
                 {"y", .2},
                 {"z", .25},
                 {"px", 1},
                 {"time", 1e-9},
                 {"q", 1.6e-19},
                 {"m", .938}}) {
        ASSERT_EQ(H5PartReadDataFloat64(file, field.first.c_str(), &value), H5_SUCCESS);
        EXPECT_NEAR(value, field.second, std::max(1e-30, std::abs(field.second) * 1e-13));
    }
    h5_int64_t readId = -1;
    EXPECT_EQ(H5PartReadDataInt64(file, "id", &readId), H5_SUCCESS);
    EXPECT_EQ(readId, expectedId);
    EXPECT_EQ(H5GetNumSteps(file), 2);
    ASSERT_EQ(H5SetStep(file, 1), H5_SUCCESS);
    EXPECT_EQ(H5PartReadDataInt64(file, "id", &readId), H5_SUCCESS);
    EXPECT_EQ(readId, expectedId);
    H5CloseFile(file);
}

TEST_F(CyclotronProbeTest, RepeatedFlushAppendsHitsAndRetainsCumulativeSummaries) {
    Options::asciidump = true;
    auto c             = config("flush");
    auto pc            = container();
    position(pc, -1, 1);
    CyclotronProbeRecorder recorder({c});
    recorder.initialize({pc});
    recorder.save();  // Empty flush must not prevent later file creation.
    recorder.observe(pc, 0, 0, 0);
    position(pc, 1, 1);
    recorder.observe(pc, 2e-9, 0, 0);
    recorder.save();
    EXPECT_EQ(rows(c.output + ".loss").size(), 1u);
    recorder.save();
    EXPECT_EQ(rows(c.output + ".loss").size(), 1u);
    position(pc, 1, -1);
    recorder.observe(pc, 3e-9, 0, 1);  // Arm reverse direction away from the plane.
    position(pc, -1, -1);
    recorder.observe(pc, 4e-9, 0, 1);
    recorder.save();
    EXPECT_EQ(rows(c.output + ".loss").size(), 2u);
    double total = 0;
    for (const auto& row : rows(c.output + ".hist"))
        total += row[0];
    EXPECT_EQ(total, 2);
    EXPECT_EQ(rows(c.output + ".peaks").size(), 2u);
}
