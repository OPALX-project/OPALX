//
// Tests of the step output of DataSink: which outputs are due at a step, the moments they read,
// and that writeStep() writes handed-in moments instead of computing them again.
//
#include <gtest/gtest.h>

#include "Ippl.h"

#include "../KeywordObjects.h"
#include "../TemporaryDirectory.h"
#include "Moments/BeamMoments.h"
#include "PartBunch/PartBunch.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"

#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

extern Inform* gmsg;

namespace {
    using opalx::moments::Moment;
    using Vector = Vector_t<double, 3>;

    class DataSinkStepOutputTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
            gmsg = new Inform(nullptr, -1);
        }
        static void TearDownTestSuite() {
            delete gmsg;
            gmsg = nullptr;
            ippl::finalize();
        }

        void SetUp() override {
            directory_m.create("opalx-data-sink-step");
            const auto beam          = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
            auto withSpin            = opalx::run::describeBeamForBunch(*beam);
            auto withoutSpin         = withSpin;
            withSpin.hasPolarization = true;
            bunch_m                  = std::make_unique<PartBunch_t>(
                    std::vector<double>(2, 1.0), std::vector<double>(2, 1.0),
                    std::vector{withoutSpin, withSpin}, std::vector<size_t>(2, 8), 1.0, "LF2",
                    opalx::spacecharge::CartesianDomainConfig3D{});
            for (const auto& pc : bunch_m->getParticleContainers())
                place(*pc, {1e-3, 2e-3, 3e-3, 6e-3});
        }

        /// Particles at x = @p x [m] moving along z.
        static void place(ParticleContainer_t& pc, const std::vector<double>& x) {
            if (pc.getLocalNum() < x.size()) pc.createParticles(x.size() - pc.getLocalNum());
            auto R = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc.R.getView());
            auto P = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc.P.getView());
            for (size_t i = 0; i < x.size(); ++i) {
                R(i) = Vector({x[i], 0.0, 0.0});
                P(i) = Vector({0.0, 0.0, 1.0});
            }
            Kokkos::deep_copy(pc.R.getView(), R);
            Kokkos::deep_copy(pc.P.getView(), P);
        }

        /// The data rows of SDDS file @p file, by column name.
        static std::vector<std::map<std::string, double>> readRows(const std::string& file) {
            std::ifstream in(file);
            std::vector<std::string> names;
            std::vector<std::map<std::string, double>> rows;
            bool column = false, data = false;
            for (std::string line; std::getline(in, line);) {
                if (line == "&column") {
                    column = true;
                } else if (const auto at = line.find("name="); column && at != std::string::npos) {
                    names.push_back(line.substr(at + 5, line.find(',') - at - 5));
                    column = false;
                } else if (line == "&data") {
                    data = true;
                } else if (data) {
                    // The parameter values before the rows have fewer fields.
                    std::istringstream tokens(line);
                    std::vector<std::string> values;
                    for (std::string token; tokens >> token;)
                        values.push_back(token);
                    if (values.size() != names.size()) continue;
                    auto& row = rows.emplace_back();
                    for (size_t i = 0; i < names.size(); ++i)
                        row[names[i]] = std::stod(values[i]);
                }
            }
            return rows;
        }

        opalx::test::TemporaryDirectory directory_m;
        std::unique_ptr<PartBunch_t> bunch_m;
    };

    TEST_F(DataSinkStepOutputTest, OutputsAreDueAtTheirCadence) {
        const DataSink sink(
                {.psDumpFreq = 4, .statDumpFreq = 3, .stepInfoFreq = 5, .printRankDistrFreq = 1},
                {}, 1, (directory_m.path() / "due").string(), "test.in");
        std::set<long long> phaseSpace, statistics, stepLine, rankLoad;
        for (long long step = 0; step <= 10; ++step) {
            const auto due = sink.outputsDue(step);
            if (due.phaseSpace) phaseSpace.insert(step);
            if (due.statistics) statistics.insert(step);
            if (due.stepLine) stepLine.insert(step);
            if (due.rankLoad) rankLoad.insert(step);
        }
        // Dumps and the rank load at the last step of each cycle, the step line at its start.
        EXPECT_EQ(phaseSpace, (std::set<long long>{3, 7}));
        EXPECT_EQ(statistics, (std::set<long long>{2, 5, 8}));
        EXPECT_EQ(stepLine, (std::set<long long>{0, 5, 10}));
        EXPECT_EQ(rankLoad, (std::set<long long>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));

        const DataSink never({}, {}, 1, (directory_m.path() / "never").string(), "test.in");
        for (long long step = 0; step <= 10; ++step) {
            const auto due = never.outputsDue(step);
            EXPECT_FALSE(due.phaseSpace || due.statistics || due.stepLine || due.rankLoad);
        }
    }

    TEST_F(DataSinkStepOutputTest, RequiredMomentsAreWhatComputeProvides) {
        const DataSink sink({}, {}, 2, (directory_m.path() / "moments").string(), "test.in");
        EXPECT_TRUE(sink.requiredMoments({}, false).empty());
        EXPECT_TRUE(sink.requiredMoments({}, true).empty());
        // No phase-space moments without HDF5.
        EXPECT_TRUE(sink.requiredMoments({.phaseSpace = true}, true).empty());
        EXPECT_EQ(sink.requiredMoments({.stepLine = true}, false), Moment::MeanMomentum);

        for (const auto& pc : bunch_m->getParticleContainers()) {
            const auto request = sink.requiredMoments({.statistics = true}, pc->hasSpin());
            EXPECT_EQ(request.contains(Moment::Polarization), pc->hasSpin());
            EXPECT_TRUE(opalx::moments::compute(*pc, request).available().contains(request));
        }
    }

    TEST_F(DataSinkStepOutputTest, WriteStepUsesHandedInMoments) {
        const std::string stem = (directory_m.path() / "write").string();
        DataSink sink({.statDumpFreq = 1}, {}, 2, stem, "test.in");
        auto pc        = bunch_m->getParticleContainer(0);
        const auto due = sink.outputsDue(0);
        ASSERT_TRUE(due.statistics);
        const auto before  = opalx::moments::compute(*pc, sink.requiredMoments(due, false));
        const auto partial = opalx::moments::compute(*pc, Moment::MeanPosition);
        place(*pc, {0.101, 0.102, 0.103, 0.106});

        int fieldCalls        = 0;
        const auto zeroFields = [&fieldCalls] {
            ++fieldCalls;
            return DataSink::ReferenceFields(2, {Vector(0.0), Vector(0.0)});
        };
        sink.writeStep(*bunch_m, {}, {}, zeroFields, 1);
        EXPECT_EQ(fieldCalls, 0);
        sink.writeStep(*bunch_m, due, std::vector{before}, zeroFields, 1);
        sink.writeStep(*bunch_m, due, {}, zeroFields, 1);
        // Moments without what the statistics read are computed again.
        sink.writeStep(*bunch_m, due, std::vector{partial}, zeroFields, 1);
        EXPECT_EQ(fieldCalls, 3);

        const auto rows = readRows(stem + "_c0.stat");
        ASSERT_EQ(rows.size(), 3u);
        EXPECT_NEAR(rows[0].at("mean_x"), 3e-3, 1e-15);
        EXPECT_NEAR(rows[0].at("max_x"), 6e-3, 1e-15);
        for (size_t row : {1, 2}) {
            EXPECT_NEAR(rows[row].at("mean_x"), 0.103, 1e-14);
            EXPECT_NEAR(rows[row].at("max_x"), 0.106, 1e-14);
        }
    }
}  // namespace
