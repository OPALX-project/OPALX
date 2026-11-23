//
// End-to-end test of ParallelTracker: a RUN of four particles through a drift, built as RUN
// builds it (freeze, validate, buildElements, buildTracking), then tracked.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "../TemporaryDirectory.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/ParallelTracker.h"
#include "Attributes/Attributes.h"
#include "Ippl.h"
#include "OpalConfigure/Configure.h"
#include "OpalParser/OpalParser.h"
#include "OpalParser/StringStream.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"
#include "Track/BuildElements.h"
#include "Track/BuildTracking.h"
#include "Track/FreezeRun.h"
#include "Track/TrackRun.h"
#include "Track/ValidateRun.h"
#include "Utility/Inform.h"

extern Inform* gmsg;

namespace {
    using Vector = Vector_t<double, 3>;

    // Centred on the origin, so the reference particle starts at the origin of its frame.
    const std::array<std::array<double, 6>, 4> particles{
            {{1e-3, 0.0, 1e-3, 0.005, 0.0, 0.1066},
             {-1e-3, 0.0, -1e-3, -0.005, 0.0, 0.1066},
             {0.0, 1e-3, 1e-3, 0.0, 0.005, 0.1066},
             {0.0, -1e-3, -1e-3, 0.0, -0.005, 0.1066}}};

    constexpr double dt = 1e-10;

    class ParallelTrackerTest : public ::testing::Test {
    protected:
        void SetUp() override {
            gmsg = new Inform("", output_m, INFORM_ALL_NODES);
            Configure::configure();
            directory_m.create("opalx-parallel-tracker");
            OpalData::getInstance()->storeInputFn((directory_m.path() / "test.in").string());

            const std::string file = (directory_m.path() / "particles.dat").string();
            std::ofstream out(file);
            out << particles.size() << "\nx y z px py pz\n";
            for (const auto& row : particles) {
                for (double value : row)
                    out << value << ' ';
                out << '\n';
            }
            out.close();

            // A 0.1 m drift, the four particles above as protons and no field solver.
            parse("OPTION, VERSION = 10900; OPTION, ENABLEHDF5 = FALSE;"
                  "D1: DRIFT, L = 0.1, ELEMEDGE = 0; L1: LINE = (D1);"
                  "DI: DISTRIBUTION, TYPE = FROMFILE, FNAME = \""
                  + file
                  + "\", NPARTDIST = 4;"
                    "ES1: EMISSIONSOURCE, DISTRIBUTION = DI; SRC: EMISSIONSOURCELIST = (ES1);"
                    "FS1: FIELDSOLVER, TYPE = NONE, NX = 8, NY = 8, NZ = 8;"
                    "B1: BEAM, PARTICLE = PROTON, PC = 0.1, NALLOC = 4, BCHARGE = 1e-15, "
                    "SOURCES = SRC;");
        }

        void TearDown() override {
            // The run objects log to gmsg when they are destroyed.
            objects_m.reset();
            OpalData::deleteInstance();
            delete gmsg;
            gmsg = nullptr;
        }

        static void parse(const std::string& text) {
            OpalParser parser;
            parser.run(new StringStream(text));
        }

        /// TRACK, LINE = L1, BEAM = B1, DT = dt, MAXSTEPS = @p maxSteps[, ZSTOP = @p zStop];
        /// RUN, METHOD = "PARALLEL", FIELDSOLVER = FS1;
        static opalx::run::RunDescription describeRun(
                unsigned long long maxSteps, std::optional<double> zStop) {
            opalx::run::TrackDescription track;
            track.line       = "L1";
            track.beamNames  = {"B1"};
            track.dT         = {dt};
            track.maxSteps   = {maxSteps};
            track.zStop      = {zStop.value_or(1e6)};
            track.zStopGiven = zStop.has_value();
            TrackRun runCommand;
            Attributes::setPredefinedString(*runCommand.findAttribute("METHOD"), "PARALLEL");
            Attributes::setString(*runCommand.findAttribute("FIELDSOLVER"), "FS1");
            return opalx::run::freeze(runCommand, track, 1);
        }

        /// Build the RUN and track it. Returns the setup's beam stops.
        std::vector<BeamStop> track(unsigned long long maxSteps, std::optional<double> zStop) {
            const opalx::run::RunDescription description = describeRun(maxSteps, zStop);
            EXPECT_TRUE(opalx::run::validate(description).empty());
            std::vector<std::string> errors;
            auto elements =
                    opalx::run::buildElements(description.line, description.context.p0, errors);
            EXPECT_TRUE(errors.empty());
            TrackingSetup setup =
                    opalx::run::buildTracking(description, std::move(elements), *objects_m);
            const std::vector<BeamStop> stops = setup.beamStops;
            ParallelTracker(std::move(setup)).execute();
            return stops;
        }

        /// Per particle, by ID, the lab position and momentum after tracking.
        std::vector<std::array<Vector, 2>> labPhaseSpace() const {
            const auto& pc = *objects_m->bunch->getParticleContainer(0);
            auto R  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc.R.getView());
            auto P  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc.P.getView());
            auto id = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc.ID.getView());
            // One create call numbers the particles of a rank consecutively, in file order.
            const size_t n     = pc.getLocalNum();
            std::int64_t first = std::numeric_limits<std::int64_t>::max();
            for (size_t i = 0; i < n; ++i)
                first = std::min(first, static_cast<std::int64_t>(id(i)));
            std::vector<std::array<Vector, 2>> result(n);
            for (size_t i = 0; i < n; ++i)
                result.at(static_cast<size_t>(id(i) - first)) = {
                        pc.getToLabTrafo().transformTo(R(i)), pc.getToLabTrafo().rotateTo(P(i))};
            return result;
        }

        static double beta(const Vector& p) {
            const double norm = std::sqrt(dot(p, p));
            return norm / std::sqrt(1.0 + norm * norm);
        }

        opalx::test::TemporaryDirectory directory_m;
        std::ostringstream output_m;
        std::unique_ptr<opalx::run::RunObjects> objects_m =
                std::make_unique<opalx::run::RunObjects>();
    };
}  // namespace

TEST_F(ParallelTrackerTest, DriftMovesEveryParticleByBetaCT) {
    track(10, std::nullopt);
    const PartBunch_t& bunch = *objects_m->bunch;
    const auto& pc           = *bunch.getParticleContainer(0);
    EXPECT_EQ(bunch.getGlobalTrackStep(), 10);
    EXPECT_NEAR(bunch.getT(), 10 * dt, 1e-23);
    ASSERT_EQ(pc.getTotalNum(), particles.size());

    const double t     = bunch.getT();
    const auto tracked = labPhaseSpace();
    for (size_t i = 0; i < particles.size(); ++i) {
        const Vector r0({particles[i][0], particles[i][1], particles[i][2]});
        const Vector p0({particles[i][3], particles[i][4], particles[i][5]});
        const Vector expected = r0 + Physics::c * t / std::sqrt(1.0 + dot(p0, p0)) * p0;
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_NEAR(tracked[i][0](d), expected(d), 1e-15) << i << " " << d;
            EXPECT_NEAR(tracked[i][1](d), p0(d), 1e-15) << i << " " << d;
        }
    }
    EXPECT_NEAR(pc.get_sPos(), beta(pc.getRefPartP()) * Physics::c * t, 1e-15);
}

TEST_F(ParallelTrackerTest, ZStopEndsTheRunOnceTheLastParticlePassed) {
    const double zStop = 0.05;
    const auto stops   = track(100, zStop);
    ASSERT_EQ(stops.size(), 1u);
    EXPECT_EQ(stops[0].zStop, zStop);

    const PartBunch_t& bunch = *objects_m->bunch;
    const auto& pc           = *bunch.getParticleContainer(0);
    EXPECT_LT(bunch.getGlobalTrackStep(), 100);
    ASSERT_EQ(pc.getTotalNum(), particles.size());

    // The reference moves along z: the plane is z = ZSTOP. Every particle is past it, the last
    // one passed it in the last step.
    double last = std::numeric_limits<double>::max(), step = 0.0;
    for (const auto& [r, p] : labPhaseSpace()) {
        last = std::min(last, r(2));
        step = std::max(step, Physics::c * dt * p(2) / std::sqrt(1.0 + dot(p, p)));
    }
    EXPECT_GE(last, zStop);
    EXPECT_LT(last, zStop + step);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ippl::initialize(argc, argv);
    const int result = RUN_ALL_TESTS();
    ippl::finalize();
    return result;
}
