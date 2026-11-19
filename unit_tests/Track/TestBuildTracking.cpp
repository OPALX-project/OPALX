//
// Test of buildTracking(): a RUN needs only its run description. The input objects are parsed
// and described with freeze(), then the parsed objects are deleted before the elements and the
// tracking setup are built.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#include <gtest/gtest.h>

#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

#include <unistd.h>

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Ippl.h"
#include "OpalConfigure/Configure.h"
#include "OpalParser/OpalParser.h"
#include "OpalParser/StringStream.h"
#include "PartBunch/PartBunch.h"
#include "Track/BuildElements.h"
#include "Track/BuildTracking.h"
#include "Track/FreezeRun.h"
#include "Track/TrackRun.h"
#include "Track/ValidateRun.h"
#include "Utility/Inform.h"

extern Inform* gmsg;

namespace {
    namespace fs = std::filesystem;

    // A drift and a quadrupole, a GAUSS beam of 10 particles and no field solver.
    const std::string definitions = R"(
        OPTION, VERSION = 10900;
        D1: DRIFT, L = 0.5, ELEMEDGE = 0;
        Q1: QUADRUPOLE, L = 0.1, K1 = 1.5, ELEMEDGE = 0.6;
        L1: LINE = (D1, Q1);
        DI: DISTRIBUTION, TYPE = GAUSS, SIGMAX = 1e-3, SIGMAY = 1e-3, SIGMAZ = 1e-3,
            NPARTDIST = 10;
        ES1: EMISSIONSOURCE, DISTRIBUTION = DI;
        SRC: EMISSIONSOURCELIST = (ES1);
        FS1: FIELDSOLVER, TYPE = NONE, NX = 8, NY = 8, NZ = 8;
        B1: BEAM, PARTICLE = PROTON, PC = 0.1, NALLOC = 10, BCHARGE = 1e-15, SOURCES = SRC;
    )";

    class BuildTrackingTest : public ::testing::Test {
    protected:
        void SetUp() override {
            gmsg = new Inform("", output_m, INFORM_ALL_NODES);
            Configure::configure();
            directory_m = fs::temp_directory_path()
                          / ("opalx-build-tracking-" + std::to_string(::getpid()));
            fs::create_directories(directory_m);
            OpalData::getInstance()->storeInputFn((directory_m / "test.in").string());
        }

        void TearDown() override {
            OpalData::deleteInstance();
            delete gmsg;
            gmsg = nullptr;
            fs::remove_all(directory_m);
        }

        /// Parse @p text as an input file.
        static void parse(const std::string& text) {
            std::string line = text;
            for (char& c : line) {
                if (c == '\n') c = ' ';
            }
            OpalParser parser;
            parser.run(new StringStream(line));
        }

        /**
         * Describe a RUN of L1 with beam B1 and field solver FS1, as TRACK and RUN would:
         * TRACK, LINE = L1, BEAM = B1, MAXSTEPS = 5, DT = 1e-12, ZSTOP = 1;
         * RUN, METHOD = "PARALLEL", FIELDSOLVER = FS1;
         */
        static opalx::run::RunDescription describeRun() {
            opalx::run::TrackDescription track;
            track.line       = "L1";
            track.beamNames  = {"B1"};
            track.dT         = {1e-12};
            track.maxSteps   = {5};
            track.zStop      = {1.0};
            track.zStopGiven = true;
            TrackRun runCommand;
            Attributes::setPredefinedString(*runCommand.findAttribute("METHOD"), "PARALLEL");
            Attributes::setString(*runCommand.findAttribute("FIELDSOLVER"), "FS1");
            return opalx::run::freeze(runCommand, track, 1);
        }

        fs::path directory_m;
        std::ostringstream output_m;
    };
}  // namespace

TEST_F(BuildTrackingTest, RunNeedsOnlyTheDescription) {
    parse(definitions);
    const opalx::run::RunDescription description = describeRun();
    ASSERT_TRUE(opalx::run::validate(description).empty());

    // Delete every parsed object; only the run description is left.
    OpalData::deleteInstance();
    OpalData::getInstance()->storeInputFn((directory_m / "test.in").string());
    ASSERT_EQ(OpalData::getInstance()->find("B1"), nullptr);
    ASSERT_EQ(OpalData::getInstance()->find("L1"), nullptr);

    std::vector<std::shared_ptr<ElementBase>> elements =
            opalx::run::buildElements(description.line, description.context.p0);
    opalx::run::RunObjects objects;
    TrackingSetup setup = opalx::run::buildTracking(description, std::move(elements), objects);

    // #S D1 Q1 #E were built; the beamline holds D1 and Q1, markers are not tracked.
    EXPECT_EQ(setup.beamline.getElements().size(), 2u);
    EXPECT_EQ(setup.dt, std::vector<double>{1e-12});
    EXPECT_EQ(setup.maxSteps, std::vector<unsigned long long>{5});
    EXPECT_EQ(setup.sStop, std::vector<double>{1.0});
    ASSERT_NE(setup.bunch, nullptr);
    EXPECT_EQ(setup.bunch->getParticleContainer()->getTotalNum(), 10u);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ippl::initialize(argc, argv);
    const int result = RUN_ALL_TESTS();
    ippl::finalize();
    return result;
}
