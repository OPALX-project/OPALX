//
// Tests for freeze() and the run description text: a full input is parsed with --check, so each
// RUN writes what it would track with and does not track.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include <unistd.h>

#include "AbsBeamline/EndFieldModel/EndFieldModelManager.h"
#include "AbstractObjects/OpalData.h"
#include "Ippl.h"
#include "OpalConfigure/Configure.h"
#include "OpalParser/OpalParser.h"
#include "OpalParser/StringStream.h"
#include "Utility/Inform.h"

extern Inform* gmsg;

namespace {
    namespace fs = std::filesystem;

    // A beam, a few elements and a line that uses -NAME, N*NAME and a nested LINE.
    const std::string common      = R"(
        OPTION, PSDUMPFREQ = 7, VERSION = 10900;
        REAL KV := 3;
        Q1: QUADRUPOLE, L = 0.1, K1 = 1.5, ELEMEDGE = 0;
        Q2: Q1, K1 = 2.5;
        Q3: QUADRUPOLE, L = 0.1, K1 := KV, ELEMEDGE = 0.8;
        D1: DRIFT, L = 0.5, ELEMEDGE = 0.2;
        SUB: LINE = (D1, Q3);
        L1: LINE = (Q1, SUB, -SUB, 2*Q2);
        B1: BEAM, PARTICLE = PROTON, PC = 0.1, NALLOC = 10, BCHARGE = 1e-15;
    )";
    const std::string trackAndRun = R"(
        TRACK, LINE = L1, BEAM = B1, MAXSTEPS = 5, DT = 1e-12, ZSTOP = 1;
        RUN, METHOD = "PARALLEL";
        ENDTRACK;
    )";

    class RunDescriptionTest : public ::testing::Test {
    protected:
        void SetUp() override {
            gmsg = new Inform("", output_m, INFORM_ALL_NODES);
            Configure::configure();
            directory_m = fs::temp_directory_path()
                          / ("opalx-run-description-" + std::to_string(::getpid()));
            fs::create_directories(directory_m);
            auto* opal = OpalData::getInstance();
            opal->storeInputFn((directory_m / "test.in").string());
            opal->checkOnly = true;
        }

        void TearDown() override {
            OpalData::deleteInstance();
            delete gmsg;
            gmsg = nullptr;
            fs::remove_all(directory_m);
        }

        void run(const std::string& text) {
            std::string line = text;
            for (char& c : line) {
                if (c == '\n') c = ' ';
            }
            OpalParser parser;
            parser.run(new StringStream(line));
        }

        /// Text of the dump the n-th RUN wrote.
        std::string dump(int runIndex) const {
            std::ifstream in(directory_m / ("test_run" + std::to_string(runIndex) + ".txt"));
            std::stringstream buffer;
            buffer << in.rdbuf();
            return buffer.str();
        }

        /// Text of section "[header]" up to the next section.
        static std::string section(const std::string& text, const std::string& header) {
            const auto begin = text.find("[" + header + "]");
            if (begin == std::string::npos) return "";
            const auto end = text.find("\n[", begin + 1);
            return text.substr(begin, end == std::string::npos ? std::string::npos : end - begin);
        }

        fs::path directory_m;
        std::ostringstream output_m;
    };
}  // namespace

TEST_F(RunDescriptionTest, LineOrderFollowsTheTracker) {
    run(common + trackAndRun);
    const std::string text = dump(1);
    ASSERT_FALSE(text.empty());
    // #S Q1 (#S D1 Q3 #E) (#S D1 Q3 #E) Q2 Q2 #E: each LINE brings its own markers, and the
    // reflected -SUB keeps its order but is flagged.
    EXPECT_NE(text.find("elements = 13"), std::string::npos) << section(text, "line L1");
    EXPECT_NE(text.find("[element 1 Q1]"), std::string::npos);
    EXPECT_NE(
            text.find("[element 3 D1]\nelementType = Drift\nreflected = false"), std::string::npos);
    EXPECT_NE(
            text.find("[element 7 D1]\nelementType = Drift\nreflected = true"), std::string::npos);
    EXPECT_NE(text.find("[element 10 Q2]"), std::string::npos);
    EXPECT_NE(text.find("[element 11 Q2]"), std::string::npos);
}

TEST_F(RunDescriptionTest, DerivedElementKeepsKeywordAndOwnValues) {
    run(common + trackAndRun);
    const std::string q2 = section(dump(1), "element 10 Q2");
    EXPECT_NE(q2.find("keyword = QUADRUPOLE"), std::string::npos) << q2;
    EXPECT_NE(q2.find("attr K1 (real, set) = 2.5"), std::string::npos) << q2;
    EXPECT_NE(q2.find("attr L (real, set) = 0.10000000000000001"), std::string::npos) << q2;
}

TEST_F(RunDescriptionTest, DeferredValueIsReadAtEachRun) {
    run(common + trackAndRun + "KV = 5;" + trackAndRun);
    EXPECT_NE(section(dump(1), "element 4 Q3").find("attr K1 (real, set) = 3"), std::string::npos);
    EXPECT_NE(section(dump(2), "element 4 Q3").find("attr K1 (real, set) = 5"), std::string::npos);
    EXPECT_NE(section(dump(2), "element 8 Q3").find("attr K1 (real, set) = 5"), std::string::npos);
}

TEST_F(RunDescriptionTest, OptionsTrackAndBeam) {
    run(common + trackAndRun);
    const std::string text = dump(1);
    EXPECT_NE(section(text, "options").find("psDumpFreq = 7"), std::string::npos);
    const std::string track = section(text, "track");
    EXPECT_NE(track.find("beams = B1"), std::string::npos) << track;
    EXPECT_NE(track.find("maxSteps = 5"), std::string::npos) << track;
    EXPECT_NE(track.find("dT = 1.0000000000000002e-12"), std::string::npos) << track;
    const std::string beam = section(text, "beam B1");
    EXPECT_NE(beam.find("keyword = BEAM"), std::string::npos) << beam;
    EXPECT_NE(beam.find("reference.momentum = 99999999.99999998"), std::string::npos) << beam;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ippl::initialize(argc, argv);
    const int result = RUN_ALL_TESTS();
    // Configure::configure() registers end-field models that hold Kokkos views; release them
    // before Kokkos is finalized, as main() does.
    endfieldmodel::EndFieldModelManager::clearEFMManager();
    ippl::finalize();
    return result;
}
