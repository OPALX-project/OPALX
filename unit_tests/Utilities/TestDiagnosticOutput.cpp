#include <gtest/gtest.h>

#include "Utilities/BeamBeamWindowAnimation.h"
#include "Utilities/DiagnosticOutput.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>

TEST(DiagnosticOutputTest, StepLimitIsInclusiveAndMissingOrNegativeDisablesOutput) {
    using opalx::diagnostics::spaceChargeFieldDumpRequested;
    EXPECT_FALSE(spaceChargeFieldDumpRequested(nullptr, 0));
    EXPECT_FALSE(spaceChargeFieldDumpRequested("", 0));
    EXPECT_FALSE(spaceChargeFieldDumpRequested("-1", 0));
    EXPECT_TRUE(spaceChargeFieldDumpRequested("0", 0));
    EXPECT_FALSE(spaceChargeFieldDumpRequested("0", 1));
    EXPECT_TRUE(spaceChargeFieldDumpRequested("12", 12));
    EXPECT_FALSE(spaceChargeFieldDumpRequested("12", 13));
    EXPECT_FALSE(spaceChargeFieldDumpRequested("12", -1));
    EXPECT_TRUE(spaceChargeFieldDumpRequested("9223372036854775807", 13));
}

TEST(DiagnosticOutputTest, InvalidStepLimitsNeverEnableUnboundedOutput) {
    for (const auto* value :
         {"yes", "12junk", "1.5", " 12", "12 ", "+12", "-", "9223372036854775808",
          "-9223372036854775809"}) {
        EXPECT_THROW(opalx::diagnostics::spaceChargeFieldDumpRequested(value, 0), OpalException)
                << value;
    }
}

TEST(DiagnosticOutputTest, AnimationCreatesParentAndReportsOpenFailure) {
    char pattern[]   = "animation-output-XXXXXX";
    const char* path = mkdtemp(pattern);
    ASSERT_NE(path, nullptr);
    const std::filesystem::path directory(path);
    struct Cleanup {
        std::filesystem::path path;
        ~Cleanup() {
            std::error_code error;
            std::filesystem::remove_all(path, error);
        }
    } cleanup{directory};
    const auto frames = directory / "nested" / "frames.txt";
    {
        BeamBeamWindowAnimation animation(frames.string());
        animation.render(0.0, -0.01, 0.01, -0.01, 0.01, 0.0, true, true);
    }
    std::ifstream input(frames);
    const std::string content{std::istreambuf_iterator<char>(input), {}};
    EXPECT_NE(content.find("FRAME 0"), std::string::npos);
    // A regular file cannot serve as the parent of another output file.
    EXPECT_THROW(BeamBeamWindowAnimation((frames / "child.txt").string()), OpalException);
}
