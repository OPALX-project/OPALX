//
// Unit tests for the BinningCmd definition command.
//

#include <gtest/gtest.h>

#include "Ippl.h"

#include "Structure/BinningCmd.h"
#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"

#include <memory>
#include <string>

// Test-only helper exposing minimal setters for internal attributes.
class TestableBinningCmd : public BinningCmd {
public:
    void setMaxBins(int value) {
        Attributes::setReal(itsAttr[BINNING::MAXBINS], static_cast<double>(value));
    }

    void setDesiredWidth(double value) {
        Attributes::setReal(itsAttr[BINNING::DESIREDWIDTH], value);
    }

    void setBinningAlpha(double value) {
        Attributes::setReal(itsAttr[BINNING::BINNINGALPHA], value);
    }

    void setBinningBeta(double value) {
        Attributes::setReal(itsAttr[BINNING::BINNINGBETA], value);
    }

    void setParameterString(const std::string& value) {
        // PARAMETER is declared as a predefined string; this helper mirrors normal usage.
        Attributes::setPredefinedString(itsAttr[BINNING::PARAMETER], value);
    }
};

class BinningCmdTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() {
        ippl::finalize();
    }
};

// ConstructionDefaults:
// Verify that a freshly constructed exemplar has the documented defaults.
TEST_F(BinningCmdTest, ConstructionDefaults) {
    BinningCmd cmd;

    EXPECT_EQ(cmd.getMaxBins(), 128);
    EXPECT_DOUBLE_EQ(cmd.getDesiredWidth(), 0.1);
    EXPECT_DOUBLE_EQ(cmd.getBinningAlpha(), 1.0);
    EXPECT_DOUBLE_EQ(cmd.getBinningBeta(), 1.5);

    EXPECT_EQ(cmd.getParameter(), "VELOCITYZ");
    EXPECT_EQ(cmd.getParameterType(), BinningParameter::VELOCITYZ);
}

// GettersReflectAttributes:
// Changing the underlying attributes via the helper should be reflected by the getters.
TEST_F(BinningCmdTest, GettersReflectAttributes) {
    TestableBinningCmd cmd;

    cmd.setMaxBins(42);
    cmd.setDesiredWidth(0.25);
    cmd.setBinningAlpha(0.5);
    cmd.setBinningBeta(1.75);
    cmd.setParameterString("PZ");

    // update() refreshes the internal cached parameter string.
    cmd.update();

    EXPECT_EQ(cmd.getMaxBins(), 42);
    EXPECT_DOUBLE_EQ(cmd.getDesiredWidth(), 0.25);
    EXPECT_DOUBLE_EQ(cmd.getBinningAlpha(), 0.5);
    EXPECT_DOUBLE_EQ(cmd.getBinningBeta(), 1.75);

    EXPECT_EQ(cmd.getParameter(), "PZ");
}

// ExecuteMapsKnownParameters:
// For each supported PARAMETER option, execute() should set the correct enum.
TEST_F(BinningCmdTest, ExecuteMapsKnownParameters) {
    struct Case {
        const char* name;
        BinningParameter expected;
    } cases[] = {
        {"VELOCITYZ", BinningParameter::VELOCITYZ},
        {"POSITIONZ", BinningParameter::POSITIONZ},
        {"PZ",        BinningParameter::PZ},
        {"GAMMAZ",    BinningParameter::GAMMAZ}
    };

    for (const auto& c : cases) {
        TestableBinningCmd cmd;
        cmd.setParameterString(c.name);

        // execute() calls setParameterType() and then update().
        cmd.execute();

        EXPECT_EQ(cmd.getParameter(), std::string(c.name));
        EXPECT_EQ(cmd.getParameterType(), c.expected);
    }
}

// ExecuteThrowsOnUnknownParameter:
// An unsupported PARAMETER string should surface as an OpalException
// when setParameterType() (via execute()) cannot map it.
TEST_F(BinningCmdTest, ExecuteThrowsOnUnknownParameter) {
    TestableBinningCmd cmd;

    // This value is not present in the BinningParameter mapping table.
    cmd.setParameterString("FOO");

    EXPECT_THROW(
        cmd.execute(),
        OpalException
    );
}

// CloneCopiesState:
// clone(name) should copy attribute state and cached parameter type/name.
TEST_F(BinningCmdTest, CloneCopiesState) {
    TestableBinningCmd original;
    original.setMaxBins(7);
    original.setDesiredWidth(0.3);
    original.setBinningAlpha(0.8);
    original.setBinningBeta(1.2);
    original.setParameterString("POSITIONZ");

    original.execute(); // ensure parameterType_m is updated

    std::unique_ptr<BinningCmd> copy(original.clone("BINNING_COPY"));

    ASSERT_NE(copy, nullptr);
    EXPECT_EQ(copy->getOpalName(), "BINNING_COPY");

    EXPECT_EQ(copy->getMaxBins(), original.getMaxBins());
    EXPECT_DOUBLE_EQ(copy->getDesiredWidth(), original.getDesiredWidth());
    EXPECT_DOUBLE_EQ(copy->getBinningAlpha(), original.getBinningAlpha());
    EXPECT_DOUBLE_EQ(copy->getBinningBeta(), original.getBinningBeta());

    EXPECT_EQ(copy->getParameter(), original.getParameter());
    EXPECT_EQ(copy->getParameterType(), original.getParameterType());
}

// PrintInfoDoesNotThrow:
// Smoke test that printInfo can be called without throwing.
TEST_F(BinningCmdTest, PrintInfoDoesNotThrow) {
    BinningCmd cmd;

    // Inform is defined in Ippl; using a simple instance here is sufficient
    // to exercise the printing path.
    Inform os("BinningCmd::printInfo");

    EXPECT_NO_THROW({
        cmd.printInfo(os);
    });
}

