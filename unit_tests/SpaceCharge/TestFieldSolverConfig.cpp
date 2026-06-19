#include <gtest/gtest.h>

#include <string>

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "SpaceCharge/FieldSolverConfigBuilder.hpp"
#include "SpaceCharge/Solvers/PoissonBackendRegistry.hpp"
#include "Structure/BinningCmd.h"
#include "Structure/FieldSolverCmd.h"

namespace {

    class TestableFieldSolverCmd : public FieldSolverCmd {
    public:
        void setType(const std::string& value) {
            Attributes::setPredefinedString(itsAttr[FIELDSOLVER::TYPE], value);
        }

        void setBinsName(const std::string& value) {
            Attributes::setString(itsAttr[FIELDSOLVER::BINS], value);
        }

        void setBoxIncrease(double value) {
            Attributes::setReal(itsAttr[FIELDSOLVER::BBOXINCR], value);
        }

        void setDomainDecomposition(bool x, bool y, bool z) {
            Attributes::setBool(itsAttr[FIELDSOLVER::PARFFTX], x);
            Attributes::setBool(itsAttr[FIELDSOLVER::PARFFTY], y);
            Attributes::setBool(itsAttr[FIELDSOLVER::PARFFTZ], z);
        }
    };

    class TestableBinningCmd : public BinningCmd {
    public:
        void setMaxBins(int value) {
            Attributes::setReal(itsAttr[BINNING::MAXBINS], static_cast<double>(value));
        }

        void setDesiredWidth(double value) {
            Attributes::setReal(itsAttr[BINNING::DESIREDWIDTH], value);
        }

        void setAlpha(double value) {
            Attributes::setReal(itsAttr[BINNING::BINNINGALPHA], value);
        }

        void setBeta(double value) {
            Attributes::setReal(itsAttr[BINNING::BINNINGBETA], value);
        }

        void setAdaptiveBinning(bool value) {
            Attributes::setBool(itsAttr[BINNING::ADAPTIVEBINNING], value);
        }

        void setParameterString(const std::string& value) {
            Attributes::setPredefinedString(itsAttr[BINNING::PARAMETER], value);
        }

        void setDumpFileName(const std::string& value) {
            Attributes::setString(itsAttr[BINNING::DUMPBINSFILE], value);
        }

        void setDumpFrequency(int value) {
            Attributes::setReal(itsAttr[BINNING::DUMPBINSFREQ], static_cast<double>(value));
        }

        void setTablePrintFrequency(int value) {
            Attributes::setReal(itsAttr[BINNING::TABLEPRINTFREQ], static_cast<double>(value));
        }
    };

    void configureBaseFieldSolver(TestableFieldSolverCmd& command, const std::string& solverType) {
        command.setType(solverType);
        command.setNX(8);
        command.setNY(9);
        command.setNZ(10);
        command.setBoxIncrease(3.5);
        command.setDomainDecomposition(true, false, true);
        command.setBinsName("NONE");
    }

    void defineBinningCommand(
            const std::string& name, BinningParameter parameter, bool adaptiveBinning) {
        auto* command = new TestableBinningCmd();
        command->setOpalName(name);
        command->setMaxBins(11);
        command->setDesiredWidth(0.25);
        command->setAlpha(0.75);
        command->setBeta(1.25);
        command->setAdaptiveBinning(adaptiveBinning);
        command->setParameterString(parameter == BinningParameter::GAMMAZ ? "GAMMAZ" : "VELOCITYZ");
        command->setDumpFileName("stage6_bins");
        command->setDumpFrequency(5);
        command->setTablePrintFrequency(7);
        command->execute();
        OpalData::getInstance()->define(command);
    }

    using Registry = PoissonBackendRegistry<double, 3>;

    TEST(FieldSolverConfigTest, NormalizesNoOpSolverWithoutBinning) {
        TestableFieldSolverCmd command;
        configureBaseFieldSolver(command, "NONE");
        const auto config =
                opalx::makeFieldSolverConfig<3>(command, Registry::capabilitiesFor("NONE"));

        EXPECT_EQ(config.solverType, "NONE");
        EXPECT_TRUE(config.isNoOp());
        EXPECT_EQ(config.meshCells[0], 8);
        EXPECT_EQ(config.meshCells[1], 9);
        EXPECT_EQ(config.meshCells[2], 10);
        EXPECT_TRUE(config.domainDecomposition[0]);
        EXPECT_FALSE(config.domainDecomposition[1]);
        EXPECT_TRUE(config.domainDecomposition[2]);
        EXPECT_DOUBLE_EQ(config.boundingBoxIncreasePercent, 3.5);
        EXPECT_TRUE(config.normalizeRhoByCellVolume);
        EXPECT_TRUE(config.subtractNeutralizingBackground);
        EXPECT_FALSE(config.binning.enabled);
    }

    TEST(FieldSolverConfigTest, CopiesOpenSolverCapabilitiesForShiftedGreens) {
        TestableFieldSolverCmd command;
        configureBaseFieldSolver(command, "OPEN");
        const auto config =
                opalx::makeFieldSolverConfig<3>(command, Registry::capabilitiesFor("OPEN"));

        EXPECT_FALSE(config.isNoOp());
        EXPECT_TRUE(config.capabilities.supportsShiftedGreens);
        EXPECT_TRUE(config.capabilities.providesPotential);
        EXPECT_FALSE(config.usesSeparatePotentialField());
        EXPECT_FALSE(config.subtractNeutralizingBackground);
    }

    TEST(FieldSolverConfigTest, CopiesDensityPolicyFromCapabilities) {
        TestableFieldSolverCmd command;
        configureBaseFieldSolver(command, "OPEN");

        SolverCapabilities capabilities = Registry::capabilitiesFor("OPEN");
        capabilities.normalizeRhoByCellVolume       = false;
        capabilities.subtractNeutralizingBackground = true;

        const auto config = opalx::makeFieldSolverConfig<3>(command, capabilities);

        EXPECT_FALSE(config.normalizeRhoByCellVolume);
        EXPECT_TRUE(config.subtractNeutralizingBackground);
    }

    TEST(FieldSolverConfigTest, CapturesBinningAndDiagnosticSettings) {
        const std::string binningName = "STAGE6_CONFIG_BINNING";
        defineBinningCommand(binningName, BinningParameter::GAMMAZ, false);

        TestableFieldSolverCmd command;
        configureBaseFieldSolver(command, "OPEN");
        command.setBinsName(binningName);

        const auto config =
                opalx::makeFieldSolverConfig<3>(command, Registry::capabilitiesFor("OPEN"));

        ASSERT_TRUE(config.binning.enabled);
        EXPECT_EQ(config.binning.commandName, binningName);
        EXPECT_EQ(config.binning.parameter, opalx::BinningParameterConfig::GammaZ);
        EXPECT_EQ(config.binning.parameterName, "GAMMAZ");
        EXPECT_EQ(config.binning.maxBins, 11);
        EXPECT_DOUBLE_EQ(config.binning.desiredWidth, 0.25);
        EXPECT_DOUBLE_EQ(config.binning.alpha, 0.75);
        EXPECT_DOUBLE_EQ(config.binning.beta, 1.25);
        EXPECT_FALSE(config.binning.adaptive);
        EXPECT_EQ(config.binning.dumpFileName, "stage6_bins.json");
        EXPECT_TRUE(config.binning.dumpEnabled);
        EXPECT_EQ(config.binning.dumpFrequency, 5);
        EXPECT_EQ(config.binning.tablePrintFrequency, 7);
        EXPECT_TRUE(config.binning.dumpToFile());
    }

}  // namespace
