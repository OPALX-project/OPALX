/**
 * @file SpaceChargeConfigBuilder.cpp
 * @brief Implements one-time parser-to-space-charge configuration conversion.
 */

#include "SpaceCharge/SpaceChargeConfigBuilder.h"

#include "AbstractObjects/OpalData.h"
#include "Distribution/Distribution.h"
#include "PartBunch/BCHandler.hpp"
#include "Structure/BinningCmd.h"
#include "Structure/EmissionSource.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include <algorithm>
#include <array>
#include <limits>
#include <string>

namespace opalx::spacecharge {
    namespace {

        PoissonSolverType convertPoissonSolverType(FieldSolverCmdType type) {
            switch (type) {
                case FieldSolverCmdType::NONE:
                    return PoissonSolverType::None;
                case FieldSolverCmdType::FFT:
                    return PoissonSolverType::PeriodicFFT;
                case FieldSolverCmdType::OPEN:
                    return PoissonSolverType::Open;
                case FieldSolverCmdType::CG:
                    return PoissonSolverType::ConjugateGradient;
                case FieldSolverCmdType::P3M:
                    return PoissonSolverType::P3M;
                case FieldSolverCmdType::FFT2D5:
                    break;
            }
            throw OpalException(
                    "SpaceChargeConfigBuilder::build", "Unknown FIELDSOLVER TYPE value.");
        }

        std::string fieldSolverTypeName(FieldSolverCmdType type) {
            switch (type) {
                case FieldSolverCmdType::NONE:
                    return "NONE";
                case FieldSolverCmdType::FFT:
                    return "FFT";
                case FieldSolverCmdType::OPEN:
                    return "OPEN";
                case FieldSolverCmdType::CG:
                    return "CG";
                case FieldSolverCmdType::P3M:
                    return "P3M";
                case FieldSolverCmdType::FFT2D5:
                    return "FFT2D5";
            }
            return "(unknown)";
        }

        GreenFunctionType convertGreenFunction(const std::string& name) {
            if (name == "STANDARD") {
                return GreenFunctionType::Standard;
            }
            if (name == "INTEGRATED") {
                return GreenFunctionType::Integrated;
            }
            throw OpalException(
                    "SpaceChargeConfigBuilder::build",
                    "Unknown FIELDSOLVER GREENSF value '" + name + "'.");
        }

        BinningVariable convertBinningParameter(BinningParameter parameter) {
            switch (parameter) {
                case BinningParameter::VELOCITYZ:
                    return BinningVariable::VelocityZ;
                case BinningParameter::POSITIONZ:
                    return BinningVariable::PositionZ;
                case BinningParameter::PZ:
                    return BinningVariable::MomentumZ;
                case BinningParameter::GAMMAZ:
                    return BinningVariable::GammaZ;
            }
            throw OpalException(
                    "SpaceChargeConfigBuilder::build", "Unknown BINNING PARAMETER value.");
        }

        FFT2D5LongitudinalFieldMode convertFFT2D5Mode(const std::string& name) {
            if (name == "OPEN") {
                return FFT2D5LongitudinalFieldMode::Open;
            }
            if (name == "CIRCULAR") {
                return FFT2D5LongitudinalFieldMode::Cylindrical;
            }
            if (name == "PLATES") {
                return FFT2D5LongitudinalFieldMode::Plates;
            }
            if (name == "NONE") {
                return FFT2D5LongitudinalFieldMode::None;
            }
            throw OpalException(
                    "SpaceChargeConfigBuilder::build",
                    "Unknown FIELDSOLVER PIPEMODE value '" + name + "'.");
        }

        std::string resolveReferencePath(const FieldSolverCmd& fieldSolver) {
            std::string path = fieldSolver.getRefPathFileName();
            if (!path.empty()) {
                return path;
            }
            OpalData* opal = OpalData::getInstance();
            return Util::combineFilePath(
                    {opal->getAuxiliaryOutputDirectory(),
                     opal->getInputBasename() + "_DesignPath.dat"});
        }

        std::array<FieldBoundaryCondition, 3> convertBoundaryConditions(
                const BCHandler<3>& boundaryConditions) {
            FieldBoundaryCondition kind = FieldBoundaryCondition::Open;
            if (boundaryConditions.isAll(BCHandler<3>::PERIODIC)) {
                kind = FieldBoundaryCondition::Periodic;
            } else if (boundaryConditions.isAll(BCHandler<3>::DIRICHLET)) {
                kind = FieldBoundaryCondition::Dirichlet;
            } else if (!boundaryConditions.isAll(BCHandler<3>::OPEN)) {
                throw OpalException(
                        "SpaceChargeConfigBuilder::build",
                        "Only uniform boundary conditions are currently supported.");
            }
            return {kind, kind, kind};
        }

        std::size_t convertMeshSize(double value, const char* attribute) {
            if (value <= 0.0
                || value > static_cast<double>(std::numeric_limits<std::size_t>::max())) {
                throw OpalException(
                        "SpaceChargeConfigBuilder::build",
                        std::string("FIELDSOLVER ") + attribute + " must be positive.");
            }
            return static_cast<std::size_t>(value);
        }

        bool usesLongitudinalResizeDecomposition(
                const std::vector<std::vector<EmissionSource*>>& emissionSources) {
            for (const auto& sourceList : emissionSources) {
                for (const EmissionSource* source : sourceList) {
                    if (source == nullptr) {
                        continue;
                    }
                    const DistributionType type =
                            Distribution::find(source->getDistributionName())->getType();
                    if (type == DistributionType::FLATTOP
                        || type == DistributionType::OPALFLATTOP) {
                        return true;
                    }
                }
            }
            return false;
        }

        std::optional<BinningConfig> buildBinningConfig(const BinningCmd* command) {
            if (command == nullptr) {
                return std::nullopt;
            }

            const int maximumBins = command->getMaxBins();
            if (maximumBins < 1) {
                throw OpalException(
                        "SpaceChargeConfigBuilder::build", "BINNING MAXBINS must be positive.");
            }

            BinningConfig values;
            values.name         = command->getOpalName();
            values.maximumBins  = static_cast<std::size_t>(maximumBins);
            values.desiredWidth = command->getDesiredWidth();
            values.alpha        = command->getBinningAlpha();
            values.beta         = command->getBinningBeta();
            values.parameter    = convertBinningParameter(command->getParameterType());
            values.adaptive     = command->getAdaptiveBinning();
            values.tablePrintFrequency =
                    static_cast<std::size_t>(command->getTablePrintFrequency());

            if (command->dumpBinsToFile()) {
                values.dumpFile      = command->getDumpBinsFileName();
                values.dumpFrequency = static_cast<std::size_t>(command->getDumpBinsFrequency());
            } else {
                values.dumpFile.clear();
                values.dumpFrequency = 0;
            }
            return values;
        }

        CorrectionConfig buildCorrectionConfig(
                const std::vector<std::vector<EmissionSource*>>& emissionSources,
                FieldSolverCmdType solverType) {
            bool enableImageCharge      = false;
            bool enableShiftedGreens    = false;
            double planeZ               = 0.0;
            int dumpFrequency           = 0;
            int maximumSteps            = 0;
            std::size_t zeroFaceSources = 0;
            std::size_t shiftedSources  = 0;

            for (const auto& sourceList : emissionSources) {
                for (const EmissionSource* source : sourceList) {
                    if (source == nullptr) {
                        continue;
                    }

                    const bool zeroFace           = source->getZeroFaceR0Z();
                    const bool shifted            = source->getShiftedGreensFunction();
                    const int sourceDumpFrequency = source->getZeroFacePlaneDumpFrequency();

                    if (zeroFace && shifted) {
                        throw OpalException(
                                "SpaceChargeConfigBuilder::build",
                                "ZEROFACE_R0Z and SHIFTED_GREENS_FUNCTION are mutually exclusive "
                                "on the same EMISSIONSOURCE. Enable exactly one.");
                    }

                    if (!zeroFace && !shifted) {
                        if (sourceDumpFrequency > 0) {
                            throw OpalException(
                                    "SpaceChargeConfigBuilder::build",
                                    "ZEROFACEPLANEDUMP > 0 requires ZEROFACE_R0Z=true on the "
                                    "same EMISSIONSOURCE. (Dumping is not supported for "
                                    "SHIFTED_GREENS_FUNCTION since the computational domain may "
                                    "be far from R0Z.)");
                        }
                        continue;
                    }

                    if (zeroFace) {
                        ++zeroFaceSources;
                        enableImageCharge = true;
                        planeZ            = source->getR0()[2];
                        dumpFrequency     = sourceDumpFrequency;
                        maximumSteps      = source->getZerofaceMaxSteps();
                    } else {
                        ++shiftedSources;
                        enableShiftedGreens = true;
                        planeZ              = source->getR0()[2];
                        if (sourceDumpFrequency > 0) {
                            throw OpalException(
                                    "SpaceChargeConfigBuilder::build",
                                    "ZEROFACEPLANEDUMP > 0 is not supported with "
                                    "SHIFTED_GREENS_FUNCTION=true (the computational domain may "
                                    "be far from R0Z, making the interpolated plane dump "
                                    "meaningless).");
                        }
                        maximumSteps = source->getZerofaceMaxSteps();
                    }
                }
            }

            if (zeroFaceSources > 1) {
                throw OpalException(
                        "SpaceChargeConfigBuilder::build",
                        "Cannot have more than one emission source with ZEROFACE_R0Z=true, since "
                        "image charge computation is only implemented for one plane.");
            }
            if (shiftedSources > 1) {
                throw OpalException(
                        "SpaceChargeConfigBuilder::build",
                        "Cannot have more than one emission source with "
                        "SHIFTED_GREENS_FUNCTION=true, since the shifted Green's function "
                        "correction is only implemented for one plane.");
            }
            if (enableImageCharge && enableShiftedGreens) {
                throw OpalException(
                        "SpaceChargeConfigBuilder::build",
                        "Cannot have ZEROFACE_R0Z=true on one EMISSIONSOURCE and "
                        "SHIFTED_GREENS_FUNCTION=true on another; the two "
                        "Dirichlet-correction paths are mutually exclusive at the run level.");
            }
            if (enableShiftedGreens && solverType != FieldSolverCmdType::OPEN) {
                throw OpalException(
                        "SpaceChargeConfigBuilder::build",
                        "SHIFTED_GREENS_FUNCTION=true requires FIELDSOLVER TYPE=OPEN (got '"
                                + fieldSolverTypeName(solverType) + "').");
            }

            CorrectionConfig values;
            values.kind               = enableImageCharge ? SpaceChargeCorrectionType::ImageCharge
                                        : enableShiftedGreens ? SpaceChargeCorrectionType::ShiftedGreen
                                                              : SpaceChargeCorrectionType::None;
            values.planeZ             = planeZ;
            values.planeDumpFrequency = static_cast<std::size_t>(dumpFrequency);
            values.maximumSteps       = static_cast<std::size_t>(maximumSteps);
            return values;
        }

    }  // namespace

    SpaceChargeConfig buildSpaceChargeConfig(
            const FieldSolverCmd& fieldSolver,
            const std::vector<std::vector<EmissionSource*>>& emissionSources) {
        const FieldSolverCmdType solverType = fieldSolver.getFieldSolverCmdType();
        if (solverType == FieldSolverCmdType::CG) {
            throw OpalException(
                    "SpaceChargeConfigBuilder::build",
                    "FIELDSOLVER TYPE=CG is recognized but not implemented.");
        }

        const std::array<std::size_t, 3> meshSize{
                convertMeshSize(fieldSolver.getNX(), "NX"),
                convertMeshSize(fieldSolver.getNY(), "NY"),
                convertMeshSize(fieldSolver.getNZ(), "NZ")};
        const auto decomposition = fieldSolver.getDomainDecomposition();
        CartesianGridConfig grid;
        grid.meshSize                   = meshSize;
        grid.decomposition              = {decomposition[0], decomposition[1], decomposition[2]};
        grid.boundingBoxIncreasePercent = fieldSolver.getBoxIncr();

        if (solverType == FieldSolverCmdType::FFT2D5) {
            const CorrectionConfig correction = buildCorrectionConfig(emissionSources, solverType);
            if (correction.enabled()) {
                throw OpalException(
                        "SpaceChargeConfigBuilder::build",
                        "FFT2D5 does not support source-plane corrections.");
            }
            FFT2D5Config values;
            values.grid                  = grid;
            values.longitudinalFieldMode = convertFFT2D5Mode(fieldSolver.getPipeMode());
            values.pipeSizeX             = fieldSolver.getPipeSizeX();
            values.pipeSizeY             = fieldSolver.getPipeSizeY();
            values.beamRadius            = fieldSolver.getBeamRadius();
            values.closedRing            = fieldSolver.getClosedRing();
            values.scatterLongitudinally = fieldSolver.getScatterLongitudinally();
            values.referencePathFile     = resolveReferencePath(fieldSolver);
            SpaceChargeConfig config     = std::move(values);
            validateSpaceChargeConfig(config);
            return config;
        }

        CartesianPICConfig values;
        values.grid                       = grid;
        values.backend                    = convertPoissonSolverType(solverType);
        values.layoutRebuildDecomposition = usesLongitudinalResizeDecomposition(emissionSources)
                                                    ? std::array<bool, 3>{false, false, true}
                                                    : values.grid.decomposition;
        values.boundaryConditions = convertBoundaryConditions(fieldSolver.constructBCHandler());
        values.greenFunction      = convertGreenFunction(fieldSolver.getGreensFunction());
        values.p3mCutoff = solverType == FieldSolverCmdType::P3M ? fieldSolver.getP3MCutoff() : 0.0;
        values.binning   = buildBinningConfig(fieldSolver.getBinningCmd());
        values.repartitionFrequency =
                Options::repartFreq > 0 ? static_cast<std::size_t>(Options::repartFreq) : 0;
        values.loadBalancingThreshold = Options::loadBalancingThreshold;
        values.correction             = buildCorrectionConfig(emissionSources, solverType);

        SpaceChargeConfig config = std::move(values);
        validateSpaceChargeConfig(config);
        return config;
    }

}  // namespace opalx::spacecharge
