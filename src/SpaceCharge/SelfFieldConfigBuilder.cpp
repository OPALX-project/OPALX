/**
 * @file SelfFieldConfigBuilder.cpp
 * @brief Implements one-time parser-to-self-field configuration conversion.
 */

#include "SpaceCharge/SelfFieldConfigBuilder.h"

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

        PoissonBackendKind convertBackend(FieldSolverCmdType type) {
            switch (type) {
                case FieldSolverCmdType::NONE:
                    return PoissonBackendKind::None;
                case FieldSolverCmdType::FFT:
                    return PoissonBackendKind::FftPeriodic;
                case FieldSolverCmdType::OPEN:
                    return PoissonBackendKind::Open;
                case FieldSolverCmdType::CG:
                    return PoissonBackendKind::ConjugateGradient;
                case FieldSolverCmdType::P3M:
                    return PoissonBackendKind::P3M;
                case FieldSolverCmdType::FFT2D5:
                    break;
            }
            throw OpalException("SelfFieldConfigBuilder::build", "Unknown FIELDSOLVER TYPE value.");
        }

        std::string backendName(FieldSolverCmdType type) {
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

        GreenFunctionKind convertGreenFunction(const std::string& name) {
            if (name == "STANDARD") {
                return GreenFunctionKind::Standard;
            }
            if (name == "INTEGRATED") {
                return GreenFunctionKind::Integrated;
            }
            throw OpalException(
                    "SelfFieldConfigBuilder::build",
                    "Unknown FIELDSOLVER GREENSF value '" + name + "'.");
        }

        BinningParameterKind convertBinningParameter(BinningParameter parameter) {
            switch (parameter) {
                case BinningParameter::VELOCITYZ:
                    return BinningParameterKind::VelocityZ;
                case BinningParameter::POSITIONZ:
                    return BinningParameterKind::PositionZ;
                case BinningParameter::PZ:
                    return BinningParameterKind::MomentumZ;
                case BinningParameter::GAMMAZ:
                    return BinningParameterKind::GammaZ;
            }
            throw OpalException(
                    "SelfFieldConfigBuilder::build", "Unknown BINNING PARAMETER value.");
        }

        Pic2d5LongitudinalFieldMode convertPic2d5Mode(const std::string& name) {
            if (name == "OPEN") {
                return Pic2d5LongitudinalFieldMode::Open;
            }
            if (name == "CIRCULAR") {
                return Pic2d5LongitudinalFieldMode::Cylindrical;
            }
            if (name == "PLATES") {
                return Pic2d5LongitudinalFieldMode::Plates;
            }
            if (name == "NONE") {
                return Pic2d5LongitudinalFieldMode::None;
            }
            throw OpalException(
                    "SelfFieldConfigBuilder::build",
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

        std::array<BoundaryConditionKind, 3> convertBoundaryConditions(
                const BCHandler<3>& boundaryConditions) {
            BoundaryConditionKind kind = BoundaryConditionKind::Open;
            if (boundaryConditions.isAll(BCHandler<3>::PERIODIC)) {
                kind = BoundaryConditionKind::Periodic;
            } else if (boundaryConditions.isAll(BCHandler<3>::DIRICHLET)) {
                kind = BoundaryConditionKind::Dirichlet;
            } else if (!boundaryConditions.isAll(BCHandler<3>::OPEN)) {
                throw OpalException(
                        "SelfFieldConfigBuilder::build",
                        "Only uniform boundary conditions are currently supported.");
            }
            return {kind, kind, kind};
        }

        std::size_t convertMeshSize(double value, const char* attribute) {
            if (value <= 0.0
                || value > static_cast<double>(std::numeric_limits<std::size_t>::max())) {
                throw OpalException(
                        "SelfFieldConfigBuilder::build",
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
                        "SelfFieldConfigBuilder::build", "BINNING MAXBINS must be positive.");
            }

            BinningConfigValues values;
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
            return BinningConfig(std::move(values));
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
                                "SelfFieldConfigBuilder::build",
                                "ZEROFACE_R0Z and SHIFTED_GREENS_FUNCTION are mutually exclusive "
                                "on the same EMISSIONSOURCE. Enable exactly one.");
                    }

                    if (!zeroFace && !shifted) {
                        if (sourceDumpFrequency > 0) {
                            throw OpalException(
                                    "SelfFieldConfigBuilder::build",
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
                                    "SelfFieldConfigBuilder::build",
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
                        "SelfFieldConfigBuilder::build",
                        "Cannot have more than one emission source with ZEROFACE_R0Z=true, since "
                        "image charge computation is only implemented for one plane.");
            }
            if (shiftedSources > 1) {
                throw OpalException(
                        "SelfFieldConfigBuilder::build",
                        "Cannot have more than one emission source with "
                        "SHIFTED_GREENS_FUNCTION=true, since the shifted Green's function "
                        "correction is only implemented for one plane.");
            }
            if (enableImageCharge && enableShiftedGreens) {
                throw OpalException(
                        "SelfFieldConfigBuilder::build",
                        "Cannot have ZEROFACE_R0Z=true on one EMISSIONSOURCE and "
                        "SHIFTED_GREENS_FUNCTION=true on another; the two "
                        "Dirichlet-correction paths are mutually exclusive at the run level.");
            }
            if (enableShiftedGreens && solverType != FieldSolverCmdType::OPEN) {
                throw OpalException(
                        "SelfFieldConfigBuilder::build",
                        "SHIFTED_GREENS_FUNCTION=true requires FIELDSOLVER TYPE=OPEN (got '"
                                + backendName(solverType) + "').");
            }

            CorrectionConfigValues values;
            values.kind               = enableImageCharge     ? CorrectionKind::ImageCharge
                                        : enableShiftedGreens ? CorrectionKind::ShiftedGreen
                                                              : CorrectionKind::None;
            values.planeZ             = planeZ;
            values.planeDumpFrequency = static_cast<std::size_t>(dumpFrequency);
            values.maximumSteps       = static_cast<std::size_t>(maximumSteps);
            return CorrectionConfig(values);
        }

    }  // namespace

    SelfFieldConfig SelfFieldConfigBuilder::build(
            const FieldSolverCmd& fieldSolver,
            const std::vector<std::vector<EmissionSource*>>& emissionSources) {
        const FieldSolverCmdType solverType = fieldSolver.getFieldSolverCmdType();
        if (solverType == FieldSolverCmdType::CG) {
            throw OpalException(
                    "SelfFieldConfigBuilder::build",
                    "FIELDSOLVER TYPE=CG is recognized but not implemented.");
        }

        const std::array<std::size_t, 3> meshSize{
                convertMeshSize(fieldSolver.getNX(), "NX"),
                convertMeshSize(fieldSolver.getNY(), "NY"),
                convertMeshSize(fieldSolver.getNZ(), "NZ")};
        const auto decomposition = fieldSolver.getDomainDecomposition();
        ParticleStorageConfig3d storage;
        storage.meshSize                   = meshSize;
        storage.decomposition              = {decomposition[0], decomposition[1], decomposition[2]};
        storage.boundingBoxIncreasePercent = fieldSolver.getBoxIncr();

        if (solverType == FieldSolverCmdType::FFT2D5) {
            const CorrectionConfig correction = buildCorrectionConfig(emissionSources, solverType);
            if (correction.enabled()) {
                throw OpalException(
                        "SelfFieldConfigBuilder::build",
                        "FFT2D5 does not support source-plane corrections.");
            }
            Pic2d5ConfigValues values;
            values.meshSize              = meshSize;
            values.longitudinalFieldMode = convertPic2d5Mode(fieldSolver.getPipeMode());
            values.pipeSizeX             = fieldSolver.getPipeSizeX();
            values.pipeSizeY             = fieldSolver.getPipeSizeY();
            values.beamRadius            = fieldSolver.getBeamRadius();
            values.closedRing            = fieldSolver.getClosedRing();
            values.scatterLongitudinally = fieldSolver.getScatterLongitudinally();
            values.referencePathFile     = resolveReferencePath(fieldSolver);
            return SelfFieldConfig(Pic2d5Config(std::move(values)), std::move(storage));
        }

        Pic3DConfigValues values;
        values.backend            = convertBackend(solverType);
        values.meshSize           = meshSize;
        values.parallelDimensions = {decomposition[0], decomposition[1], decomposition[2]};
        values.layoutRebuildParallelDimensions =
                usesLongitudinalResizeDecomposition(emissionSources)
                        ? std::array<bool, 3>{false, false, true}
                        : values.parallelDimensions;
        values.boundaryConditions = convertBoundaryConditions(fieldSolver.constructBCHandler());
        values.greenFunction      = convertGreenFunction(fieldSolver.getGreensFunction());
        values.p3mCutoff = solverType == FieldSolverCmdType::P3M ? fieldSolver.getP3MCutoff() : 0.0;
        values.boundingBoxIncreasePercent = fieldSolver.getBoxIncr();
        values.binning                    = buildBinningConfig(fieldSolver.getBinningCmd());
        values.repartitionFrequency =
                Options::repartFreq > 0 ? static_cast<std::size_t>(Options::repartFreq) : 0;
        values.loadBalancingThreshold = Options::loadBalancingThreshold;
        values.correction             = buildCorrectionConfig(emissionSources, solverType);

        storage.periodicParticleBoundary = std::all_of(
                values.boundaryConditions.begin(), values.boundaryConditions.end(),
                [](BoundaryConditionKind kind) {
                    return kind == BoundaryConditionKind::Periodic;
                });
        if (solverType == FieldSolverCmdType::P3M) {
            storage.layoutKind    = ParticleLayoutKind::SpatialOverlap;
            storage.overlapCutoff = values.p3mCutoff;
        }

        return SelfFieldConfig(Pic3DConfig(std::move(values)), std::move(storage));
    }

}  // namespace opalx::spacecharge
