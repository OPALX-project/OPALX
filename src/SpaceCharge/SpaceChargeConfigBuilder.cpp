/**
 * @file SpaceChargeConfigBuilder.cpp
 * @brief Implements one-time parser-to-space-charge configuration conversion.
 */

#include "SpaceCharge/SpaceChargeConfigBuilder.h"

#include "AbstractObjects/OpalData.h"
#include "PartBunch/BCHandler.hpp"
#include "Structure/BinningCmd.h"
#include "Structure/EmissionSource.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include <algorithm>
#include <array>
#include <cmath>
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

        std::string resolveReferencePath(const ObjectValues& fieldSolver) {
            std::string path = FieldSolverCmd::getRefPathFileName(fieldSolver);
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
            std::array<FieldBoundaryCondition, 3> result;
            for (unsigned d = 0; d < 3; ++d) {
                switch (boundaryConditions[d]) {
                    case BCHandler<3>::OPEN:
                        result[d] = FieldBoundaryCondition::Open;
                        break;
                    case BCHandler<3>::PERIODIC:
                        result[d] = FieldBoundaryCondition::Periodic;
                        break;
                    case BCHandler<3>::DIRICHLET:
                        result[d] = FieldBoundaryCondition::Dirichlet;
                        break;
                }
            }
            return result;
        }

        std::size_t convertMeshSize(double value, const char* attribute) {
            if (!std::isfinite(value) || value <= 0.0 || std::trunc(value) != value
                || value > static_cast<double>(std::numeric_limits<int>::max())) {
                throw OpalException(
                        "SpaceChargeConfigBuilder::build",
                        std::string("FIELDSOLVER ") + attribute
                                + " must be a positive integer fitting an IPPL Index.");
            }
            return static_cast<std::size_t>(value);
        }

        bool usesLongitudinalResizeDecomposition(const std::vector<ObjectValues>& distributions) {
            for (const ObjectValues& distribution : distributions) {
                const std::string type = distribution.getString("TYPE");
                if (type == "FLATTOP" || type == "OPALFLATTOP") {
                    return true;
                }
            }
            return false;
        }

        std::optional<BinningConfig> buildBinningConfig(
                const std::optional<ObjectValues>& command) {
            if (!command) {
                return std::nullopt;
            }

            const int maximumBins = BinningCmd::getMaxBins(*command);
            if (maximumBins < 1) {
                throw OpalException(
                        "SpaceChargeConfigBuilder::build", "BINNING MAXBINS must be positive.");
            }

            BinningConfig values;
            values.name         = command->getName();
            values.maximumBins  = static_cast<std::size_t>(maximumBins);
            values.desiredWidth = BinningCmd::getDesiredWidth(*command);
            values.alpha        = BinningCmd::getBinningAlpha(*command);
            values.beta         = BinningCmd::getBinningBeta(*command);
            values.parameter    = convertBinningParameter(BinningCmd::getParameterType(*command));
            values.adaptive     = BinningCmd::getAdaptiveBinning(*command);
            values.tablePrintFrequency =
                    static_cast<std::size_t>(BinningCmd::getTablePrintFrequency(*command));

            if (BinningCmd::dumpBinsToFile(*command)) {
                values.dumpFile = BinningCmd::getDumpBinsFileName(*command);
                values.dumpFrequency =
                        static_cast<std::size_t>(BinningCmd::getDumpBinsFrequency(*command));
            } else {
                values.dumpFile.clear();
                values.dumpFrequency = 0;
            }
            return values;
        }

        DirichletPlaneConfig buildDirichletPlaneConfig(
                const std::vector<ObjectValues>& emissionSources) {
            bool enableImageCharge      = false;
            bool enableShiftedGreens    = false;
            double planeZ               = 0.0;
            int dumpFrequency           = 0;
            int maximumSteps            = 0;
            std::size_t zeroFaceSources = 0;
            std::size_t shiftedSources  = 0;

            for (const ObjectValues& source : emissionSources) {
                const bool zeroFace = EmissionSource::getZeroFaceR0Z(source);
                const bool shifted  = EmissionSource::getShiftedGreensFunction(source);
                const int sourceDumpFrequency =
                        EmissionSource::getZeroFacePlaneDumpFrequency(source);

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
                    planeZ            = EmissionSource::getR0(source)[2];
                    dumpFrequency     = sourceDumpFrequency;
                    maximumSteps      = EmissionSource::getZerofaceMaxSteps(source);
                } else {
                    ++shiftedSources;
                    enableShiftedGreens = true;
                    planeZ              = EmissionSource::getR0(source)[2];
                    if (sourceDumpFrequency > 0) {
                        throw OpalException(
                                "SpaceChargeConfigBuilder::build",
                                "ZEROFACEPLANEDUMP > 0 is not supported with "
                                "SHIFTED_GREENS_FUNCTION=true (the computational domain may "
                                "be far from R0Z, making the interpolated plane dump "
                                "meaningless).");
                    }
                    maximumSteps = EmissionSource::getZerofaceMaxSteps(source);
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
                        "Dirichlet boundary condition is only implemented for one plane.");
            }
            if (enableImageCharge && enableShiftedGreens) {
                throw OpalException(
                        "SpaceChargeConfigBuilder::build",
                        "Cannot have ZEROFACE_R0Z=true on one EMISSIONSOURCE and "
                        "SHIFTED_GREENS_FUNCTION=true on another; the two "
                        "Dirichlet-plane methods are mutually exclusive at the run level.");
            }
            if (dumpFrequency < 0 || maximumSteps < 0) {
                throw OpalException(
                        "buildSpaceChargeConfig",
                        "Dirichlet-plane dump frequency and maximum steps must not be negative.");
            }
            DirichletPlaneConfig values;
            values.kind               = enableImageCharge     ? DirichletPlaneType::ImageCharge
                                        : enableShiftedGreens ? DirichletPlaneType::ShiftedGreen
                                                              : DirichletPlaneType::None;
            values.planeZ             = planeZ;
            values.planeDumpFrequency = static_cast<std::size_t>(dumpFrequency);
            values.maximumSteps       = static_cast<std::size_t>(maximumSteps);
            return values;
        }

    }  // namespace

    SpaceChargeConfig buildSpaceChargeConfig(
            const ObjectValues& fieldSolver, const std::optional<ObjectValues>& binning,
            const std::vector<ObjectValues>& emissionSources,
            const std::vector<ObjectValues>& distributions) {
        const FieldSolverCmdType solverType = FieldSolverCmd::getFieldSolverCmdType(fieldSolver);
        if (solverType == FieldSolverCmdType::CG) {
            throw OpalException(
                    "SpaceChargeConfigBuilder::build",
                    "FIELDSOLVER TYPE=CG is recognized but not implemented.");
        }

        const std::array<std::size_t, 3> meshSize{
                convertMeshSize(FieldSolverCmd::getNX(fieldSolver), "NX"),
                convertMeshSize(FieldSolverCmd::getNY(fieldSolver), "NY"),
                convertMeshSize(FieldSolverCmd::getNZ(fieldSolver), "NZ")};
        const auto decomposition = FieldSolverCmd::getDomainDecomposition(fieldSolver);
        CartesianGridConfig grid;
        grid.meshSize                   = meshSize;
        grid.decomposition              = {decomposition[0], decomposition[1], decomposition[2]};
        grid.boundingBoxIncreasePercent = FieldSolverCmd::getBoxIncr(fieldSolver);

        if (FieldSolverCmd::hasFFT2D5Model(fieldSolver)
            && solverType != FieldSolverCmdType::FFT2D5) {
            throw OpalException("buildSpaceChargeConfig", "FFT2D5MODEL requires TYPE=FFT2D5.");
        }
        if (solverType == FieldSolverCmdType::FFT2D5) {
            if (FieldSolverCmd::getFFT2D5Model(fieldSolver) != "ORIGINAL") {
                throw OpalException(
                        "buildSpaceChargeConfig", "Only FFT2D5MODEL=ORIGINAL is supported.");
            }
            if (ippl::Comm->size() > 1 && !grid.decomposition[2]) {
                throw OpalException(
                        "buildSpaceChargeConfig", "Multiple FFT2D5 ranks require PARFFTZ=TRUE.");
            }
            const auto binsName = FieldSolverCmd::getBinsName(fieldSolver);
            if (!binsName.empty() && binsName != "NONE") {
                throw OpalException("buildSpaceChargeConfig", "FFT2D5 does not support BINS.");
            }
            const auto boundaries =
                    convertBoundaryConditions(FieldSolverCmd::constructBCHandler(fieldSolver));
            if (!std::all_of(boundaries.begin(), boundaries.end(), [](auto boundary) {
                    return boundary == FieldBoundaryCondition::Open;
                })) {
                throw OpalException(
                        "buildSpaceChargeConfig",
                        "FFT2D5 transverse Poisson slices are OPEN. Use CLOSEDRING for "
                        "longitudinal periodicity and PIPEMODE for the longitudinal-field model.");
            }
            const DirichletPlaneConfig dirichletPlane = buildDirichletPlaneConfig(emissionSources);
            if (dirichletPlane.enabled()) {
                throw OpalException(
                        "SpaceChargeConfigBuilder::build",
                        "FFT2D5 does not support Dirichlet planes.");
            }
            if (Options::loadBalancingPairWeight > 0.0) {
                throw OpalException(
                        "buildSpaceChargeConfig",
                        "LOADBALANCINGPAIRWEIGHT above 0 requires FIELDSOLVER TYPE=P3M.");
            }
            FFT2D5Config values;
            values.grid = grid;
            values.longitudinalFieldMode =
                    convertFFT2D5Mode(FieldSolverCmd::getPipeMode(fieldSolver));
            values.pipeSizeX             = FieldSolverCmd::getPipeSizeX(fieldSolver);
            values.pipeSizeY             = FieldSolverCmd::getPipeSizeY(fieldSolver);
            values.beamRadius            = FieldSolverCmd::getBeamRadius(fieldSolver);
            values.closedRing            = FieldSolverCmd::getClosedRing(fieldSolver);
            values.scatterLongitudinally = FieldSolverCmd::getScatterLongitudinally(fieldSolver);
            values.referencePathFile     = resolveReferencePath(fieldSolver);
            SpaceChargeConfig config     = std::move(values);
            validateSpaceChargeConfig(config);
            return config;
        }

        CartesianPIC3DConfig values;
        values.grid                       = grid;
        values.backend                    = convertPoissonSolverType(solverType);
        values.layoutRebuildDecomposition = usesLongitudinalResizeDecomposition(distributions)
                                                    ? std::array<bool, 3>{false, false, true}
                                                    : values.grid.decomposition;
        values.boundaryConditions =
                convertBoundaryConditions(FieldSolverCmd::constructBCHandler(fieldSolver));
        values.greenFunction = convertGreenFunction(FieldSolverCmd::getGreensFunction(fieldSolver));
        if (FieldSolverCmd::hasP3MCutoff(fieldSolver)
            && FieldSolverCmd::hasP3MCutoffCells(fieldSolver)) {
            throw OpalException(
                    "buildSpaceChargeConfig",
                    "FIELDSOLVER sets both RCUT and RCUT_CELLS. Set RCUT for a fixed P3M cutoff "
                    "in metres or RCUT_CELLS for a cutoff in units of the largest mesh spacing, "
                    "not both.");
        }
        values.p3mCutoff = FieldSolverCmd::getP3MCutoff(fieldSolver);
        // Without RCUT, P3M ties its cutoff to the mesh (RCUT_CELLS, default 2).
        const bool meshTiedCutoff = FieldSolverCmd::hasP3MCutoffCells(fieldSolver)
                                    || (values.backend == PoissonSolverType::P3M
                                        && !FieldSolverCmd::hasP3MCutoff(fieldSolver));
        values.p3mCutoffCells =
                meshTiedCutoff ? FieldSolverCmd::getP3MCutoffCells(fieldSolver) : 0.0;
        values.periodicBox = FieldSolverCmd::getPeriodicBox(fieldSolver);
        values.binning     = buildBinningConfig(binning);
        values.repartitionFrequency =
                Options::repartFreq > 0 ? static_cast<std::size_t>(Options::repartFreq) : 0;
        values.loadBalancingThreshold  = Options::loadBalancingThreshold;
        values.loadBalancingPairWeight = Options::loadBalancingPairWeight;
        values.dirichletPlane          = buildDirichletPlaneConfig(emissionSources);

        SpaceChargeConfig config = std::move(values);
        validateSpaceChargeConfig(config);
        return config;
    }

}  // namespace opalx::spacecharge
