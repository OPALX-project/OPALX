#ifndef OPALX_FIELD_SOLVER_CONFIG_BUILDER_HPP
#define OPALX_FIELD_SOLVER_CONFIG_BUILDER_HPP

#include "PartBunch/FieldSolverConfig.hpp"
#include "Structure/BinningCmd.h"
#include "Structure/FieldSolverCmd.h"

namespace opalx {

/**
 * @brief Convert parser-level binning enum values into runtime config enum values.
 *
 * This keeps `FieldSolverConfig` independent from parser headers while preserving
 * every parser spelling for validation and diagnostics.
 */
inline BinningParameterConfig makeBinningParameterConfig(const BinningParameter parameter) {
    switch (parameter) {
        case BinningParameter::VELOCITYZ:
            return BinningParameterConfig::VelocityZ;
        case BinningParameter::POSITIONZ:
            return BinningParameterConfig::PositionZ;
        case BinningParameter::PZ:
            return BinningParameterConfig::Pz;
        case BinningParameter::GAMMAZ:
            return BinningParameterConfig::GammaZ;
    }

    return BinningParameterConfig::VelocityZ;
}

/**
 * @brief Build the normalized field-solver config from OPAL parser commands.
 *
 * @param command      Parsed FIELD_SOLVER command.
 * @param capabilities Backend capabilities from `PoissonBackendRegistry`.
 */
template <unsigned Dim>
FieldSolverConfig<Dim> makeFieldSolverConfig(
        FieldSolverCmd& command, const SolverCapabilities& capabilities) {
    static_assert(Dim == 3, "FieldSolverConfig currently normalizes 3D OPAL field solvers.");

    FieldSolverConfig<Dim> config;
    config.commandName = command.getOpalName();
    config.solverType  = command.getType();
    config.meshCells   = ippl::Vector<int, Dim>(
            static_cast<int>(command.getNX()), static_cast<int>(command.getNY()),
            static_cast<int>(command.getNZ()));

    const ippl::Vector<bool, 3> domDec = command.getDomDec();
    for (unsigned d = 0; d < Dim; ++d) {
        config.domainDecomposition[d] = domDec[d];
    }

    config.boundingBoxIncreasePercent = command.getBoxIncr();
    config.capabilities               = capabilities;

    // Copy density policies from backend capabilities once at setup time.
    config.normalizeRhoByCellVolume       = capabilities.normalizeRhoByCellVolume;
    config.subtractNeutralizingBackground = capabilities.subtractNeutralizingBackground;

    BinningCmd* binningCmd = command.getBinningCmd();
    if (!binningCmd) {
        return config;
    }

    config.binning.enabled             = true;
    config.binning.commandName         = binningCmd->getOpalName();
    config.binning.parameter           = makeBinningParameterConfig(binningCmd->getParameterType());
    config.binning.parameterName       = binningCmd->getParameter();
    config.binning.maxBins             = binningCmd->getMaxBins();
    config.binning.desiredWidth        = binningCmd->getDesiredWidth();
    config.binning.alpha               = binningCmd->getBinningAlpha();
    config.binning.beta                = binningCmd->getBinningBeta();
    config.binning.adaptive            = binningCmd->getAdaptiveBinning();
    config.binning.dumpEnabled         = binningCmd->dumpBinsToFile();
    config.binning.dumpFileName        = binningCmd->getDumpBinsFileName();
    config.binning.dumpFrequency       = binningCmd->getDumpBinsFrequency();
    config.binning.tablePrintFrequency = binningCmd->getTablePrintFrequency();

    return config;
}

}  // namespace opalx

#endif  // OPALX_FIELD_SOLVER_CONFIG_BUILDER_HPP
