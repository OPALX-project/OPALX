#ifndef OPALX_FIELD_SOLVER_CONFIG_HPP
#define OPALX_FIELD_SOLVER_CONFIG_HPP

#include <string>

#include "Ippl.h"
#include "SpaceCharge/Solvers/SolverCapabilities.hpp"

namespace opalx {

/**
 * @brief Supported normalized binning parameters.
 */
enum class BinningParameterConfig : short {
    VelocityZ = 0,  //!< Bin by longitudinal velocity.
    PositionZ = 1,  //!< Reserved for parser compatibility; not supported by PartBunch yet.
    Pz = 2,         //!< Reserved for parser compatibility; not supported by PartBunch yet.
    GammaZ = 3      //!< Bin by longitudinal relativistic gamma.
};

/**
 * @brief Normalized BINNING command values used by the self-field runtime.
 *
 * This is a value snapshot of parser state. Runtime code should read this instead of
 * reaching back into `BinningCmd` for every solver step.
 */
struct BinningConfig {
    bool enabled             = false;
    std::string commandName  = "";
    BinningParameterConfig parameter = BinningParameterConfig::VelocityZ;
    std::string parameterName  = "VELOCITYZ";
    int maxBins              = 0;
    double desiredWidth      = 0.0;
    double alpha             = 0.0;
    double beta              = 0.0;
    bool adaptive            = true;
    bool dumpEnabled         = false;
    std::string dumpFileName = "NONE";
    int dumpFrequency        = 0;
    int tablePrintFrequency  = 0;

    /**
     * @brief Whether bin diagnostics should be written to a file.
     */
    bool dumpToFile() const { return enabled && dumpEnabled; }
};

/**
 * @brief Normalized FIELD_SOLVER command state plus backend capabilities.
 *
 * `FieldSolverCmd` remains the user-facing parser object. This config is the internal
 * handoff object used by `PartBunch`, field storage, and self-field drivers.
 */
template <unsigned Dim>
struct FieldSolverConfig {
    std::string commandName = "";
    std::string solverType  = "";
    ippl::Vector<int, Dim> meshCells;
    ippl::Vector<bool, Dim> domainDecomposition;
    double boundingBoxIncreasePercent = 0.0;
    SolverCapabilities capabilities;
    bool normalizeRhoByCellVolume       = true;
    bool subtractNeutralizingBackground = true;
    BinningConfig binning;

    /**
     * @brief Whether the selected backend is a no-op solver.
     */
    bool isNoOp() const { return capabilities.isNoOp; }

    /**
     * @brief Whether storage must allocate a separate potential field.
     */
    bool usesSeparatePotentialField() const {
        return capabilities.usesSeparatePotentialField;
    }
};

}  // namespace opalx

#endif  // OPALX_FIELD_SOLVER_CONFIG_HPP
