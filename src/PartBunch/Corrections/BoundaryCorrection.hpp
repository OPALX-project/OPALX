#ifndef OPALX_BOUNDARY_CORRECTION_HPP
#define OPALX_BOUNDARY_CORRECTION_HPP

#include <cstddef>
#include <optional>

#include "PartBunch/BinContext.hpp"
#include "PartBunch/Solvers/SolveRequest.hpp"

/**
 * @brief Particle scatter source used for a boundary-correction solve pass.
 */
enum class ImageScatterMode { PrimaryAndImage, PrimaryOnly, ImageOnly };

/**
 * @brief Identifies the configured Dirichlet boundary-correction policy.
 */
enum class BoundaryCorrectionKind { ImageCharge, ShiftedGreens };

/**
 * @brief Executable description of one extra boundary-correction solve.
 *
 * The correction policy decides what should happen. `BinnedFieldSolver` still executes the
 * mechanics: scatter, solve, accumulate, optional field flip, and optional plane diagnostics.
 */
template <typename T, unsigned Dim>
struct BoundaryCorrectionPass {
    ImageScatterMode scatterMode = ImageScatterMode::PrimaryOnly;
    SolveRequest<T, Dim> solveRequest;
    double bFieldSign = -1.0;
    int flipAxis = -1;
    bool dumpDirichletPlane = false;
    bool forceSkipFieldDump = true;
};

/**
 * @brief Interface for optional Dirichlet boundary corrections.
 *
 * Implementations own correction-specific configuration and translate it into a
 * `BoundaryCorrectionPass`. They do not scatter particles or touch fields directly, which keeps
 * the field-solver loop independent from the concrete correction policy.
 */
template <typename T, unsigned Dim>
class BoundaryCorrection {
public:
    virtual ~BoundaryCorrection() = default;

    virtual const char* name() const = 0;
    virtual BoundaryCorrectionKind kind() const = 0;
    virtual bool isEnabled() const = 0;
    virtual double planeZ() const = 0;

    /**
     * @brief Return whether the correction is enabled and within the configured step budget.
     */
    bool isActiveForStep(const std::size_t step, const int maxSteps) const {
        if (!isEnabled()) {
            return false;
        }
        if (maxSteps <= 0) {
            return true;
        }
        return step < static_cast<std::size_t>(maxSteps);
    }

    virtual std::optional<BoundaryCorrectionPass<T, Dim>> makePass(
            const BinContext<T, Dim>& context, const Vector_t<double, Dim>& meshOrigin,
            std::size_t longitudinalCellCount) const = 0;
};

#endif  // OPALX_BOUNDARY_CORRECTION_HPP
