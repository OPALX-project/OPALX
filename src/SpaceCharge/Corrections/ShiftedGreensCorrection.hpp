#ifndef OPALX_SHIFTED_GREENS_CORRECTION_HPP
#define OPALX_SHIFTED_GREENS_CORRECTION_HPP

#include <cstddef>
#include <optional>

#include "SpaceCharge/Corrections/BoundaryCorrection.hpp"
#include "SpaceCharge/Solvers/SolverCapabilities.hpp"
#include "Utilities/OpalException.h"

/**
 * @brief Shifted Green's-function Dirichlet correction policy.
 *
 * The policy requires backend support for shifted Green's functions and describes the correction
 * as a normal primary-particle scatter with a shifted solve request. The accumulated field is
 * mirrored along the longitudinal axis by the caller.
 */
template <typename T, unsigned Dim>
class ShiftedGreensCorrection final : public BoundaryCorrection<T, Dim> {
public:
    /**
     * @brief Enable or disable the correction after validating backend capabilities.
     */
    void configure(bool enabled, double zPlane, const SolverCapabilities& capabilities) {
        if (enabled && !capabilities.supportsShiftedGreens) {
            throw OpalException(
                    "ShiftedGreensCorrection::configure",
                    "SHIFTED_GREENS_FUNCTION requires a field solver backend with shifted "
                    "Green's-function support.");
        }

        enabled_ = enabled;
        zPlane_  = zPlane;
    }

    const char* name() const override { return "shifted-greens"; }
    BoundaryCorrectionKind kind() const override { return BoundaryCorrectionKind::ShiftedGreens; }
    bool isEnabled() const override { return enabled_; }
    double planeZ() const override { return zPlane_; }

    std::optional<BoundaryCorrectionPass<T, Dim>> makePass(
            const BinContext<T, Dim>& context, const Vector_t<double, Dim>& meshOrigin,
            std::size_t longitudinalCellCount) const override {
        if (!enabled_) {
            return std::nullopt;
        }
        if (longitudinalCellCount == 0) {
            throw OpalException(
                    "ShiftedGreensCorrection::makePass",
                    "Cannot build shifted Green's-function request for an empty z domain.");
        }

        BoundaryCorrectionPass<T, Dim> pass;
        pass.scatterMode = ImageScatterMode::PrimaryOnly;
        pass.bFieldSign = -1.0;
        pass.flipAxis = static_cast<int>(Dim) - 1;
        pass.dumpDirichletPlane = false;
        pass.forceSkipFieldDump = false;

        ippl::Vector<T, Dim> shift(0.0);
        const double zCenterRest =
                meshOrigin[Dim - 1]
                + 0.5 * static_cast<double>(longitudinalCellCount)
                          * context.solveSpacing[Dim - 1];
        shift[Dim - 1] = static_cast<T>(2.0 * (zCenterRest - zPlane_));
        pass.solveRequest.greensShift = shift;

        return pass;
    }

private:
    bool enabled_ = false;
    double zPlane_ = 0.0;
};

#endif  // OPALX_SHIFTED_GREENS_CORRECTION_HPP
