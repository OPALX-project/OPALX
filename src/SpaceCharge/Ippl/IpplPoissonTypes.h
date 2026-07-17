/**
 * @file IpplPoissonTypes.h
 * @brief Typed host-side records used by the IPPL Poisson adapter.
 */

#ifndef OPALX_SPACE_CHARGE_IPPL_POISSON_TYPES_H
#define OPALX_SPACE_CHARGE_IPPL_POISSON_TYPES_H

#include <algorithm>
#include <optional>

#include "Ippl.h"
#include "Manager/datatypes.h"

namespace opalx::spacecharge {

    /** @brief Borrowed field bindings required by one IPPL Poisson backend. */
    struct IpplPoissonFields {
        Field_t<3>* chargeDensity          = nullptr;
        VField_t<double, 3>* electricField = nullptr;
        Field_t<3>* potential              = nullptr;
    };

    /** @brief Options for one backend solve, independent of a concrete IPPL solver class. */
    struct IpplPoissonSolveRequest {
        std::optional<ippl::Vector<double, 3>> greenFunctionShift;

        [[nodiscard]] bool hasShiftedGreenFunction() const {
            return greenFunctionShift.has_value();
        }
    };

    /** @brief Backend storage and normalization policies consumed by 3D PIC orchestration. */
    struct IpplPoissonCapabilities {
        bool isNoOp                              = false;
        bool supportsShiftedGreenFunction        = false;
        bool usesSeparatePotentialField          = false;
        bool requiresPotentialBoundaryConditions = false;
        bool normalizeChargeByCellVolume         = true;
        bool subtractNeutralizingBackground      = true;
        bool debugDumpChargeBeforeSolve          = false;
        bool debugDumpScalarAfterSolve           = false;
        bool debugDumpVectorAfterSolve           = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_IPPL_POISSON_TYPES_H
