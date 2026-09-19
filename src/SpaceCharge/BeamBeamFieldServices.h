/** @file BeamBeamFieldServices.h
 * @brief Internal Cartesian field services borrowed by collective element interactions.
 */
#ifndef OPALX_SPACE_CHARGE_BEAM_BEAM_FIELD_SERVICES_H
#define OPALX_SPACE_CHARGE_BEAM_BEAM_FIELD_SERVICES_H

#include "PartBunch/CartesianDomain.h"
#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/SpaceChargeConfig.h"

#include <optional>
#include <string>
#include <vector>

namespace opalx::spacecharge {

    /**
     * @brief Sources represented by a BeamBeam solve in a fixed, symmetric z window.
     *
     * The copy has the primary's charge, reflected position and reversed longitudinal momentum.
     * These switches select field contributions, independently of whether the primary receives
     * a collective kick. Clearing the policy restores the ordinary space-charge algorithm.
     */
    struct BeamBeamSolvePolicy {
        bool primaryActive = true;
        bool copyActive    = false;
        /** @brief Retain the last bin's normalized primary rho before the in-place Poisson solve.
         */
        bool captureChargeDensity = false;
    };

    class BeamBeamFieldServices {
    public:
        using Domain            = CartesianDomain<double, 3>;
        using ParticleContainer = ::ParticleContainer<double, 3>;

        virtual ~BeamBeamFieldServices()                                        = default;
        virtual void configure(std::optional<BeamBeamSolvePolicy> policy)       = 0;
        [[nodiscard]] virtual const CartesianPIC3DConfig& configuration() const = 0;
        [[nodiscard]] virtual Domain& domain()                                  = 0;

        /**
         * @brief Replace target E/B from the most recent solve, in Cartesian solve axes.
         * Caller first transforms R and migrates target onto domain() on every MPI rank.
         */
        virtual void gatherFields(ParticleContainer& target) = 0;
        /** @brief Global represented charge [C], measured before normalization; absent otherwise.
         */
        [[nodiscard]] virtual std::optional<double> depositedCharge() const = 0;
        /**
         * @brief Dump captured normalized primary rho on one MPI rank; other layouts skip output.
         * Requires captureChargeDensity on the preceding solve. The copy is reconstructed in
         * field space and is not present in this primary-density diagnostic.
         */
        virtual void dumpChargeDensity(
                const std::string& prefix, const std::vector<std::string>& headers) = 0;
    };

}  // namespace opalx::spacecharge
#endif
