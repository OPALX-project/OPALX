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
        /// Include the deposited physical primary's field contribution.
        bool primaryActive = true;
        /// Include its same-charge, longitudinally reflected virtual primary.
        bool copyActive = false;
        /** @brief Retain the last bin's normalized primary rho before the in-place Poisson solve.
         */
        bool captureChargeDensity = false;
    };

    /**
     * @brief Borrowed access to the Cartesian solver's BeamBeam state and fields.
     *
     * The solver owns this service and cached fields, and borrows the domain from
     * PartBunch. References must not outlive either owner. Only OPEN without a
     * Dirichlet plane is supported.
     * Configuration/domain changes do not solve fields; callers establish the
     * fixed domain, solve, then gather before the next solve overwrites the cache.
     */
    class BeamBeamFieldServices {
    public:
        /// Shared Cartesian mesh/decomposition type.
        using Domain = CartesianDomain<double, 3>;
        /// Particle storage accepted by the gather service.
        using ParticleContainer = ::ParticleContainer<double, 3>;

        virtual ~BeamBeamFieldServices() = default;
        /**
         * @brief Select sources for subsequent solves; nullopt restores ordinary space charge.
         * Call consistently on all ranks before solving. Disabling capture clears
         * the rho snapshot, but does not clear the previous charge measurement or
         * solve fields. Unsupported backends/Dirichlet planes throw OpalException.
         * @param policy Source/capture policy, or nullopt to disable BeamBeam mode.
         */
        virtual void configure(std::optional<BeamBeamSolvePolicy> policy) = 0;
        /// @return Borrowed immutable configuration, valid for the solver lifetime.
        [[nodiscard]] virtual const CartesianPIC3DConfig& configuration() const = 0;
        /**
         * @brief Borrow the mutable Cartesian solve domain (positions and spacing in metres).
         * Domain changes must agree on all ranks and precede particle migration
         * and a fresh solve. Changing it invalidates the previous gather geometry.
         * @return Domain borrowed from PartBunch, which must outlive the service.
         */
        [[nodiscard]] virtual Domain& domain() = 0;

        /**
         * @brief Replace target E/B from the most recent solve, in Cartesian solve axes.
         * Caller first transforms R and migrates target onto domain() on every MPI rank.
         * @param target Passive particles; E [V/m] and B [T] are overwritten, not added.
         * Requires a preceding successful BeamBeam solve in the same fixed domain.
         * Every rank participates, including ranks with no local target particles.
         */
        virtual void gatherFields(ParticleContainer& target) = 0;
        /**
         * @brief Signed global charge [C] represented by the last BeamBeam solve.
         * Sums deposited dt*Q weights divided by the solve time step before density
         * normalization, across bins, weighted by the enabled primary/copy count.
         * This includes the virtual copy, not witness charge. Read only after a
         * successful solve: a new solve resets the measurement (ordinary solves
         * leave nullopt); configure() alone leaves the last value unchanged.
         * @return Last measurement, or nullopt before a BeamBeam solve/after an ordinary solve.
         */
        [[nodiscard]] virtual std::optional<double> depositedCharge() const = 0;
        /**
         * @brief Dump captured normalized primary rho on one MPI rank; other layouts skip output.
         * Requires captureChargeDensity on the preceding solve. The copy is reconstructed in
         * field space and is not present in this primary-density diagnostic.
         * The captured last-bin rho uses backend normalization, not necessarily
         * physical C/m^3. Call before changing the domain or starting another solve.
         * @param prefix Diagnostic quantity/file suffix.
         * @param headers Metadata in key=value form, preserved in ASCII or parsed for HDF5.
         */
        virtual void dumpChargeDensity(
                const std::string& prefix, const std::vector<std::string>& headers) = 0;
    };

}  // namespace opalx::spacecharge
#endif
