#ifndef BUNCH_STATE_HANDLER_H
#define BUNCH_STATE_HANDLER_H

#include <array>
#include <memory>
#include <optional>
#include <vector>

/**
 * @class BunchStateHandler
 * @brief Centralised tracking of mutable bunch-level and per-container status flags.
 *
 * One instance is owned by PartBunch. Per-container coordinate and moment-cache flags live in
 * registered slots whose lifetime follows each ParticleContainer. Fixed Cartesian bounds are
 * bunch-wide and persist independently of container lifetimes.
 *
 * `unitlessPositions` describes the representation @f$R/(c\Delta t)@f$ and does not dirty moments.
 * `momentsDirty` is set after position or momentum changes and cleared after recomputation.
 * `Options::aggressiveStateSync` optionally converges these booleans across MPI ranks.
 */
class BunchStateHandler {
public:
    /**
     * @brief Persistent fixed bounds for Cartesian PIC solves.
     *
     * Bounds use metres in the Cartesian PIC solve frame. This state is bunch-wide rather than
     * tied to a particle container, and remains active until explicitly cleared.
     */
    struct FixedCartesianDomainState {
        std::array<double, 3> lower;
        std::array<double, 3> upper;

        bool operator==(const FixedCartesianDomainState&) const = default;
    };

    /**
     * @brief Per-container coordinate-representation and moment-cache state.
     */
    struct ContainerState {
        bool unitlessPositions = false;
        bool momentsDirty      = true;

        void setUnitlessPositions(bool v);
        void markMomentsDirty();
        void markMomentsClean();
    };

    BunchStateHandler() = default;

    /** @brief Allocate a slot whose strong ownership is transferred to the container. */
    std::shared_ptr<ContainerState> registerContainer();

    /**
     * @brief Activate fixed Cartesian PIC bounds.
     *
     * Reapplying identical bounds is idempotent. Active bounds must be cleared before they can be
     * replaced. This host-side operation is collective by contract: every MPI rank must call it
     * with identical values.
     *
     * @throws OpalException if a bound is non-finite, not strictly increasing, or attempts to
     * replace different active bounds.
     */
    void setFixedCartesianDomain(std::array<double, 3> lower, std::array<double, 3> upper);

    /**
     * @brief Return to domain-following Cartesian PIC behavior on the next solve.
     *
     * This host-side operation is collective by contract and must be called on every MPI rank.
     */
    void clearFixedCartesianDomain() noexcept { fixedCartesianDomain_m.reset(); }

    [[nodiscard]] const std::optional<FixedCartesianDomainState>& fixedCartesianDomain()
            const noexcept {
        return fixedCartesianDomain_m;
    }

private:
    std::vector<std::weak_ptr<ContainerState>> registered_m;
    std::optional<FixedCartesianDomainState> fixedCartesianDomain_m;
};

#endif
