#ifndef BUNCH_STATE_HANDLER_H
#define BUNCH_STATE_HANDLER_H

#include <memory>
#include <vector>

/**
 * @class BunchStateHandler
 * @brief Centralised tracking of mutable bunch-level and per-container status flags.
 *
 * Single instance per PartBunch, shared with every component that needs
 * access (ParticleContainer, DistributionMoments, ...).
 *
 * Two classes of state live here:
 *
 * - **Bunch-wide flags** (`firstRepartition`, future `emittingNow`) apply to the
 *   bunch as a whole and are stored as plain members.
 *
 * - **Per-container flags** (`momentsDirty`, `unitlessPositions`) live in the
 *   nested `ContainerState` struct. Each `ParticleContainer` registers itself
 *   via `registerContainer()` at `setBunchStateHandler` time and receives a
 *   `std::shared_ptr<ContainerState>` slot. The container holds the only strong
 *   reference; the handler keeps a `std::weak_ptr` so the slot is released
 *   automatically when the container is destroyed (no unregister needed).
 *
 * ### Invariants
 *
 * - **unitlessPositions**: true while the container's positions are in the
 *   dimensionless form R' = R / (c * dt). Toggled only by
 *   ParticleContainer::switchToUnitlessPositions / switchOffUnitlessPositions.
 *   Does NOT mark the container's moments dirty (coordinate representation
 *   change only).
 *
 * - **momentsDirty**: set whenever a physics operation mutates this
 *   container's particle positions (R) or momenta (P) -- push, kick,
 *   emission, particle deletion. Cleared by DistributionMoments::computeMoments
 *   once the moments cache is consistent with the particle state.
 *
 * - **firstRepartition**: true until the first ORB-style repartition has run
 *   for the bunch. Bunch-wide because the load balancer is shared.
 *
 * ### MPI consistency
 *
 * Every flag is conceptually consistent across MPI ranks. The OPALX option
 * `AGGRESSIVE_STATE_SYNC` (see `Options::aggressiveStateSync`) forces every
 * setter below to perform an `ippl::Comm->allreduce` with
 * `std::logical_or<bool>` so ranks converge to the same (conservative) value
 * even if a caller mutated the flag on only a subset of ranks. Default off:
 * the allreduce on every mutation adds up, and correctly-written callers
 * already set the same value on every rank.
 */
class BunchStateHandler {
public:
    /**
     * @brief Per-container slot of mutable flags. One per registered
     *        ParticleContainer. Lifetime is tied to the owning container
     *        (see registerContainer() for details).
     *
     * Setters here own the MPI consistency logic directly; there is no need
     * to route back through the handler, because `Options::aggressiveStateSync`
     * is a namespace flag and `ippl::Comm` is a singleton - i.e. the allreduce
     * helper is stateless, so the struct can call it without holding any
     * handler reference.
     */
    struct ContainerState {
        bool unitlessPositions = false;
        bool momentsDirty      = true;

        void setUnitlessPositions(bool v);
        void markMomentsDirty();
        void clearMomentsDirty();
    };

    BunchStateHandler() = default;

    // -- per-container registration ---------------------------------------

    /**
     * @brief Allocate a new per-container slot and return a strong reference
     *        to it. The handler itself retains only a `std::weak_ptr`, so the
     *        slot is destroyed automatically once the caller (typically a
     *        ParticleContainer) releases its shared_ptr.
     */
    std::shared_ptr<ContainerState> registerContainer();

    // -- bunch-wide flags --------------------------------------------------

    bool isFirstRepartition() const { return firstRepartition_m; }
    /**
     * @brief Non-const reference escape hatch used by LoadBalancer::repartition.
     *
     * Bypasses the aggressive-sync setter path. The load balancer is the single
     * well-defined caller that flips this flag from true to false on all ranks
     * after the first repartition.
     */
    bool& isFirstRepartitionRef() { return firstRepartition_m; }
    void setFirstRepartition(bool v);

    // -- emission liveness -------------------------------------------------
    // TBD: not wired up yet. Re-introduce a bunch-wide bool (plus impl in .cpp
    // and a unit test) once an emitting distribution actually needs it.
    // bool isEmittingNow() const;
    // void setEmittingNow(bool v);

private:
    // Weak refs to every slot handed out by registerContainer(). Used only by
    // handler-internal operations that iterate over all containers; pruned
    // lazily on iteration. Never exposed to callers.
    std::vector<std::weak_ptr<ContainerState>> registered_m;

    bool firstRepartition_m = true;
};

#endif
