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
 * Per-container flags (`momentsDirty`, `unitlessPositions`) live in the
 *   nested `ContainerState` struct. Each `ParticleContainer` registers itself
 *   via `registerContainer()` at `setBunchStateHandler` time and receives a
 *   `std::shared_ptr<ContainerState>` slot. The container holds the only strong
 *   reference; the handler keeps a `std::weak_ptr` so the slot is released
 *   automatically when the container is destroyed (no unregister needed).
 *
 * - **emittingNow** (aka "isEmitting"): managed by the emitting
 *   distribution itself (e.g. FlatTop sets it when t > t0 and
 *   particles are being created).
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
        void markMomentsClean();
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

private:
    // Weak refs to every slot handed out by registerContainer(). Used only by
    // handler-internal operations that iterate over all containers; pruned
    // lazily on iteration. Never exposed to callers.
    std::vector<std::weak_ptr<ContainerState>> registered_m;
};

#endif
