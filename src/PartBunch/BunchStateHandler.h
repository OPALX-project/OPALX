#ifndef BUNCH_STATE_HANDLER_H
#define BUNCH_STATE_HANDLER_H

/**
 * @class BunchStateHandler
 * @brief Centralised tracking of mutable bunch-level status flags.
 *
 * Owned by PartBunch as a std::shared_ptr and passed to every class
 * that needs access (ParticleContainer, SamplingBase subclasses, etc.).
 *
 * Invariants:
 *
 * - **unitlessPositions**: true while particle positions are in the
 *   dimensionless coordinate system R' = R / (c * dt).  Only
 *   switchToUnitlessPositions / switchOffUnitlessPositions toggle this.
 *   Does NOT mark the bunch dirty (coordinate representation change only).
 *
 * - **momentsDirty** (aka "isDirty"): set whenever a physics operation
 *   mutates particle positions (R) or momenta (P) -- push, kick,
 *   emission, particle deletion. Cleared by
 *   ParticleContainer::updateMoments() after
 *   DistributionMoments::computeMoments() completes.
 *
 * - **emittingNow** (aka "isEmitting"): managed by the emitting
 *   distribution itself (e.g. FlatTop sets it when t > t0 and
 *   particles are being created).
 */
class BunchStateHandler {
public:
    BunchStateHandler() = default;

    // -- coordinate representation -----------------------------------------
    bool isUnitlessPositions() const { return unitlessPositions_m; }
    void setUnitlessPositions(bool v) { unitlessPositions_m = v; }

    // -- moments cache invalidation ----------------------------------------
    bool isMomentsDirty() const { return momentsDirty_m; }
    void markMomentsDirty() { momentsDirty_m = true; }
    void clearMomentsDirty() { momentsDirty_m = false; }

    // -- emission liveness -------------------------------------------------
    bool isEmittingNow() const { return emittingNow_m; }
    void setEmittingNow(bool v) { emittingNow_m = v; }

    // -- first repartition guard -------------------------------------------
    bool isFirstRepartition() const { return firstRepartition_m; }
    bool& isFirstRepartitionRef() { return firstRepartition_m; }
    void setFirstRepartition(bool v) { firstRepartition_m = v; }

private:
    bool unitlessPositions_m = false;
    bool momentsDirty_m      = true;
    bool emittingNow_m       = false;
    bool firstRepartition_m  = true;
};

#endif
