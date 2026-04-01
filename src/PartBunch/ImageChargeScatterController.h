#ifndef OPALX_IMAGE_CHARGE_SCATTER_CONTROLLER_H
#define OPALX_IMAGE_CHARGE_SCATTER_CONTROLLER_H

#include "Ippl.h"
#include "PartBunch/Binning/AdaptBins.h"
#include "PartBunch/ParticleContainer.hpp"

/**
 * @brief Orchestrates primary and image-charge scatter deposition.
 *
 * The controller performs:
 * 1) primary scatter using the existing `dt * Q` deposition workflow, and
 * 2) optional image scatter by reflecting particle z positions around an xy plane
 *    and flipping the sign of charge `Q`.
 *
 * The particle state is restored after the image deposition pass.
 *
 * @tparam T   Scalar type (typically `double`).
 * @tparam Dim Spatial dimension (currently constrained to 3).
 */
template <typename T, unsigned Dim>
class ImageChargeScatterController {
    static_assert(Dim == 3, "ImageChargeScatterController currently supports Dim == 3 only.");

public:
    using ParticleCtr_t  = ParticleContainer<T, Dim>;
    using PositionAttr_t = typename ippl::ParticleBase<
            ippl::ParticleSpatialLayout<T, Dim>>::particle_position_type;
    using RhoField_t  = Field_t<Dim>;
    using BinPolicy_t = Kokkos::RangePolicy<>;
    using Hash_t      = typename ParticleBinning::AdaptBinsBase<ParticleCtr_t>::hash_type;

public:
    /**
     * @brief Construct disabled controller at `zPlane=0`.
     */
    ImageChargeScatterController() = default;

    /**
     * @brief Construct controller with explicit image-charge settings.
     * @param enabled Enable image scatter pass when true.
     * @param zPlane  Mirror plane location in z (meters).
     */
    ImageChargeScatterController(bool enabled, double zPlane)
        : enabled_m(enabled), zPlane_m(zPlane) {}

    /**
     * @brief Update runtime settings.
     * @param enabled Enable image scatter pass when true.
     * @param zPlane  Mirror plane location in z (meters).
     */
    void configure(bool enabled, double zPlane) {
        enabled_m = enabled;
        zPlane_m  = zPlane;
    }

    /// @brief Returns whether image-charge scatter is enabled.
    bool isEnabled() const { return enabled_m; }
    /// @brief Returns the mirror plane position in z.
    double getZPlane() const { return zPlane_m; }

    /**
     * @brief Test helper: apply image transform to all local particles.
     *
     * Mirrors z around the configured plane and flips Q sign.
     */
    void applyMirrorForTesting(std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions) const {
        applyMirrorTransformAll(pc, positions);
    }

    /**
     * @brief Test helper: restore image transform for all local particles.
     *
     * This is the inverse of `applyMirrorForTesting`.
     */
    void restoreMirrorForTesting(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions) const {
        restoreMirrorTransformAll(pc, positions);
    }

    /**
     * @brief Scatter all local particles with optional image pass.
     *
     * Primary pass deposits with the nominal charge.
     * If enabled, image pass deposits mirrored particles with flipped charge sign.
     *
     * @param pc        Shared particle container.
     * @param positions Particle positions attribute (typically `pc->R`).
     * @param rho       Target charge-density mesh field.
     */
    void scatterPrimaryAndImage(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho) const;

    /**
     * @brief Scatter a bin-restricted particle subset with optional image pass.
     *
     * @param pc        Shared particle container.
     * @param positions Particle positions attribute (typically `pc->R`).
     * @param rho       Target charge-density mesh field.
     * @param policy    Bin iteration policy (index-space for hash access).
     * @param hash      Mapping from policy index to particle index.
     */
    void scatterPrimaryAndImage(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho,
            const BinPolicy_t& policy, const Hash_t& hash) const;

private:
    /// @brief Apply z-mirror + charge sign flip for all local particles.
    void applyMirrorTransformAll(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions) const;

    /// @brief Restore z-mirror + charge sign flip for all local particles.
    void restoreMirrorTransformAll(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions) const;

    /// @brief Apply z-mirror + charge sign flip for a hashed subset.
    void applyMirrorTransformSubset(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, const BinPolicy_t& policy,
            const Hash_t& hash) const;

    /// @brief Restore z-mirror + charge sign flip for a hashed subset.
    void restoreMirrorTransformSubset(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, const BinPolicy_t& policy,
            const Hash_t& hash) const;

    /// @brief Flip charge sign for all local particles (or shared scalar Q in single mode).
    void flipChargeSignAll(std::shared_ptr<ParticleCtr_t> pc) const;

    /// @brief Flip charge sign for subset particles (or shared scalar Q in single mode).
    void flipChargeSignSubset(
            std::shared_ptr<ParticleCtr_t> pc, const BinPolicy_t& policy, const Hash_t& hash) const;

private:
    bool enabled_m  = false;
    double zPlane_m = 0.0;
};

#include "PartBunch/ImageChargeScatterController.tpp"

#endif
