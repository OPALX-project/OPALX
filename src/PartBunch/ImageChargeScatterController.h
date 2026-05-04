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
            ippl::ParticleSpatialLayout<T, Dim>,
            Kokkos::DefaultExecutionSpace::memory_space>::particle_position_type;
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

    /**
     * @brief Scatter only the primary (real) charges for all local particles.
     */
    void scatterPrimaryOnly(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho) const;

    /**
     * @brief Scatter only the primary (real) charges for a bin-restricted subset.
     */
    void scatterPrimaryOnly(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho,
            const BinPolicy_t& policy, const Hash_t& hash) const;

    /**
     * @brief Scatter only the image (mirrored) charges for all local particles.
     *
     * Applies the mirror transform, scatters, and restores. Does nothing if image
     * charges are disabled.
     */
    void scatterImageOnly(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho) const;

    /**
     * @brief Scatter only the image (mirrored) charges for a bin-restricted subset.
     *
     * Applies the mirror transform, scatters, and restores. Does nothing if image
     * charges are disabled.
     */
    void scatterImageOnly(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho,
            const BinPolicy_t& policy, const Hash_t& hash) const;

    // ==== All of the following functions could be private, but need to be public for device
    // compilation ====

    /**
     * @brief Scatter all local particles using the standard `dt*Q` workflow.
     *
     * This helper applies `scaleDtByCharge()`, scatters `dt`, and restores `dt`.
     *
     * @param pc Shared particle container.
     * @param positions Particle positions attribute.
     * @param rho Target charge-density field.
     */
    void scatterScaledDtAll(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho) const;

    /**
     * @brief Scatter a hashed particle subset using the standard `dt*Q` workflow.
     *
     * @param pc Shared particle container.
     * @param positions Particle positions attribute.
     * @param rho Target charge-density field.
     * @param policy Bin iteration policy.
     * @param hash Mapping from policy index to particle index.
     */
    void scatterScaledDtSubset(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho,
            const BinPolicy_t& policy, const Hash_t& hash) const;

    /**
     * @brief Apply z-mirror and charge sign flip for all local particles.
     * @param pc Shared particle container.
     * @param positions Particle positions attribute.
     */
    void applyMirrorTransformAll(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions) const;

    /**
     * @brief Restore z-mirror and charge sign flip for all local particles.
     *
     * The operation is its own inverse.
     *
     * @param pc Shared particle container.
     * @param positions Particle positions attribute.
     */
    void restoreMirrorTransformAll(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions) const;

    /**
     * @brief Apply z-mirror and charge sign flip for a hashed particle subset.
     * @param pc Shared particle container.
     * @param positions Particle positions attribute.
     * @param policy Bin iteration policy.
     * @param hash Mapping from policy index to particle index.
     */
    void applyMirrorTransformSubset(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, const BinPolicy_t& policy,
            const Hash_t& hash) const;

    /**
     * @brief Restore z-mirror and charge sign flip for a hashed particle subset.
     *
     * The operation is its own inverse.
     *
     * @param pc Shared particle container.
     * @param positions Particle positions attribute.
     * @param policy Bin iteration policy.
     * @param hash Mapping from policy index to particle index.
     */
    void restoreMirrorTransformSubset(
            std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, const BinPolicy_t& policy,
            const Hash_t& hash) const;

    /**
     * @brief Flip charge sign for all local particles.
     *
     * In single-value Q mode, only the shared scalar is flipped.
     *
     * @param pc Shared particle container.
     */
    void flipChargeSignAll(std::shared_ptr<ParticleCtr_t> pc) const;

    /**
     * @brief Flip charge sign for a hashed particle subset.
     *
     * In single-value Q mode, only the shared scalar is flipped.
     *
     * @param pc Shared particle container.
     * @param policy Bin iteration policy.
     * @param hash Mapping from policy index to particle index.
     */
    void flipChargeSignSubset(
            std::shared_ptr<ParticleCtr_t> pc, const BinPolicy_t& policy, const Hash_t& hash) const;

private:
    bool enabled_m  = false;
    double zPlane_m = 0.0;
};

#include "PartBunch/ImageChargeScatterController.tpp"

#endif
