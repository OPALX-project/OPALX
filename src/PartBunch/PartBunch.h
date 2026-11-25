/**
 * @file PartBunch.h
 * @brief All beams of one RUN: one particle container per BEAM and what they share.
 */

#ifndef PARTBUNCH_H
#define PARTBUNCH_H

#include <cstddef>
#include <functional>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "Moments/BeamMoments.h"
#include "PartBunch/BunchStateHandler.h"
#include "PartBunch/CartesianDomain.h"
#include "PartBunch/CartesianDomainConfig.h"
#include "PartBunch/ParticleContainer.hpp"

namespace opalx::run {
    struct BeamDescription;
}  // namespace opalx::run

/**
 * @brief All beams of one RUN.
 *
 * PartBunch owns one ParticleContainer per BEAM and what the beams share: the mesh and field
 * layout (CartesianDomain), one time, one time step and one step counter. Per container it keeps
 * whether the container is active, whether its beam has passed its ZSTOP, which container's
 * reference particle it follows (a decay daughter follows the first container of its decay
 * chain), and the macro charge and mass.
 *
 * The particles, the reference particle, the path length and the bunch-to-lab transform are not
 * held here; they live in each ParticleContainer.
 *
 * The step functions (pushParticles(), kickAndPushParticles(), ...) loop over the containers and
 * skip those that are not active.
 *
 * @tparam T   Floating-point type for positions/fields (typically double).
 * @tparam Dim Spatial dimension (3 for OPALX).
 */
template <typename T, unsigned Dim>
class PartBunch {
public:
    using ParticleContainer_t = ParticleContainer<T, Dim>;
    using CartesianDomain_t   = opalx::spacecharge::CartesianDomain<T, Dim>;
    using binIndex_t          = typename ParticleContainer_t::bin_index_type;

    /**
     * @brief Construct a multi-beam bunch and its initial Cartesian particle layout.
     *
     * @param qi                     Macrocharge per container (C).
     * @param mi                     Macromass per container (GeV/c²).
     * @param beams                  One beam per container: reference particle, PARTICLE name
     *                               and whether it tracks spin (POLARIZATION).
     * @param totalParticlesPerBeam  Target macroparticle count per beam (for local allocation).
     * @param lbt                    OPTION LOADBALANCINGTHRESHOLD; only printed by print().
     * @param domainConfig           Immutable Cartesian domain and particle-layout setup.
     * @param qmAttributes           Store Q and M per particle (OPTION QM_MODE = ATTRIBUTES).
     * @param aggressiveStateSync    OPTION AGGRESSIVE_STATE_SYNC, see BunchStateHandler.
     */
    PartBunch(
            std::vector<double> qi, std::vector<double> mi,
            const std::vector<opalx::run::BeamDescription>& beams,
            std::vector<size_t> totalParticlesPerBeam, double lbt,
            opalx::spacecharge::CartesianDomainConfig<T, Dim> domainConfig,
            bool qmAttributes = false, bool aggressiveStateSync = false);

    /**
     * @brief Sum of @c getTotalNum() over all particle containers.
     */
    size_t getTotalNumAllContainers() const {
        size_t total = 0;
        for (const auto& pc : pcontainers_m) {
            total += pc->getTotalNum();
        }
        return total;
    }

    /// @brief Make every non-empty container active and every empty one inactive, and clear the
    ///        ZSTOP flag of every container (see pcAtSStop()). Called at the start of each DT
    ///        segment; BeamStops::keepDoneBeamsStopped() then sets the flag again for beams that
    ///        are done.
    void resetPcActive();

    /// @param i Container index.
    /// @return Whether container @p i participates in the current segment.
    /// @throws std::out_of_range if @p i is not a container index.
    bool isPcActive(size_t i) const { return pcActive_m.at(i); }

    /// @brief Force container @p i active (e.g. for containers with pending emission).
    /// @throws std::out_of_range if @p i is not a container index.
    void setPcActive(size_t i) { pcActive_m.at(i) = true; }

    /// @param i Container index.
    /// @return Whether the beam of container @p i has passed its ZSTOP.
    /// @throws std::out_of_range if @p i is not a container index.
    bool pcAtSStop(size_t i) const { return pcAtSStop_m.at(i); }

    /// @brief Mark the beam of container @p i as past its ZSTOP and make the container inactive.
    ///        Called by BeamStops and by a checkpoint restore.
    /// @throws std::out_of_range if @p i is not a container index.
    void setPcAtSStop(size_t i);

    /// @brief Make every non-empty container active whose beam has not passed its ZSTOP. Called
    ///        after emission and after a decay.
    void refreshPcActiveAfterEmit();

    /// @brief Make container @p i use the reference particle of container @p root.
    /// @throws std::out_of_range if @p i or @p root is not a container index.
    void setPcRoot(size_t i, size_t root) {
        if (root >= pcRoot_m.size()) {
            throw std::out_of_range("PartBunch::setPcRoot: root index out of range");
        }
        pcRoot_m.at(i) = root;
    }

    /// @param i Container index.
    /// @return Index of the container whose reference particle container @p i uses.
    /// @throws std::out_of_range if @p i is not a container index.
    size_t getPcRoot(size_t i) const { return pcRoot_m.at(i); }

    /// @param i Container index.
    /// @return Whether container @p i is a decay daughter container that uses another
    ///   container's reference particle, frame, path length and orbit.
    bool isFollower(size_t i) const { return getPcRoot(i) != i; }

    /// @param i Container index.
    /// @return Whether the reference particle of container @p i must be advanced: it is
    ///   active, or a container that uses its reference particle is active.
    bool referenceMoves(size_t i) const {
        if (isPcActive(i)) {
            return true;
        }
        for (size_t f = 0; f < pcRoot_m.size(); ++f) {
            if (f != i && pcRoot_m[f] == i && isPcActive(f)) {
                return true;
            }
        }
        return false;
    }

    /// @throws OpalException from @p where if a container's particles are not in its reference
    ///         frame.
    void requireReferenceFrame(const std::string& where) const {
        for (const auto& pc : pcontainers_m) {
            pc->requireReferenceFrame(where);
        }
    }

    // --- Particle push, fields and loss over the containers ---

    /// @brief First half of the leapfrog step: half drift of every active container.
    void pushParticles();

    /// @brief Second half of the leapfrog step for every active container: kick, second half of
    ///        the spin rotation with the momentum after the kick, half drift; then each
    ///        container's particle dt is set to the bunch time step.
    void kickAndPushParticles();

    /// @brief First half of the Thomas-BMT spin rotation for every active container that stores
    ///        Pol, with the momentum before the kick. kickAndPushParticles() does the second
    ///        half, so the spin update is second order in dt. Call after the external and
    ///        space-charge fields at the half-step position have been accumulated.
    void evolveSpinHalfStep();

    /// @brief Zero E and B on every active container.
    void resetFields();

    /// @brief Apply the global processes of every container (e.g. decay) over the bunch time
    ///        step getdT() and switch on daughter containers that received particles.
    ///        Collective.
    /// @return Global number of particles the processes marked invalid.
    size_t applyGlobalProcesses();

    /**
     * @brief Mark the particles of each active, non-empty container outside the box of
     *        @p sigmas rms sizes around its mean position (OPTION BOUNDPDESTROY). Collective.
     *
     * Deletion waits for deleteInvalidParticles().
     * @param request Moment groups to compute per container in addition to those of the box,
     *        so the caller can reuse the moments.
     * @param[in,out] nMarked Gets the global number of newly marked particles added.
     * @return Per container the moments computed for the box; empty for a container without box
     *         and for all if @p sigmas is not positive.
     */
    std::vector<opalx::moments::BeamMoments> markOutliers(
            double sigmas,
            const std::function<opalx::moments::Request(const ParticleContainer_t&)>& request,
            size_t& nMarked);

    /// @brief Delete the particles marked invalid in every container, or in the active ones
    ///        only, and log the count per container with @p reason. Collective.
    void deleteInvalidParticles(bool activeOnly, Inform& m, const std::string& reason);

    /// @return True if any container is active.
    bool anyPcActive() const {
        for (bool a : pcActive_m) {
            if (a) {
                return true;
            }
        }
        return false;
    }

    std::shared_ptr<BunchStateHandler> getBunchStateHandler() { return bunchState_m; }
    std::shared_ptr<const BunchStateHandler> getBunchStateHandler() const { return bunchState_m; }

    /** @brief Domain whose mesh and layout back every particle container in this bunch. */
    [[nodiscard]] CartesianDomain_t& cartesianDomain() { return *cartesianDomain_m; }
    [[nodiscard]] const CartesianDomain_t& cartesianDomain() const { return *cartesianDomain_m; }

    /// @return The first particle container.
    [[nodiscard]] std::shared_ptr<ParticleContainer_t> getParticleContainer() const {
        return pcontainers_m.front();
    }

    /// @return Particle container @p index.
    /// @throws std::out_of_range if @p index is not a container index.
    [[nodiscard]] std::shared_ptr<ParticleContainer_t> getParticleContainer(size_t index) const {
        return pcontainers_m.at(index);
    }

    [[nodiscard]] const std::vector<std::shared_ptr<ParticleContainer_t>>& getParticleContainers()
            const {
        return pcontainers_m;
    }

    [[nodiscard]] size_t getNumParticleContainers() const { return pcontainers_m.size(); }

    /**
     * @brief Set the per-particle charge for each particle container.
     * @note Copies values from `qi_m` into each particle container via `setQ`.
     */
    void setCharge() {
        for (size_t i = 0; i < pcontainers_m.size(); ++i) {
            pcontainers_m[i]->setQ(qi_m[i]);
        }
    }

    /**
     * @brief Set the per-particle mass for each particle container.
     * @note Copies values from `mi_m` into each particle container via `setM`.
     */
    void setMass() {
        for (size_t i = 0; i < pcontainers_m.size(); ++i) {
            pcontainers_m[i]->setM(mi_m[i]);
        }
    }

    /**
     * @brief Human-readable dump of each container to @p os.
     * @param moments One snapshot per container holding at least printedMoments.
     */
    Inform& print(Inform& os, std::span<const opalx::moments::BeamMoments> moments) const;

    /// Moment groups print() reads.
    static constexpr opalx::moments::Request printedMoments =
            opalx::moments::Moment::Bounds | opalx::moments::Moment::MeanPosition
            | opalx::moments::Moment::MeanMomentum | opalx::moments::Moment::Spread;

    /// @brief Particle species name for container @p i (from BEAM PARTICLE input).
    /// @throws std::out_of_range if @p i is not a container index.
    std::string getParticleName(size_t i) const { return particleNames_m.at(i); }

    /// @brief Set the time step of the bunch only; the particle dt in the containers is not
    ///        changed (see setTimeStep()).
    void setdT(double dt) { dt_m = dt; }

    /// @brief Set the time step of the bunch and the particle dt of every particle in every
    ///        container to @p dt (s).
    void setTimeStep(double dt);

    /// @brief Get the global time step.
    double getdT() const { return dt_m; }

    /// @brief Set the current simulation time.
    void setT(double t) { t_m = t; }

    /// @brief Advance time by one global time step.
    void incrementT() { t_m += dt_m; }

    /// @brief Get the current simulation time.
    double getT() const { return t_m; }

    /// @param n Global step counter to store.
    void setGlobalTrackStep(long long n) { globalTrackStep_m = n; }

    /// @brief Current global tracking step.
    long long getGlobalTrackStep() const { return globalTrackStep_m; }

    /// @brief Increment @c globalTrackStep_m by one.
    void incTrackSteps() { globalTrackStep_m++; }

private:
    double dt_m;   ///< Global time step @f$\Delta t@f$ (s).
    double lbt_m;  ///< OPTION LOADBALANCINGTHRESHOLD; only printed by print().

    std::vector<bool> pcActive_m;   ///< Per-container: participate in this track segment.
    std::vector<bool> pcAtSStop_m;  ///< Per-container: the beam has passed its ZSTOP.
    std::vector<std::string> particleNames_m;  ///< Per-container beam particle names.
    /// Per-container index of the container whose reference particle it uses: its own index,
    /// or for a decay daughter container the first container of its decay chain.
    std::vector<size_t> pcRoot_m;

    std::vector<double> qi_m;  ///< Charge per macroparticle [C], one entry per container.
    std::vector<double> mi_m;  ///< Mass per macroparticle [GeV], one entry per container.

    std::shared_ptr<BunchStateHandler>
            bunchState_m;  ///< Shared per-container coordinate state and fixed bounds.

    std::unique_ptr<CartesianDomain_t>
            cartesianDomain_m;  ///< Outlives every particle layout that borrows it.
    std::vector<std::shared_ptr<ParticleContainer_t>> pcontainers_m;  ///< One per BEAM.

    double boundingBoxIncreasePercent_m;  ///< Domain bounding-box increase (%); only printed.

    double t_m;  ///< Current simulation time (s).

    long long globalTrackStep_m;  ///< Global integration step counter.
};

extern template class PartBunch<double, 3>;  ///< Explicit instantiation for 3D double-precision.

#endif  // PARTBUNCH_H
