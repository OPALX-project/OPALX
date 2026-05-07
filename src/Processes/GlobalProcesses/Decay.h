#ifndef OPAL_GLOBALPROCESS_DECAY_H
#define OPAL_GLOBALPROCESS_DECAY_H

#include "Processes/GlobalProcesses/GlobalProcess.h"

#include "Ippl.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include <cstddef>
#include <memory>

/**
 * @brief Abstract base class for particle decay processes.
 *
 * @note Handles the shared mechanics of relativistic decay:
 *   1. Mark particles as decayed (exponential law with time dilation)
 *   2. Collect parent kinematics into compact temporary views
 *   3. Delegate daughter creation to the subclass via createDaughterParticles()
 *   4. Leave parent deletion to ParticleContainer::deleteInvalidParticles()
 *
 *   Subclasses implement the physics-specific daughter momentum sampling
 *   (e.g. Michel spectrum for muon decay, fixed two-body kinematics for pion
 *   decay). If no daughter container is set, decayed particles are simply
 *   marked for deletion.
 */
class Decay : public GlobalProcess {
public:
    Decay(double restLifetimeSeconds, std::size_t containerIndex, double parentMassGeV);

    ~Decay() override = default;

    /**
     * @brief Set the daughter particle container and its rest mass.
     * @param daughterPC  Container that receives decay products.
     * @param daughterMassGeV  Rest mass of the daughter particle [GeV].
     *
     * @note When set, decayed particles spawn a daughter in @p daughterPC.
     *   When not set (default), decayed particles are marked for deletion.
     */
    void setDaughterContainer(
            std::shared_ptr<ParticleContainer<double, 3>> daughterPC, double daughterMassGeV);

    size_t apply(
            ParticleContainer<double, 3>& pc, double dt, long long globalTrackStep,
            size_t containerIdx) override;

    /**
     * @brief Create daughter particles from collected parent data.
     *
     * @note Called only when daughterPC_m is set and localDestroyNum > 0.
     *   Subclasses implement the physics-specific momentum sampling here.
     *   Must be public because device compilers require that a Kokkos
     *   lambda's enclosing member function has public access.
     */
    virtual void createDaughterParticles(
            std::size_t localDestroyNum, std::size_t oldDaughterLocal,
            const Kokkos::View<ippl::Vector<double, 3>*>& parentR,
            const Kokkos::View<ippl::Vector<double, 3>*>& parentP,
            const Kokkos::View<double*>& parentDt) = 0;

    /// Compact views of the kinematics of parents marked for decay.
    struct DecayedParentViews {
        Kokkos::View<ippl::Vector<double, 3>*> R;
        Kokkos::View<ippl::Vector<double, 3>*> P;
        Kokkos::View<double*> dt;
    };

    /**
     * @brief Mark particles for decay using the relativistic exponential law.
     *
     * @note `tau0_m` and `randPool_m` are copied into locals before the device
     *   kernel so the KOKKOS_LAMBDA captures values — never `this`.
     */
    ippl::detail::size_type markDecayedParticles(
            ippl::detail::size_type nLocal, double dt, Kokkos::View<ippl::Vector<double, 3>*> Pview,
            Kokkos::View<bool*> decayed);

    /**
     * @brief Gather R/P/dt of parents marked for decay into compact views.
     *
     * @note Uses no class state, so no capture shields are needed; still a
     *   member for organizational grouping with @ref markDecayedParticles.
     */
    DecayedParentViews collectDecayedParents(
            ippl::detail::size_type nLocal, ippl::detail::size_type localDestroyNum,
            Kokkos::View<bool*> decayed, Kokkos::View<ippl::Vector<double, 3>*> Rview,
            Kokkos::View<ippl::Vector<double, 3>*> Pview, Kokkos::View<double*> dtView);

protected:
    /// Mean lifetime at rest [s].
    double tau0_m;

    /// Random pool for decay sampling.
    Kokkos::Random_XorShift64_Pool<> randPool_m;

    /// Daughter container for decay products (nullptr = destroy-only mode).
    std::shared_ptr<ParticleContainer<double, 3>> daughterPC_m;

    /// Physically allowed daughter species (short value of ParticleType enum).
    /// Subclasses should set this in their constructor; default = ParticleType::UNNAMED (-1).
    short allowedDaughterSpecies_m = -1;

    /// Rest mass of the daughter particle [GeV].
    double daughterMassGeV_m = 0.0;

    /// Rest mass of the parent particle [GeV].
    double parentMassGeV_m;
};

#endif
