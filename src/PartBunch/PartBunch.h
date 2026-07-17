/**
 * @file PartBunch.h
 * @brief Template PIC bunch: IPPL PicManager, shared field mesh/solver, and multiple particle
 * containers.
 */

#ifndef PARTBUNCH_H
#define PARTBUNCH_H

#include <cstdint>
#include <memory>
#include <vector>

#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/Matrix.h"
#include "Algorithms/PartData.h"
#include "Attributes/Attributes.h"
#include "BCHandler.hpp"
#include "Manager/BaseManager.h"
#include "Manager/PicManager.h"
#include "PartBunch/BunchStateHandler.h"
#include "PartBunch/FieldContainer.hpp"
#include "PartBunch/FieldSolver.hpp"
#include "PartBunch/LoadBalancer.hpp"
#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"
#include "Random/Distribution.h"
#include "Random/InverseTransformSampling.h"
#include "Random/NormalDistribution.h"
#include "Random/Randn.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/OpalException.h"

class DataSink;
class Beam;

namespace opalx::spacecharge {
    class BinConfigurationObserver;
    template <typename T, unsigned Dim>
    class IterationPlan;
    class SelfFieldDiagnostics;
}  // namespace opalx::spacecharge

extern Inform* gmsg;

using view_type = typename ippl::detail::ViewType<ippl::Vector<double, 3>, 1>::view_type;

/**
 * @brief OPAL particle bunch: field container, solver, load balancer, and one or more beams.
 *
 * @tparam T   Floating-point type for positions/fields (typically double).
 * @tparam Dim Spatial dimension (3 for OPALX).
 */
template <typename T, unsigned Dim>
class BinnedFieldSolver;

template <typename T, unsigned Dim>
class PartBunch
    : public ippl::PicManager<
              T, Dim, ParticleContainer<T, Dim>, FieldContainer<T, Dim>, LoadBalancer<T, Dim>> {
public:
    using ParticleContainer_t = ParticleContainer<T, Dim>;
    using FieldContainer_t    = FieldContainer<T, Dim>;
    using PicWorkspace_t      = opalx::spacecharge::PicWorkspace<T, Dim>;
    using BinnedFieldSolver_t = BinnedFieldSolver<T, Dim>;
    using LoadBalancer_t      = LoadBalancer<T, Dim>;
    using Base                = ippl::ParticleBase<
                           ippl::ParticleSpatialLayout<T, Dim, ippl::UniformCartesian<T, Dim>>,
                           Kokkos::DefaultExecutionSpace::memory_space>;

    using BCHandler_t = BCHandler<Dim>;

public:
    // --- Shared state (all containers / mesh) ---

    double dt_m;                       ///< Global time step @f$\Delta t@f$ (s).
    int it_m;                          ///< Iteration counter (legacy / diagnostics).
    std::string integration_method_m;  ///< Integrator name (e.g. leapfrog).
    std::string solver_m;              ///< Field solver type string from input.
    Vector_t<int, Dim> nr_m;           ///< Mesh cell count per dimension.
    Vector_t<double, Dim> origin_m;    ///< Mesh origin (lab coordinates).
    Vector_t<double, Dim> rmin_m;  ///< Current bunch spatial minimum (from primary container stats;
                                   ///< see calcBeamParameters).
    Vector_t<double, Dim>
            rmax_m;              ///< Current bunch spatial maximum (from primary container stats).
    Vector_t<double, Dim> hr_m;  ///< Mesh spacing (m).

    double lbt_m;                    ///< Load-balancer timescale parameter.
    ippl::NDIndex<Dim> domain_m;     ///< Global mesh index extent per dimension.
    std::array<bool, Dim> decomp_m;  ///< Domain decomposition flags (per axis).

private:
    std::vector<bool> pcActive_m;   ///< Per-container: participate in this track segment.
    std::vector<bool> pcAtSStop_m;  ///< Per-container: frozen at current s-stop until next segment.
    std::vector<std::string> particleNames_m;  ///< Per-container beam particle names.

    std::vector<double> qi_m;  ///< Charge per macroparticle [C], one entry per container.
    std::vector<double> mi_m;  ///< Mass per macroparticle [GeV], one entry per container.

    const PartData* reference_m = nullptr;  ///< Reference particle data (set by TrackRun::execute).

    std::shared_ptr<BunchStateHandler>
            bunchState_m;  ///< Shared per-container coordinate and moment state.

    std::shared_ptr<BCHandler_t> bcHandler_m;  ///< Field boundary conditions.
    FieldSolverCmd* OPALFieldSolver_m;         ///< Borrowed parsed FIELD_SOLVER command.
    DataSink* dataSink_m;                      ///< Borrowed diagnostics and dump output sink.

    double t_m;  ///< Current simulation time (s).

    long long globalTrackStep_m;  ///< Global integration step counter.

    std::unique_ptr<size_t[]>
            globalPartPerNode_m;  ///< Per-rank particle counts for load-balance stats.

    double rmsDensity_m;  ///< Legacy RMS density placeholder (may still appear in stats).

public:
    /**
     * @brief Construct a multi-beam bunch: mesh, solver, containers, and capacity.
     *
     * @param qi                     Macrocharge per container (C).
     * @param mi                     Macromass per container (GeV/c²).
     * @param beams                  One @c Beam per container (reference particle, species).
     * @param totalParticlesPerBeam  Target macroparticle count per beam (for local allocation).
     * @param lbt                    Load-balancer timescale.
     * @param integration_method     Integrator label (e.g. leapfrog).
     * @param OPALFieldSolver        Borrowed field solver command (mesh, BCs, optional binning).
     * @param dataSink               Borrowed non-null diagnostics output sink.
     */
    PartBunch(
            std::vector<double> qi, std::vector<double> mi, const std::vector<Beam*>& beams,
            std::vector<size_t> totalParticlesPerBeam, double lbt, std::string integration_method,
            FieldSolverCmd* OPALFieldSolver, DataSink* dataSink);

    /**
     * @brief Recompute moments for every particle container without changing the PIC domain.
     *
     * Mesh geometry, field layout, particle migration, and backend refresh belong exclusively to
     * the concrete self-field algorithm. Tracking calls this after particle mutations when only
     * current statistics are required.
     */
    void updateAllParticleMoments();

    /**
     * @brief Set the image-charge configuration for the field solver.
     *
     * @param enabled Enable image-charge scatter mirror when true.
     * @param zPlane Mirror plane position in z [m].
     */
    void setImageChargeConfiguration(bool enabled, double zPlane);

    /**
     * @brief Set the shifted Green's function Dirichlet-correction configuration.
     *
     * Alternative to @c setImageChargeConfiguration. Mutually exclusive with it.
     * Requires the OPEN field solver (checked at runtime in the correction pass).
     *
     * @param enabled Enable the shifted-Green's-function correction when true.
     * @param zPlane  Dirichlet plane position in z [m].
     */
    void setShiftedGreensConfiguration(bool enabled, double zPlane);

    /**
     * @brief Configure diagnostic dump frequency for the ZEROFACE plane potential.
     *
     * @param frequency Dump every n-th global timestep. `0` disables dumping.
     */
    void setZeroFacePlaneDumpFrequency(int frequency);

    /// @brief Set the maximum number of timesteps for which image charges are active (0 =
    /// unlimited).
    void setZerofaceMaxSteps(int maxSteps);

    /**
     * @brief Sum of @c getTotalNum() over all particle containers.
     */
    size_t getTotalNumAllContainers() const {
        size_t total = 0;
        for (const auto& pc : this->getParticleContainers()) {
            if (pc) {
                total += pc->getTotalNum();
            }
        }
        return total;
    }

    /// @brief Build field solver and load balancer from @c OPALFieldSolver_m.
    void setSolver();

    /// @brief Warm-up: zero rho and run the field solver once (skip full dumps).
    void pre_run() override;

    /**
     * @brief Validate BC handler, solver wiring, field pointers, and layout extents.
     * @throw OpalException if initialization is inconsistent.
     */
    void performBunchSanityChecks() const;

    /// @brief At segment start: active if container is non-empty; inactive if empty.
    void resetPcActive();

    /// @param i Container index.
    /// @return Whether container @p i participates in the current segment.
    bool isPcActive(size_t i) const { return i < pcActive_m.size() && pcActive_m[i]; }

    /// @brief Force container @p i active (e.g. for containers with pending emission).
    void setPcActive(size_t i) {
        if (i < pcActive_m.size()) {
            pcActive_m[i] = true;
        }
    }

    /// @param i Container index.
    /// @return Whether container @p i is frozen at the current s-stop.
    bool pcAtSStop(size_t i) const { return i < pcAtSStop_m.size() && pcAtSStop_m[i]; }

    /// @brief Deactivate container @p i until the next step-size segment (s-stop reached).
    void setPcAtSStop(size_t i);

    /// @brief After emission: reactivate non-empty containers not marked at s-stop.
    void refreshPcActiveAfterEmit();

    /// @return True if any container is active.
    bool anyPcActive() const {
        for (bool a : pcActive_m) {
            if (a) {
                return true;
            }
        }
        return false;
    }

    /// @brief PicManager hook; throws (tracking does not use this path).
    void advance() override {
        throw OpalException(
                "PartBunch::advance",
                "Not used: just exists because ippl::PicManager wants it that way.");
    }

    /// @brief PicManager hook; throws (scatter handled elsewhere).
    void par2grid() override {
        throw OpalException(
                "PartBunch::par2grid",
                "Not used: just exists because ippl::PicManager wants it that way.");
    }

    /// @brief PicManager hook; throws (gather handled elsewhere).
    void grid2par() override {
        throw OpalException(
                "PartBunch::grid2par",
                "Not used: just exists because ippl::PicManager wants it that way.");
    }

    /// @brief Scratch E field used by the binned solver path.
    std::shared_ptr<VField_t<T, Dim>> getTempEField() {
        return this->fcontainer_m ? this->fcontainer_m->getTempEField() : nullptr;
    }

    /// @brief Scratch B field used by the binned solver path.
    std::shared_ptr<VField_t<T, Dim>> getTempBField() {
        return this->fcontainer_m ? this->fcontainer_m->getTempBField() : nullptr;
    }

    /// @param bcHandler Boundary-condition handler for the mesh.
    void setBCHandler(std::shared_ptr<BCHandler_t> bcHandler) { bcHandler_m = bcHandler; }

    /// @brief Current boundary-condition handler.
    std::shared_ptr<BCHandler_t> getBCHandler() const { return bcHandler_m; }

    DataSink* getDataSink() const { return dataSink_m; }

    std::shared_ptr<BunchStateHandler> getBunchStateHandler() { return bunchState_m; }
    std::shared_ptr<const BunchStateHandler> getBunchStateHandler() const { return bunchState_m; }

    /** @brief Share the stable PIC workspace with the configured self-field algorithm. */
    [[nodiscard]] std::shared_ptr<PicWorkspace_t> sharedPicWorkspace() {
        if (!this->fcontainer_m) {
            throw OpalException(
                    "PartBunch::sharedPicWorkspace", "FieldContainer is not initialized.");
        }
        return this->fcontainer_m->sharedWorkspace();
    }

    void updateMoments() { this->pcontainer_m->updateMoments(); }

    /**
     * @brief Update moments and set @c rmin_m / @c rmax_m from the primary (first) container.
     * @note Space-charge and mesh updates currently use container 0 only.
     */
    void calcBeamParameters();

    /**
     * @brief Set the per-particle charge for each particle container.
     * @note Copies values from `qi_m` into each particle container via `setQ`.
     * @note Throws if the number of particle containers and `qi_m` entries do not match.
     */
    void setCharge() {
        const auto& containers = this->getParticleContainers();
        if (containers.size() != qi_m.size()) {
            throw OpalException(
                    "PartBunch::setCharge",
                    "Number of particle containers and qi values do not match.");
        }
        for (size_t i = 0; i < containers.size(); ++i) {
            containers[i]->setQ(qi_m[i]);
        }
    }

    /**
     * @brief Set the per-particle mass for each particle container.
     * @note Copies values from `mi_m` into each particle container via `setM`.
     * @note Throws if the number of particle containers and `mi_m` entries do not match.
     */
    void setMass() {
        const auto& containers = this->getParticleContainers();
        if (containers.size() != mi_m.size()) {
            throw OpalException(
                    "PartBunch::setMass",
                    "Number of particle containers and mi values do not match.");
        }
        for (size_t i = 0; i < containers.size(); ++i) {
            containers[i]->setM(mi_m[i]);
        }
    }

    /**
     * @brief Get the total charge for a given particle container.
     * @param containerIndex Index of the particle container.
     * @returns `qi_m[containerIndex] * getParticleContainers()[containerIndex]->getTotalNum()`.
     * @note Throws if the number of particle containers and `qi_m` entries do not match, or if
     *       `containerIndex` is out of range.
     */
    double getCharge(size_t containerIndex = 0) const {
        const auto& containers = this->getParticleContainers();
        if (containers.size() != qi_m.size()) {
            throw OpalException(
                    "PartBunch::getCharge",
                    "Number of particle containers and qi values do not match.");
        }
        if (containerIndex >= containers.size()) {
            throw OpalException("PartBunch::getCharge", "Container index out of range.");
        }
        return qi_m[containerIndex] * containers[containerIndex]->getTotalNum();
    }

    /**
     * @brief Get the charge per particle for a given particle container.
     * @param containerIndex Index of the particle container.
     * @returns `qi_m[containerIndex]`.
     * @note Throws if `containerIndex` is out of range.
     */
    double getChargePerParticle(size_t containerIndex = 0) const {
        if (containerIndex >= qi_m.size()) {
            throw OpalException("PartBunch::getChargePerParticle", "Container index out of range.");
        }
        return qi_m[containerIndex];
    }

    /**
     * @brief Get the mass per particle for a given particle container.
     * @param containerIndex Index of the particle container.
     * @returns `mi_m[containerIndex]`.
     * @note Throws if `containerIndex` is out of range.
     */
    double getMassPerParticle(size_t containerIndex = 0) const {
        if (containerIndex >= mi_m.size()) {
            throw OpalException("PartBunch::getMassPerParticle", "Container index out of range.");
        }
        return mi_m[containerIndex];
    }

    /**
     * @brief Alias for `getCharge(containerIndex)`.
     * @param containerIndex Index of the particle container.
     * @returns Equivalent to `getCharge(containerIndex)`.
     */
    double getQ(size_t containerIndex = 0) const { return this->getCharge(containerIndex); }

    /**
     * @brief Get the total mass for a given particle container.
     * @param containerIndex Index of the particle container.
     * @returns `mi_m[containerIndex] * getParticleContainers()[containerIndex]->getTotalNum()`.
     * @note Throws if the number of particle containers and `mi_m` entries do not match, or if
     *       `containerIndex` is out of range.
     */
    double getM(size_t containerIndex = 0) const {
        const auto& containers = this->getParticleContainers();
        if (containers.size() != mi_m.size()) {
            throw OpalException(
                    "PartBunch::getM", "Number of particle containers and mi values do not match.");
        }
        if (containerIndex >= containers.size()) {
            throw OpalException("PartBunch::getM", "Container index out of range.");
        }
        return mi_m[containerIndex] * containers[containerIndex]->getTotalNum();
    }

    /**
     * @brief Get the total charge across all particle containers.
     * @returns `sum_i(qi_m[i] * containers[i]->getTotalNum())`.
     * @note Throws if the number of particle containers and `qi_m` entries do not match.
     */
    double getTotalCharge() const {
        const auto& containers = this->getParticleContainers();
        if (containers.size() != qi_m.size()) {
            throw OpalException(
                    "PartBunch::getTotalCharge",
                    "Number of particle containers and qi values do not match.");
        }
        double charge = 0.0;
        for (size_t i = 0; i < containers.size(); ++i) {
            charge += qi_m[i] * containers[i]->getTotalNum();
        }
        return charge;
    }

    /**
     * @brief Get the total mass across all particle containers.
     * @returns `sum_i(mi_m[i] * containers[i]->getTotalNum())`.
     * @note Throws if the number of particle containers and `mi_m` entries do not match.
     */
    double getTotalMass() const {
        const auto& containers = this->getParticleContainers();
        if (containers.size() != mi_m.size()) {
            throw OpalException(
                    "PartBunch::getTotalMass",
                    "Number of particle containers and mi values do not match.");
        }
        double mass = 0.0;
        for (size_t i = 0; i < containers.size(); ++i) {
            mass += mi_m[i] * containers[i]->getTotalNum();
        }
        return mass;
    }

    double getdE() {
        return this->pcontainer_m->getStdKineticEnergy();  // Unit: MeV
    }

    const PartData* getReference() const { return reference_m; }

    /// Set inside TrackRun::execute
    void setReference(const PartData* ref) {
        reference_m = ref;
        if (reference_m && this->pcontainer_m) {
            // Ensure mean/std kinetic energy in DistributionMoments are computed using reference
            // mass. PartData mass is stored in eV; DistributionMoments expects GeV for its energy
            // computation.
            this->pcontainer_m->setEnergyReferenceMass(reference_m->getM() * Units::eV2GeV, true);
        }
    }

    void gatherLoadBalanceStatistics();

    size_t getLoadBalance(int p) { return globalPartPerNode_m[p]; }

    /// @brief Stub; logs and returns 0.
    size_t calcNumPartsOutside(Vector_t<double, Dim> /*x*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__
              << " function: " << __func__ << endl;
        return 0;
    }

    /// @brief Stub; logs only.
    void calcLineDensity(
            unsigned int /*nBins*/, std::vector<double>& /*lineDensity*/,
            std::pair<double, double>& /*meshInfo*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__
              << " function: " << __func__ << endl;
    }

    /// @brief Stub; logs and returns zero vector.
    Vector_t<double, Dim> getEExtrema() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__
              << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0);
    }

    /**
     * @brief Compute the bunch self-fields (binned when available).
     *
     * The actual implementation lives in the solver object (see `BinnedFieldSolver`).
     * `ParallelTracker` only orchestrates reference/beam-frame transforms and calls
     * this delegator once per step.
     *
     * @param iterationPlan Solver-owned bin and pass plan for the current field evaluation.
     * @param particleGeneration Generation of the current particle population.
     * @param binConfigurationObserver Optional sink for bin configuration snapshots.
     * @param diagnostics Self-field timing and event diagnostics.
     */
    void computeSelfFields(
            opalx::spacecharge::IterationPlan<T, Dim>& iterationPlan,
            std::uint64_t particleGeneration,
            opalx::spacecharge::BinConfigurationObserver* binConfigurationObserver,
            opalx::spacecharge::SelfFieldDiagnostics& diagnostics);

    /// @brief Human-readable dump of each container to @p os.
    /// @param os Output stream wrapper.
    /// @return Reference to @p os.
    Inform& print(Inform& os);

    /// @return True if a field solver instance is installed.
    bool hasFieldSolver() const { return this->fsolver_m != nullptr; }

    /**
     * @brief Non-const pointer to the concrete BinnedFieldSolver.
     * @note Definition in PartBunch.cpp; requires complete BinnedFieldSolver type.
     */
    BinnedFieldSolver_t* getFieldSolver();

    /// @brief Const overload of getFieldSolver().
    const BinnedFieldSolver_t* getFieldSolver() const;

    /// @brief Backend type string (e.g. FFT, OPEN, CG, NONE).
    std::string getFieldSolverType();

    /**
     * @brief Effective bin count reported by the field solver for legacy diagnostics.
     */
    int getCurrentNBins() const;

    /// @brief Compatibility stub; logs and returns 0.
    double calcMeanPhi() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__
              << " function: " << __func__ << endl;
        return 0.0;
    }

    /// @brief Particle species name for container @p i (from BEAM PARTICLE input).
    std::string getParticleName(size_t i) const {
        if (i >= particleNames_m.size()) {
            return "";
        }
        return particleNames_m[i];
    }

    /// @brief Do not use; throws (access positions via @c ParticleContainer::R).
    Vector_t<double, Dim> R(size_t) {
        throw OpalException(
                "PartBunch::R",
                "Not implemented: shouldn't be called, since this is not the correct way to access "
                "particle positions.");
        return Vector_t<double, Dim>(0.0);
    }

    /// @brief Do not use; throws (access momenta via @c ParticleContainer::P).
    Vector_t<double, Dim> P(size_t) {
        throw OpalException(
                "PartBunch::P",
                "Not implemented: shouldn't be called, since this is not the correct way to access "
                "particle momenta.");
        return Vector_t<double, Dim>(0.0);
    }

    /**
     * @brief Copy cached bunch extent (@c rmin_m, @c rmax_m) from calcBeamParameters.
     * @param[out] rmin Lower corner.
     * @param[out] rmax Upper corner.
     */
    void get_bounds(Vector_t<double, Dim>& rmin, Vector_t<double, Dim>& rmax) {
        rmin = rmin_m;
        rmax = rmax_m;
    }

    /// @brief Set the global time step.
    void setdT(double dt) { dt_m = dt; }

    /// @brief Get the global time step.
    double getdT() const { return dt_m; }

    /// @brief Set the current simulation time.
    void setT(double t) { t_m = t; }

    /// @brief Advance time by one global time step.
    void incrementT() { t_m += dt_m; }

    /// @brief Get the current simulation time.
    double getT() const { return t_m; }

    /// @brief Cached minimum extent (@c rmin_m); prefer per-container min/max for multi-beam
    /// detail.
    Vector_t<double, Dim> get_origin() const { return rmin_m; }
    /// @brief Cached maximum extent (@c rmax_m).
    Vector_t<double, Dim> get_maxExtent() const { return rmax_m; }

    /// @brief Stub; logs and returns zero vector.
    Vector_t<double, Dim> get_halo() const {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__
              << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }

    /// @brief Get mesh spacing
    Vector_t<double, Dim> get_hr() const { return hr_m; }

    /// @param n Global step counter to store.
    void setGlobalTrackStep(long long n) { globalTrackStep_m = n; }

    /// @brief Current global tracking step.
    long long getGlobalTrackStep() const { return globalTrackStep_m; }

    /// @brief Increment @c globalTrackStep_m by one.
    void incTrackSteps() { globalTrackStep_m++; }

    /// @brief Legacy RMS density field (may be unused).
    double get_rmsDensity() const { return rmsDensity_m; }
};

extern template class PartBunch<double, 3>;  ///< Explicit instantiation for 3D double-precision.

#endif  // PARTBUNCH_H
