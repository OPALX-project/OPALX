/**
 * @file PartBunch.h
 * @brief Template PIC bunch: IPPL PicManager, shared field mesh/solver, and multiple particle containers.
 */

#ifndef PARTBUNCH_H
#define PARTBUNCH_H

#include <memory>
#include <vector>

#include "Algorithms/Matrix.h"
#include "Algorithms/CoordinateSystemTrafo.h"
#include "Attributes/Attributes.h"
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
#include "Utilities/OpalException.h"
#include "BCHandler.hpp"
#include "Structure/FieldSolverCmd.h"
#include "Algorithms/PartData.h"
#include "PartBunch/Binning/AdaptBins.h"

class DataSink;  ///< Forward declaration; full definition only required in the .cpp translation unit.
class Beam;

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
class PartBunch : public ippl::PicManager<
          T, Dim, ParticleContainer<T, Dim>, FieldContainer<T, Dim>, LoadBalancer<T, Dim>> {
public:
    using ParticleContainer_t = ParticleContainer<T, Dim>;
    using FieldContainer_t    = FieldContainer<T, Dim>;
    using BinnedFieldSolver_t = BinnedFieldSolver<T, Dim>;
    using LoadBalancer_t      = LoadBalancer<T, Dim>;
    using Base                = ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>;

    using CoordinateSelector_t = typename ParticleBinning::CoordinateSelector<ParticleContainer_t>;
    using GammaSelector_t      = typename ParticleBinning::GammaSelector<ParticleContainer_t>;
    using AdaptBins_t          = typename ParticleBinning::AdaptBinsBase<ParticleContainer_t>;
    using binIndex_t           = typename ParticleContainer_t::bin_index_type;

    using BCHandler_t          = BCHandler<Dim>;

public:
    // --- Shared state (all containers / mesh) ---

    double dt_m;                         ///< Global time step @f$\Delta t@f$ (s).
    int it_m;                            ///< Iteration counter (legacy / diagnostics).
    std::string integration_method_m;    ///< Integrator name (e.g. leapfrog).
    std::string solver_m;                ///< Field solver type string from input.
    Vector_t<int, Dim> nr_m;             ///< Mesh cell count per dimension.
    int nrZBase_m = 0;                   ///< Base z grid count before any image-charge doubling; used to reset nr_m.
    Vector_t<double, Dim> origin_m;      ///< Mesh origin (lab coordinates).
    Vector_t<double, Dim> rmin_m;        ///< Current bunch spatial minimum (from primary container stats; see calcBeamParameters).
    Vector_t<double, Dim> rmax_m;        ///< Current bunch spatial maximum (from primary container stats).
    Vector_t<double, Dim> hr_m;          ///< Mesh spacing (m).

    double lbt_m;                       ///< Load-balancer timescale parameter.
    bool isFirstRepartition_m;          ///< True until the first ORB-style repartition completes.
    ippl::NDIndex<Dim> domain_m;        ///< Global mesh index extent per dimension.
    std::array<bool, Dim> decomp_m;     ///< Domain decomposition flags (per axis).

private:
    std::vector<bool> pcActive_m;        ///< Per-container: participate in this track segment.
    std::vector<bool> pcAtZStop_m;       ///< Per-container: frozen at current z-stop until next segment.
    std::vector<std::string> particleNames_m; ///< Per-container beam particle names.

    std::vector<double> qi_m;             ///< Charge per macroparticle [C], one entry per container.
    std::vector<double> mi_m;             ///< Mass per macroparticle [GeV], one entry per container.

    const PartData* reference_m = nullptr; ///< Reference particle data (set by TrackRun::execute).

    std::shared_ptr<BunchStateHandler> bunchState_m; ///< Bunch state: unitless flag, repartition flag, etc.

    std::shared_ptr<BCHandler_t> bcHandler_m;           ///< Field boundary conditions.
    std::shared_ptr<AdaptBins_t> bins_m;                ///< Adaptive velocity/gamma binning (optional).
    FieldSolverCmd* OPALFieldSolver_m;                  ///< Borrowed parsed FIELD_SOLVER command.
    DataSink* dataSink_m;                               ///< Borrowed diagnostics and dump output sink.

    double t_m;                          ///< Current simulation time (s).

    /** Scratch E field for binned accumulation (same layout as mesh E). */
    std::shared_ptr<VField_t<T, Dim>> Etmp_m;
    /** Scratch B field for binned accumulation (same layout as mesh B). */
    std::shared_ptr<VField_t<T, Dim>> Btmp_m;

    long long globalTrackStep_m;         ///< Global integration step counter.

    std::unique_ptr<size_t[]> globalPartPerNode_m; ///< Per-rank particle counts for load-balance stats.

    double rmsDensity_m;                 ///< Legacy RMS density placeholder (may still appear in stats).

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
    PartBunch(std::vector<double> qi, 
              std::vector<double> mi,
              const std::vector<Beam*>& beams,
              std::vector<size_t> totalParticlesPerBeam,
              double lbt,
              std::string integration_method,
              FieldSolverCmd* OPALFieldSolver,
              DataSink* dataSink);

    /**
     * @brief Refresh mesh from particle extents, update layouts, and recompute moments.
     *
     * @par Typical call sites
     * - @c Track/TrackRun.cpp (initial layout)\n
     * - @c PartBunch.cpp (space-charge path; review if still appropriate)\n
     * - @c ParallelTracker.cpp (after position push / frame changes)
     */
    void bunchUpdate();

    /**
     * @brief Computes the spatial bounds for the field solver based on the current particle
     * distribution.
     *
     * Determines the minimum and maximum coordinates (`lower`, `upper`) in each dimension that
     * encompass all particles in the current container. Adjusts these bounds if image-charge
     * boundary conditions are enabled to ensure the domain includes both the real and mirrored
     * charge distributions. Guarantees a minimal span (e.g., 1e-6) in each dimension for validity
     * and applies an additional extension based on the field solver's box increment percentage.
     *
     * @param[out] lower The lowest coordinate per dimension after considering all particles and any
     *                   boundary extensions.
     * @param[out] upper The highest coordinate per dimension after considering all particles and
     * any boundary extensions.
     */
    void computeBoundsForFieldSolve(Vector_t<double, Dim>& lower, Vector_t<double, Dim>& upper);

    /**
     * @brief Updates the mesh/grid and internal data structures to match the given spatial bounds.
     *
     * Sets the mesh spacing and origin for the field container based on the difference between
     * `lower` and `upper`, updates references to domain boundaries, and applies these updates to
     * the underlying mesh and field layout. Triggers a reevaluation of particle container layouts
     * to ensure the grid matches the computed domain.
     *
     * @param[in] lower The minimum coordinates (origin) for the domain in all dimensions.
     * @param[in] upper The maximum coordinates for the domain in all dimensions.
     */
    void applyGridUpdate(
            const Vector_t<double, Dim>& lower, const Vector_t<double, Dim>& upper);

    /**
     * @brief Reinitialize the z dimension of the field grid to `nrZ` cells.
     *
     * Rebuilds the FieldLayout, all field arrays, the accumulation buffers, and the
     * IPPL Poisson solver to match the new z extent. A no-op if `nrZ` equals the
     * current z cell count. Called from `bunchUpdate` to double the z resolution
     * while image charges are active.
     *
     * @param nrZ Target number of z grid cells.
     */
    void reinitializeGridZ(int nrZ);

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

    /// @brief Set the maximum number of timesteps for which image charges are active (0 = unlimited).
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

    /// @brief Create adaptive bins from the binning command (VELOCITYZ / GAMMAZ).
    void setBins();

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
    bool isPcActive(size_t i) const {
        return i < pcActive_m.size() && pcActive_m[i];
    }

    /// @brief Force container @p i active (e.g. for containers with pending emission).
    void setPcActive(size_t i) {
        if (i < pcActive_m.size()) {
            pcActive_m[i] = true;
        }
    }

    /// @param i Container index.
    /// @return Whether container @p i is frozen at the current z-stop.
    bool pcAtZStop(size_t i) const {
        return i < pcAtZStop_m.size() && pcAtZStop_m[i];
    }

    /// @brief Deactivate container @p i until the next step-size segment (z-stop reached).
    void setPcAtZStop(size_t i);

    /// @brief After emission: reactivate non-empty containers not marked at z-stop.
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
        return this->Etmp_m; 
    }

    /// @param Etmp Scratch E field matching the mesh layout.
    void setTempEField(std::shared_ptr<VField_t<T, Dim>> Etmp) { 
        this->Etmp_m = Etmp; 
    }

    /// @brief Scratch B field used by the binned solver path.
    std::shared_ptr<VField_t<T, Dim>> getTempBField() { 
        return this->Btmp_m; 
    }

    /// @param Btmp Scratch B field matching the mesh layout.
    void setTempBField(std::shared_ptr<VField_t<T, Dim>> Btmp) { 
        this->Btmp_m = Btmp; 
    }

    /// @brief Non-const access to adaptive binning state.
    std::shared_ptr<AdaptBins_t> getBins() { 
        return bins_m; 
    }

    /// @brief Const access to adaptive binning state.
    std::shared_ptr<AdaptBins_t> getBins() const { 
        return bins_m; 
    }
    
    /// @param bins Adaptive binning object (or nullptr to clear).
    void setBins(std::shared_ptr<AdaptBins_t> bins) {
        bins_m = bins;
    }

    /// @param bcHandler Boundary-condition handler for the mesh.
    void setBCHandler(std::shared_ptr<BCHandler_t> bcHandler) {
        bcHandler_m = bcHandler;
    }

    /// @brief Current boundary-condition handler.
    std::shared_ptr<BCHandler_t> getBCHandler() const {
        return bcHandler_m;
    }

    std::shared_ptr<DataSink> getDataSink() const { return dataSink_m; }

    std::shared_ptr<BunchStateHandler> getBunchStateHandler() { return bunchState_m; }
    std::shared_ptr<const BunchStateHandler> getBunchStateHandler() const { return bunchState_m; }

    void updateMoments(){
        this->pcontainer_m->updateMoments();
    }

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
            throw OpalException("PartBunch::setCharge",
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
            throw OpalException("PartBunch::setMass",
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
            throw OpalException("PartBunch::getCharge",
                                "Number of particle containers and qi values do not match.");
        }
        if (containerIndex >= containers.size()) {
            throw OpalException("PartBunch::getCharge",
                                "Container index out of range.");
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
            throw OpalException("PartBunch::getChargePerParticle",
                                "Container index out of range.");
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
            throw OpalException("PartBunch::getMassPerParticle",
                                "Container index out of range.");
        }
        return mi_m[containerIndex];
    }

    /**
     * @brief Alias for `getCharge(containerIndex)`.
     * @param containerIndex Index of the particle container.
     * @returns Equivalent to `getCharge(containerIndex)`.
     */
    double getQ(size_t containerIndex = 0) const {
        return this->getCharge(containerIndex);
    }

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
            throw OpalException("PartBunch::getM",
                                "Number of particle containers and mi values do not match.");
        }
        if (containerIndex >= containers.size()) {
            throw OpalException("PartBunch::getM",
                                "Container index out of range.");
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
            throw OpalException("PartBunch::getTotalCharge",
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
            throw OpalException("PartBunch::getTotalMass",
                                "Number of particle containers and mi values do not match.");
        }
        double mass = 0.0;
        for (size_t i = 0; i < containers.size(); ++i) {
            mass += mi_m[i] * containers[i]->getTotalNum();
        }
        return mass;
    }

    double getdE() {
        return this->pcontainer_m->getStdKineticEnergy(); // Unit: MeV
    }
  
    const PartData* getReference() const {
        return reference_m;
    }


    /// Set inside TrackRun::execute
    void setReference (const PartData* ref) {
        reference_m = ref;
        if (reference_m && this->pcontainer_m) {
            // Ensure mean/std kinetic energy in DistributionMoments are computed using reference mass.
            // PartData mass is stored in eV; DistributionMoments expects GeV for its energy computation.
            this->pcontainer_m->setEnergyReferenceMass(reference_m->getM() * Units::eV2GeV, true);
        }
    }

    void gatherLoadBalanceStatistics();

    size_t getLoadBalance(int p) {
        return globalPartPerNode_m[p];
    }

    /**
     * @brief Transform particle positions to a unitless coordinate system.
     *
     * This converts the stored particle positions \f$ R \f$ to unitless positions
     * \f$ R' \f$ according to
     * \f[
     *   R' = \frac{R}{c \, \Delta t},
     * \f]
     * where \f$ c \f$ is the speed of light and \f$ \Delta t \f$ is a time step.
     * The resulting coordinates are dimensionless and are used by algorithms
     * that operate in this normalized coordinate system.
     *
     * By default, a single global time step \f$ \Delta t = \text{getdT()} \f$ is
     * used for all particles. If @p use_dt_per_particle is set to @c true,
     * each particle's individual time step @c dt is used instead, i.e.
     * \f$ R'_i = R_i / (c \, dt_i) \f$.
     *
     * @param use_dt_per_particle If @c true, use each particle's own @c dt
     *        value in the normalization; if @c false (default), use the
     *        global time step returned by getdT().
     *
     * @pre The bunch must not already be in unitless positions. If the internal
     *      state flag indicates that positions are already unitless,
     *      this function throws an OpalException.
     */
    void switchToUnitlessPositions(bool use_dt_per_particle = false) {
        if (bunchState_m->isUnitlessPositions()) {
            throw OpalException("PartBunch::switchToUnitlessPositions",
                                "PartBunch is already in unitless positions!");
        }

        double unitless_factor = 1.0 / (Physics::c * this->getdT());
        auto Rview             = this->getParticleContainer()->R.getView();
        auto dtview            = this->getParticleContainer()->dt.getView();
        const size_t nLocal    = this->getParticleContainer()->getLocalNum();
        Kokkos::parallel_for(
            "switchToUnitlessPositions", nLocal, KOKKOS_LAMBDA(const size_t i) {
                double fac =
                    use_dt_per_particle ? (1.0 / (Physics::c * dtview(i))) : unitless_factor;
                Rview(i) *= fac;
            });
        bunchState_m->setUnitlessPositions(true);
        *gmsg << level4 << "* Switched to unitless positions." << endl; 
    }
    /**
     * @brief Convert particle positions from unitless back to physical coordinates.
     *
     * This function undoes the transformation applied by switchToUnitlessPositions()
     * by converting positions R' in unitless coordinates back to physical positions R
     * using the relation
     *
     *   R = R' * c * dt ,
     *
     * where c is the speed of light and dt is the time step.
     *
     * If @p use_dt_per_particle is false (default), the global time step returned by
     * getdT() is used for all particles. If it is true, the conversion uses the
     * per-particle time step stored in the particle container's dt field, i.e.
     * R_i = R'_i * c * dt_i for each particle i.
     *
     * @param use_dt_per_particle Select whether to use the per-particle dt (true) or
     *                            the global dt from getdT() (false) for the scaling.
     *
     * @pre The PartBunch must currently be in unitless coordinates. If the bunch is
     *      already in physical coordinates (unitlessPositions is false), this function
     *      throws an OpalException.
     */
    void switchOffUnitlessPositions(bool use_dt_per_particle = false) {
        if (!bunchState_m->isUnitlessPositions()) {
            throw OpalException("PartBunch::switchOffUnitlessPositions",
                                "PartBunch is already in physical positions!");
        }

        double unitless_factor = Physics::c * this->getdT();
        auto Rview  = this->getParticleContainer()->R.getView();
        auto dtview = this->getParticleContainer()->dt.getView();
        const size_t nLocal    = this->getParticleContainer()->getLocalNum();
        Kokkos::parallel_for(
            "switchOffUnitlessPositions", nLocal,
            KOKKOS_LAMBDA(const size_t i) {
                double fac = use_dt_per_particle ? (Physics::c * dtview(i)) : unitless_factor;
                Rview(i) *= fac;
            });
        bunchState_m->setUnitlessPositions(false);
        *gmsg << level4 << "* Switched to physical positions." << endl;
    }

    /// @brief Stub; logs and returns 0.
    size_t calcNumPartsOutside(Vector_t<double, Dim> /*x*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0;
    }

    /// @brief Stub; logs only.
    void calcLineDensity(
        unsigned int /*nBins*/, std::vector<double>& /*lineDensity*/, std::pair<double, double>& /*meshInfo*/) {
            *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }

    /// @brief Stub; logs and returns zero vector.
    Vector_t<double, Dim> getEExtrema() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
       return Vector_t<double, Dim>(0);
    }

    /**
     * @brief Compute the bunch self-fields (binned when available).
     *
     * The actual implementation lives in the solver object (see `BinnedFieldSolver`).
     * `ParallelTracker` only orchestrates reference/beam-frame transforms and calls
     * this delegator once per step.
     */
    void computeSelfFields();

    /**
     * @brief Write bin edges/counts to the data sink when configured.
     * @param preMerge True if called before a bin-merge step.
     */
    void dumpBinConfig(bool preMerge);

    /// @brief Human-readable dump of each container to @p os.
    /// @param os Output stream wrapper.
    /// @return Reference to @p os.
    Inform& print(Inform& os);

    /// @return True if a field solver instance is installed.
    bool hasFieldSolver() const {
        return this->fsolver_m != nullptr;
    }

    /**
     * @brief Non-const pointer to the concrete BinnedFieldSolver.
     * @note Definition in PartBunch.cpp; requires complete BinnedFieldSolver type.
     */
    BinnedFieldSolver_t* getFieldSolver();

    /// @brief Const overload of getFieldSolver().
    const BinnedFieldSolver_t* getFieldSolver() const;

    /// @brief Backend type string (e.g. FFT, OPEN, CG, NONE).
    std::string getFieldSolverType();

    /// @return True if adaptive binning is configured.
    bool hasBinning() const {
        return this->bins_m != nullptr;
    }

    /**
     * @brief Effective bin count for diagnostics (1 if binning inactive or still at max bins).
     */
    int getCurrentNBins() const {
        if (!hasBinning()) {
            return 1;
        }

        int ret_bins = static_cast<int>(bins_m->getCurrentBinCount());
        // If the number of bins is the same as the maximum number of bins, we haven't merged bins
        // yet (likely because the simulation is too empty)
        if (ret_bins == this->getBins()->getMaxBinCount()) {
            Inform m("PartBunch::getCurrentNBins");
            m << level4
              << "WARNING: Number of bins is the same as the maximum number of bins, we haven't "
                 "merged bins yet (likely because the simulation is too empty). Returning 1. If "
                 "that is not the case, check e.g. binning parameters."
              << endl;
            return 1;
        } else {
            return ret_bins;
        }
    }

    /// @brief Compatibility stub; logs and returns 0.
    double calcMeanPhi() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
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
    void setdT(double dt) {
        dt_m = dt;
    }

    /// @brief Get the global time step.
    double getdT() const {
        return dt_m;
    }

    /// @brief Set the current simulation time.
    void setT(double t) {
        t_m = t;
    }

    /// @brief Advance time by one global time step.
    void incrementT() {
        t_m += dt_m;
    }

    /// @brief Get the current simulation time.
    double getT() const {
        return t_m;
    }


    /// @brief Cached minimum extent (@c rmin_m); prefer per-container min/max for multi-beam detail.
    Vector_t<double, Dim> get_origin() const {
        return rmin_m;
    }
    /// @brief Cached maximum extent (@c rmax_m).
    Vector_t<double, Dim> get_maxExtent() const {
        return rmax_m;
    }

    /// @brief Stub; logs and returns zero vector.
    Vector_t<double, Dim> get_halo() const {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }

    /// @brief Get mesh spacing
    Vector_t<double, Dim> get_hr() const {
        return hr_m;
    }

    /// @param n Global step counter to store.
    void setGlobalTrackStep(long long n) {
        globalTrackStep_m = n;
    }

    /// @brief Current global tracking step.
    long long getGlobalTrackStep() const {
        return globalTrackStep_m;
    }

    /// @brief Increment @c globalTrackStep_m by one.
    void incTrackSteps() {
        globalTrackStep_m++;
    }

    /// @brief ORB/binary repartition when the load balancer requests it (primary container).
    void do_binaryRepart();

    /// @brief Legacy RMS density field (may be unused).
    double get_rmsDensity() const {
        return rmsDensity_m;
    }

};

extern template class PartBunch<double, 3>; ///< Explicit instantiation for 3D double-precision.

#endif  // PARTBUNCH_H
