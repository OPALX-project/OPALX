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
class PartBunch : public ippl::PicManager<
          T, Dim, ParticleContainer<T, Dim>, FieldContainer<T, Dim>, LoadBalancer<T, Dim>> {
public:
    using ParticleContainer_t = ParticleContainer<T, Dim>;
    using FieldContainer_t    = FieldContainer<T, Dim>;
    using FieldSolver_t       = FieldSolver<T, Dim>;
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
     * @param dataSink               Borrowed diagnostics output sink.
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

    /**
     * @brief Update moments and set @c rmin_m / @c rmax_m from the primary (first) container.
     * @note Space-charge and mesh updates currently use container 0 only.
     */
    void calcBeamParameters();

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
     * @brief Non-const pointer to the concrete field solver.
     * @todo Prefer smart pointers once IPPL @c FieldSolverBase allows it.
     */
    FieldSolver_t* getFieldSolver() {
        return static_cast<FieldSolver_t*>(this->fsolver_m.get());
    }

    /// @brief Const overload of getFieldSolver().
    const FieldSolver_t* getFieldSolver() const {
        return static_cast<const FieldSolver_t*>(this->fsolver_m.get());
    }

    /// @brief Backend type string (e.g. FFT, OPEN, CG, NONE).
    std::string getFieldSolverType() {
        return this->getFieldSolver()->getStype();
    }

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

    /// @brief All-reduce local particle counts into @c globalPartPerNode_m.
    void gatherLoadBalanceStatistics();

    /// @param p MPI rank index.
    /// @return Particle count recorded for rank @p p after the last gather.
    size_t getLoadBalance(int p) {
        return globalPartPerNode_m[p];
    }

    /// @brief Legacy RMS density field (may be unused).
    double get_rmsDensity() const {
        return rmsDensity_m;
    }

};

extern template class PartBunch<double, 3>; ///< Explicit instantiation for 3D double-precision.

#endif  // PARTBUNCH_H
