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

class DataSink;  // forward declaration; full type needed only in .cpp

extern Inform* gmsg;

using view_type = typename ippl::detail::ViewType<ippl::Vector<double, 3>, 1>::view_type;

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
    // Per container values ====================================================
    
    // Reference particle values
    Vector_t<double, Dim> RefPartR_m; // reference particle position
    Vector_t<double, Dim> RefPartP_m; // reference particle momentum
    CoordinateSystemTrafo toLabTrafo_m; // transformation to lab frame

    // Shared values for all containers ========================================
    
    double dt_m; // time step
    int it_m; // iteration count

    double lbt_m; // load balancer threshold
    // Routed through bunchState_m.isFirstRepartition() / setFirstRepartition().
    ippl::NDIndex<Dim> domain_m;
    std::array<bool, Dim> decomp_m;

    // Solver
    std::string integration_method_m; // integration method
    std::string solver_m; // field solver type

    // Mesh
    Vector_t<int, Dim> nr_m; // number of grid points
    int nrZBase_m = 0;       // base z grid count (before any image-charge doubling). needed to reset nr_m later
    Vector_t<double, Dim> origin_m; // origin of the mesh
    Vector_t<double, Dim> rmin_m; // minimum extent of the mesh
    Vector_t<double, Dim> rmax_m; // maximum extent of the mesh
    Vector_t<double, Dim> hr_m; // mesh size [m]

    // Unused values ===========================================================
   

private:

    // Per container values ====================================================

    std::vector<double> qi_m; // charge per macroparticle [C], per container
    std::vector<double> mi_m; // mass per macroparticle [GeV], per container

    std::unique_ptr<size_t[]> globalPartPerNode_m; // reducer object for load balance statistics

    Vector_t<double, Dim> globalMeanR_m; // global mean position
    Quaternion_t globalToLocalQuaternion_m; // global to local quaternion
    
    const PartData* reference_m; // reference particle data
   
    double spos_m; // s position along design trajectory


    // Shared values for all containers ========================================

    std::shared_ptr<BCHandler_t> bcHandler_m; // field boundary handler
    std::shared_ptr<AdaptBins_t> bins_m; // adaptive binning structure

    std::shared_ptr<FieldSolverCmd> OPALFieldSolver_m; // field solver command
    std::shared_ptr<DataSink> dataSink_m; // data sink

    std::shared_ptr<BunchStateHandler> bunchState_m;

    double t_m; // time of integration

    /// Temporary E field container used to store temporary E field during binned solver
    std::shared_ptr<VField_t<T, Dim>> Etmp_m;
    /// Temporary B field container used to store temporary B field during binned solver
    std::shared_ptr<VField_t<T, Dim>> Btmp_m;

    long long globalTrackStep_m;
    // Unused values ===========================================================

    double rmsDensity_m;


public:

    /**
     * @brief Construct a PartBunch with given macro charge/mass and configuration.
     *
     * @param qi              Vector of macrocharges per species [C]
     * @param mi              Vector of macromasses per species [GeV/c^2]
     * @param num_containeres Total number of containers 
     * @param lbt             Load-balancer timescale.
     * @param integration_method Name of the integrator (e.g. "LF2").
     * @param OPALFieldSolver Field solver command providing mesh and binning configuration.
     * @param dataSink        Shared pointer to the global DataSink used for diagnostics.
     */
    PartBunch(std::vector<double> qi, 
              std::vector<double> mi,
              size_t num_containers,
              double lbt,
              std::string integration_method,
              std::shared_ptr<FieldSolverCmd> OPALFieldSolver,
              std::shared_ptr<DataSink> dataSink);
    /**
     * @brief 
     * - recomputes mesh spacing i.e. Layout
     * - repatition the domain if nessesary
     * - recomputes all moments of the particle distribution
     
     called in:
     
     Track/TrackRun.cpp             --> initial calc)
     PartBunch/PartBunch.cpp        --> in space charge which is wrong 
     Algorithms/ParallelTracker.cpp --> after push
     
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
     * @brief Configure diagnostic dump frequency for the ZEROFACE plane potential.
     *
     * @param frequency Dump every n-th global timestep. `0` disables dumping.
     */
    void setZeroFacePlaneDumpFrequency(int frequency);

    /// @brief Set the maximum number of timesteps for which image charges are active (0 = unlimited).
    void setZerofaceMaxSteps(int maxSteps);

    size_t getTotalNumAllContainers() const {
        size_t total = 0;
        for (const auto& pc : this->getParticleContainers()) {
            if (pc) {
                total += pc->getTotalNum();
            }
        }
        return total;
    }

    void setSolver();

    void setBins();

    void pre_run() override;

    void performBunchSanityChecks() const;

    void advance() override {
        throw OpalException(
            "PartBunch::advance",
            "Not used: just exists because ippl::PicManager wants it that way.");
    }
    void par2grid() override {
        throw OpalException(
            "PartBunch::par2grid",
            "Not used: just exists because ippl::PicManager wants it that way.");
    }
    void grid2par() override {
        throw OpalException(
            "PartBunch::grid2par",
            "Not used: just exists because ippl::PicManager wants it that way.");
    }

public:
    std::shared_ptr<VField_t<T, Dim>> getTempEField() { return this->Etmp_m; }
    void setTempEField(std::shared_ptr<VField_t<T, Dim>> Etmp) { this->Etmp_m = Etmp; }
    std::shared_ptr<VField_t<T, Dim>> getTempBField() { return this->Btmp_m; }
    void setTempBField(std::shared_ptr<VField_t<T, Dim>> Btmp) { this->Btmp_m = Btmp; }

    std::shared_ptr<AdaptBins_t> getBins() { return bins_m; }
    std::shared_ptr<AdaptBins_t> getBins() const { return bins_m; }
    
    void setBins(std::shared_ptr<AdaptBins_t> bins) { bins_m = bins; }

    void setBCHandler(std::shared_ptr<BCHandler_t> bcHandler) { bcHandler_m = bcHandler; }
    std::shared_ptr<BCHandler_t> getBCHandler() const { return bcHandler_m; }
    std::shared_ptr<DataSink> getDataSink() const { return dataSink_m; }

    std::shared_ptr<BunchStateHandler> getBunchStateHandler() { return bunchState_m; }
    std::shared_ptr<const BunchStateHandler> getBunchStateHandler() const { return bunchState_m; }

    void updateMoments(){
        this->pcontainer_m->updateMoments();
    }

    size_t getTotalNum() const {
        return this->pcontainer_m->getTotalNum();
    }

    size_t getLocalNum() const {
        return this->pcontainer_m->getLocalNum();
    }

    void calcBeamParameters();

    void do_binaryRepart();

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
        const size_t nLocal    = this->getLocalNum();
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
        const size_t nLocal    = this->getLocalNum();
        Kokkos::parallel_for(
            "switchOffUnitlessPositions", nLocal,
            KOKKOS_LAMBDA(const size_t i) {
                double fac = use_dt_per_particle ? (Physics::c * dtview(i)) : unitless_factor;
                Rview(i) *= fac;
            });
        bunchState_m->setUnitlessPositions(false);
        *gmsg << level4 << "* Switched to physical positions." << endl;
    }

    size_t calcNumPartsOutside(Vector_t<double, Dim> /*x*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0;
    }

    void calcLineDensity(
        unsigned int /*nBins*/, std::vector<double>& /*lineDensity*/, std::pair<double, double>& /*meshInfo*/) {
            *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }

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
    void dumpBinConfig(bool preMerge);

    Inform& print(Inform& os);

    bool hasFieldSolver() const {
        return this->fsolver_m != nullptr;
    }

    BinnedFieldSolver_t* getFieldSolver() {
        return dynamic_cast<BinnedFieldSolver_t*>(this->fsolver_m.get());
    }

    /// Const overload for better const correctness.
    const BinnedFieldSolver_t* getFieldSolver() const {
        return dynamic_cast<const BinnedFieldSolver_t*>(this->fsolver_m.get());
    }

    std::string getFieldSolverType();

    bool hasBinning() const {
        return this->bins_m != nullptr;
    }

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

    double calcMeanPhi() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }

    void get_bounds(Vector_t<double, Dim>& rmin, Vector_t<double, Dim>& rmax) {
        rmin = rmin_m;
        rmax = rmax_m;
    }

    Vector_t<double, Dim> R(size_t) {
        throw OpalException(
            "PartBunch::R",
            "Not implemented: shouldn't be called, since this is not the correct way to access "
            "particle positions.");
        return Vector_t<double, Dim>(0.0);
    }

    Vector_t<double, Dim> P(size_t) {
        throw OpalException(
            "PartBunch::P",
            "Not implemented: shouldn't be called, since this is not the correct way to access "
            "particle momenta.");
        return Vector_t<double, Dim>(0.0);
    }

    void setdT(double dt) {
        dt_m = dt;
    }

    double getdT() const {
        return dt_m;
    }

    void setT(double t) {
        t_m = t;
    }

    void incrementT() {
        t_m += dt_m;
    }

    double getT() const {
        return t_m;
    }

    void set_sPos(double s) {
        spos_m = s;
    }

    double get_sPos() const {
        return spos_m;
    }

    double get_gamma() const {
        return this->pcontainer_m->getMeanGammaZ();
    }

    /// Mean kinetic energy over particles (mean of per-particle kinetic energy), in MeV.
    double get_meanKineticEnergy();

    Vector_t<double, Dim> get_origin() const {
        return rmin_m;
    }
    Vector_t<double, Dim> get_maxExtent() const {
        return rmax_m;
    }

    // \todo in opal, MeanPosition is return for get_centroid, which I think is wrong. We already have get_rmean()
    Vector_t<double, 2*Dim> get_centroid() const {
        return this->pcontainer_m->getCentroid();
    }

    Vector_t<double, Dim> get_rrms() const {
        return this->pcontainer_m->getRmsR();
    }

    Vector_t<double, Dim> get_rprms() const {
        return this->pcontainer_m->getRmsRP();
    }

    Vector_t<double, Dim> get_prms() const {
        return this->pcontainer_m->getRmsP();
    }

    Vector_t<double, Dim> get_rmean() const {
        return this->pcontainer_m->getMeanR();
    }

    Vector_t<double, Dim> get_pmean() const {
        return this->pcontainer_m->getMeanP();
    }

    Vector_t<double, Dim> get_emit() const {
        return this->pcontainer_m->getGeometricEmit();
    }
    Vector_t<double, Dim> get_norm_emit() const {
        return this->pcontainer_m->getNormEmit();
    }

  
    Vector_t<double, Dim> get_halo() const {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }
    // Not used, but might be useful later. Commented out for now.
    /*Vector_t<double, Dim> get_68Percentile() const {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_95Percentile() const {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_99Percentile() const {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_99_99Percentile() const {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_normalizedEps_68Percentile() const {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_normalizedEps_95Percentile() const {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_normalizedEps_99Percentile() const {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_normalizedEps_99_99Percentile() const {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }*/

    Vector_t<double, Dim> get_hr() const {
        return hr_m;
    }

    double get_Dx() const {
        return this->pcontainer_m->getDx();
    }
    double get_Dy() const {
        return this->pcontainer_m->getDy();
    }
    double get_DDx() const {
        return this->pcontainer_m->getDDx();
    }
    double get_DDy() const {
        return this->pcontainer_m->getDDy();
    }

    double get_temperature() const {
        return this->pcontainer_m->getTemperature();
    }

    void calcDebyeLength() {
         this->pcontainer_m->computeDebyeLength(rmsDensity_m);
    }

    double get_debyeLength() const {
        return this->pcontainer_m->getDebyeLength();
    }

    double get_plasmaParameter() const {
        return this->pcontainer_m->getPlasmaParameter();
    }

    double get_rmsDensity() const {
        return rmsDensity_m;
    }

    /// step in multiple TRACK commands
    void setGlobalTrackStep(long long n) {
        globalTrackStep_m = n;
    }

    long long getGlobalTrackStep() const {
        return globalTrackStep_m;
    }

    void incTrackSteps() {
        globalTrackStep_m++;
    }

    void setGlobalMeanR(Vector_t<double, Dim> globalMeanR) {
        globalMeanR_m = globalMeanR;
    }

    Vector_t<double, Dim> getGlobalMeanR() {
        return globalMeanR_m;
    }

    void setGlobalToLocalQuaternion(Quaternion_t globalToLocalQuaternion) {
        globalToLocalQuaternion_m = globalToLocalQuaternion;
    }

    Quaternion_t getGlobalToLocalQuaternion() {
        return globalToLocalQuaternion_m;
    }
};

// Explicit instantiations
extern template class PartBunch<double, 3>;

#endif
