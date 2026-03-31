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

class DataSink;  // forward declaration; full type needed only in .cpp
class Beam;

extern Inform* gmsg;

using view_type = typename ippl::detail::ViewType<ippl::Vector<double, 3>, 1>::view_type;

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
// Per container values ========================================================

// Moved to ParticleContainer for now

// Shared values for all containers ============================================
    
    // Time step
    double dt_m; 

    // Iteration count
    int it_m; 

    // Integration methods : LEAPFROG 
    std::string integration_method_m; 

    // Field solver type
    std::string solver_m; 

    // Number of mesh points per dimension
    Vector_t<int, Dim> nr_m; 

    // Mesh spacial origin
    Vector_t<double, Dim> origin_m; 

    // Maximum point of the mesh 
    Vector_t<double, Dim> rmin_m; 

    // Minimum point of the mesh
    Vector_t<double, Dim> rmax_m; 

    // Mesh spacing [m]
    Vector_t<double, Dim> hr_m; 

// Load Balancing (To be properly implemented) =================================

    // Load balancing threshold
    double lbt_m; 

    // Is first repartition flat
    bool isFirstRepartition_m; 

    // 3-D array for the domain
    ippl::NDIndex<Dim> domain_m;

    // Boolean array 
    std::array<bool, Dim> decomp_m;

private:

// Per container values ========================================================

// Moved to ParticleContainer for now

// Shared values for all containers ============================================

    // Field boundary handler
    std::shared_ptr<BCHandler_t> bcHandler_m; 

    // Adaptive binning structure
    std::shared_ptr<AdaptBins_t> bins_m; 

    // Field solver command
    std::shared_ptr<FieldSolverCmd> OPALFieldSolver_m; 

    // Data sink
    std::shared_ptr<DataSink> dataSink_m;

    // unit state of PartBunch --> always false after initialization, so use this as standard flag
    bool isUnitless_m = false; 

    // Time of integration
    double t_m; 

    /// Temporary E field container used to store temporary E field during binned solver
    std::shared_ptr<VField_t<T, Dim>> Etmp_m;

    /// Temporary B field container used to store temporary B field during binned solver
    std::shared_ptr<VField_t<T, Dim>> Btmp_m;

    // Global tracking step 
    long long globalTrackStep_m;
    
// Load Balancing (To be properly implemented) =================================

    // reducer object for load balance statistics
    std::unique_ptr<size_t[]> globalPartPerNode_m; 

// Unused values ===============================================================

    // Still written to stat file for some reason? 
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
              const std::vector<Beam*>& beams,
              std::vector<size_t> totalParticlesPerBeam,
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
     * @brief Returns the total number of particle across all containers
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

    /// @brief Set field solver
    void setSolver();

    /// @brief Set bins for binned field solver
    void setBins();

    /// @brief Setup run for the field solver
    void pre_run() override;

    /// @brief Sanity check 
    void performBunchSanityChecks() const;

    // ! NOT TO BE USED: This functionality has moved elsewhere
    void advance() override {
        throw OpalException(
            "PartBunch::advance",
            "Not used: just exists because ippl::PicManager wants it that way.");
    }

    // ! NOT TO BE USED: This functionality has moved elsewhere
    void par2grid() override {
        throw OpalException(
            "PartBunch::par2grid",
            "Not used: just exists because ippl::PicManager wants it that way.");
    }

    // ! NOT TO BE USED: This functionality has moved elsewhere
    void grid2par() override {
        throw OpalException(
            "PartBunch::grid2par",
            "Not used: just exists because ippl::PicManager wants it that way.");
    }

    /// @brief Get the temporary E field
    std::shared_ptr<VField_t<T, Dim>> getTempEField() { 
        return this->Etmp_m; 
    }

    /// @brief Set the temporary E field
    void setTempEField(std::shared_ptr<VField_t<T, Dim>> Etmp) { 
        this->Etmp_m = Etmp; 
    }

    /// @brief Get temporary B field
    std::shared_ptr<VField_t<T, Dim>> getTempBField() { 
        return this->Btmp_m; 
    }

    /// @brief Set the temporary B field
    void setTempBField(std::shared_ptr<VField_t<T, Dim>> Btmp) { 
        this->Btmp_m = Btmp; 
    }

    /// @brief Get bins
    std::shared_ptr<AdaptBins_t> getBins() { 
        return bins_m; 
    }

    /// @brief Get bins const version
    std::shared_ptr<AdaptBins_t> getBins() const { 
        return bins_m; 
    }
    
    /// @brief Set bins
    void setBins(std::shared_ptr<AdaptBins_t> bins) { 
        bins_m = bins; 
    }

    /// @brief Set BC handler
    void setBCHandler(std::shared_ptr<BCHandler_t> bcHandler) { 
        bcHandler_m = bcHandler; 
    }

    /// @brief Get BC handler
    std::shared_ptr<BCHandler_t> getBCHandler() const { 
        return bcHandler_m; 
    }

    /// @brief Update Moments and calculate rmin_m and rmax_m
    void calcBeamParameters();

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
        if (isUnitless_m) {
            throw OpalException("PartBunch::switchToUnitlessPositions",
                                "PartBunch is already in unitless positions!");
        }

        // Divide by c*dt
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
        isUnitless_m = true;

        /// \todo remove later
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
     *      already in physical coordinates (isUnitless_m is false), this function
     *      throws an OpalException.
     */
    void switchOffUnitlessPositions(bool use_dt_per_particle = false) {
        if (!isUnitless_m) {
            throw OpalException("PartBunch::switchOffUnitlessPositions",
                                "PartBunch is already in physical positions!");
        }

        // Multiply by c*dt
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
        isUnitless_m = false;

        /// \todo remove later
        *gmsg << level4 << "* Switched to physical positions." << endl;
    }

    // ! NOT IMPLEMENTED
    size_t calcNumPartsOutside(Vector_t<double, Dim> /*x*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0;
    }

    // ! NOT IMPLEMENTED
    void calcLineDensity(
        unsigned int /*nBins*/, std::vector<double>& /*lineDensity*/, std::pair<double, double>& /*meshInfo*/) {
            *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }

    // ! NOT IMPLEMENTED
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

    /// @brief Dump binning config
    void dumpBinConfig(bool preMerge);

    /// @brief Print information for each bunch (container)
    Inform& print(Inform& os);

    /// @brief Checks existence of field solver
    bool hasFieldSolver() const {
        return this->fsolver_m != nullptr;
    }

    /// @brief Get field solver
    FieldSolver_t* getFieldSolver() {
        /*
        \todo this needs to change, best would be to use a smart pointer!
        However, the parent class FieldSolverBase from IPPL uses raw pointers,
        so changing this would require some changes in IPPL...
        */
        return static_cast<FieldSolver_t*>(this->fsolver_m.get());
    }

    /// Const overload for better const correctness.
    const FieldSolver_t* getFieldSolver() const {
        return static_cast<const FieldSolver_t*>(this->fsolver_m.get());
    }

    /// @brief Get field solver backend type string.
    std::string getFieldSolverType() {
        return this->getFieldSolver()->getStype();
    }

    /// @brief Check whether adaptive binning is enabled.
    bool hasBinning() const {
        return this->bins_m != nullptr;
    }

    /// @brief Get active number of bins (returns 1 when binning is effectively inactive).
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

    // ! NOT TO BE USED: COMPATIBILITY STUB
    double calcMeanPhi() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }


    // ! NOT TO BE USED: COMPATIBILITY STUB
    Vector_t<double, Dim> R(size_t) {
        throw OpalException(
            "PartBunch::R",
            "Not implemented: shouldn't be called, since this is not the correct way to access "
            "particle positions.");
        return Vector_t<double, Dim>(0.0);
    }

    // ! NOT TO BE USED: COMPATIBILITY STUB
    Vector_t<double, Dim> P(size_t) {
        throw OpalException(
            "PartBunch::P",
            "Not implemented: shouldn't be called, since this is not the correct way to access "
            "particle momenta.");
        return Vector_t<double, Dim>(0.0);
    }

    /// @brief Get the current local particle-position bounds.
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


    /// @brief Get the current lower bound of the bunch extent.
    Vector_t<double, Dim> get_origin() const {
        return rmin_m;
    }
    /// @brief Get the current upper bound of the bunch extent.
    Vector_t<double, Dim> get_maxExtent() const {
        return rmax_m;
    }

    // ! NOT IMPLEMENTED
    Vector_t<double, Dim> get_halo() const {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }

    /// @brief Get mesh spacing
    Vector_t<double, Dim> get_hr() const {
        return hr_m;
    }

    /// @brief Set tracking step 
    void setGlobalTrackStep(long long n) {
        globalTrackStep_m = n;
    }

    /// @brief Get tracking step
    long long getGlobalTrackStep() const {
        return globalTrackStep_m;
    }

    /// @brief Increment tracking step
    void incTrackSteps() {
        globalTrackStep_m++;
    }

// Load Balancing (To be properly implemented) =================================
   
    /// @brief Load balancing reparitioning
    void do_binaryRepart();

    void gatherLoadBalanceStatistics();

    size_t getLoadBalance(int p) {
        return globalPartPerNode_m[p];
    }

// Unused / Commented out ======================================================

    double get_rmsDensity() const {
        return rmsDensity_m;
    }

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

};

// Explicit instantiations
extern template class PartBunch<double, 3>;

#endif
