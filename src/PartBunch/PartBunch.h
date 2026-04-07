#ifndef PARTBUNCH_H
#define PARTBUNCH_H

#include <memory>
#include <vector>

#include "Manager/BaseManager.h"
#include "Manager/PicManager.h"
#include "PartBunch/FieldContainer.hpp"
#include "PartBunch/FieldSolver.hpp"
#include "PartBunch/LoadBalancer.hpp"
#include "PartBunch/ParticleContainer.hpp"
#include "PartBunch/Binning/AdaptBins.h"
#include "PartBunch/Binning/AdaptBins.tpp"

#include "Utilities/OpalException.h"
#include "BCHandler.hpp"
#include "Structure/FieldSolverCmd.h"

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

    /// Per-container dynamics on this segment.
    std::vector<bool> pcActive_m;
    /// Per-container: segment z-stop reached (no reactivation after emit).
    std::vector<bool> pcAtZStop_m;

// Shared values for all containers ============================================

    // Field boundary handler
    std::shared_ptr<BCHandler_t> bcHandler_m; 

    // Adaptive binning structure
    std::shared_ptr<AdaptBins_t> bins_m; 

    // Field solver command
    std::shared_ptr<FieldSolverCmd> OPALFieldSolver_m; 

    // Data sink
    std::shared_ptr<DataSink> dataSink_m;

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

    /// Segment start: non-empty containers active; empty inactive.
    void resetPcActive();

    bool isPcActive(size_t i) const {
        return i < pcActive_m.size() && pcActive_m[i];
    }

    bool pcAtZStop(size_t i) const {
        return i < pcAtZStop_m.size() && pcAtZStop_m[i];
    }

    /// Freeze dynamics for this container until next segment.
    void setPcAtZStop(size_t i);

    /// After emit: activate non-empty containers not frozen at z-stop.
    void refreshPcActiveAfterEmit();

    bool anyPcActive() const {
        for (bool a : pcActive_m) {
            if (a) {
                return true;
            }
        }
        return false;
    }

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
    /// @note For now this only considers container 0, which is fine since
    /// the field solver acts on container 0 for now. 
    void calcBeamParameters();

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
