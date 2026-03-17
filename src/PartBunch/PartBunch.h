#ifndef PARTBUNCH_H
#define PARTBUNCH_H

#include <memory>

#include "Algorithms/Matrix.h"
#include "Algorithms/CoordinateSystemTrafo.h"
#include "Attributes/Attributes.h"
#include "Distribution/Distribution.h"
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

#include "Utilities/Options.h" // Needed to define binning parameters!
#include "PartBunch/Binning/AdaptBins.h" // TODO: binning

class DataSink;  // forward declaration; full type needed only in .cpp

extern Inform* gmsg;

template <typename T>
KOKKOS_INLINE_FUNCTION typename T::value_type L2Norm(T& x) {
    return sqrt(dot(x, x).apply());
}

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

    using BinningSelector_t   = typename ParticleBinning::CoordinateSelector<ParticleContainer_t>;
    using AdaptBins_t         = typename ParticleBinning::AdaptBins<ParticleContainer_t, BinningSelector_t>;
    using binIndex_t          = typename ParticleContainer_t::bin_index_type;

    using BCHandler_t         = BCHandler<Dim>;

    double time_m;

    size_type totalP_m;

    /// \todo doesn't do anything??? 
    // int nt_m; 

    double lbt_m;

    double dt_m;

    int it_m;

    std::string integration_method_m;

    std::string solver_m;

    bool isFirstRepartition_m;

private:
    double qi_m;

    double mi_m;

    double rmsDensity_m;

    std::shared_ptr<BCHandler_t> bcHandler_m;

public:
    Vector_t<int, Dim> nr_m;

    Vector_t<double, Dim> origin_m;
    Vector_t<double, Dim> rmin_m;
    Vector_t<double, Dim> rmax_m;

    /// mesh size [m]
    Vector_t<double, Dim> hr_m;

    // Landau damping specific
    double Bext_m;
    double alpha_m;
    double DrInv_m;

    ippl::NDIndex<Dim> domain_m;
    std::array<bool, Dim> decomp_m;

    /*
      Up to here it is like the opaltest
    */

    /**
      Reference particle structures
     */

    Vector_t<double, Dim> RefPartR_m;
    Vector_t<double, Dim> RefPartP_m;

    CoordinateSystemTrafo toLabTrafo_m;

private:

    std::unique_ptr<size_t[]> globalPartPerNode_m;

    // ParticleOrigin refPOrigin_m;
    // ParticleType refPType_m;

    /// Initialize the translation vector and rotation quaternion
    /// here. Cyclotron tracker will reset these values each timestep
    /// TTracker can just use 0 translation and 0 rotation (quat[1 0 0 0]).
    // Vector_t globalMeanR_m = Vector_t(0.0, 0.0, 0.0);
    // Quaternion_t globalToLocalQuaternion_m = Quaternion_t(1.0, 0.0, 0.0, 0.0);
    Vector_t<double, Dim> globalMeanR_m;
    Quaternion_t globalToLocalQuaternion_m;

    /**
       The structure for particle binning
    */

    // PartBins* pbin_m;

    /// if larger than 0, emitt particles for tEmission_m [s]
    double tEmission_m;

    /// holds the gamma of the bin
    std::unique_ptr<double[]> bingamma_m;

    // FIXME: this should go into the Bin class!
    //  holds number of emitted particles of the bin
    //  jjyang: opal-cycl use *nBin_m of pbin_m
    //std::unique_ptr<size_t[]> binemitted_m; // liemen_a: TODO remove!
    std::shared_ptr<AdaptBins_t> bins_m; // added by liemen_a for AdaptBins class!

    /// steps per turn for OPAL-cycl
    int stepsPerTurn_m;

    /// current bunch number
    short numBunch_m;

    /// number of particles per bunch
    std::vector<size_t> bunchTotalNum_m;
    std::vector<size_t> bunchLocalNum_m;

    /// this parameter records the current steps since last bunch injection
    /// it helps to inject new bunches correctly in the restart run of OPAL-cycl
    /// it is stored during phase space dump.
    int SteptoLastInj_m;

    bool fixed_grid;

    const PartData* reference_m;

    /// step in a TRACK command
    long long localTrackStep_m;

    /// if multiple TRACK commands
    long long globalTrackStep_m;


    std::shared_ptr<FieldSolverCmd> OPALFieldSolver_m;

    std::shared_ptr<DataSink> dataSink_m;
    
    // unit state of PartBunch --> always false after initialization, so use this as standard flag
    // UnitState_t unit_state_m;
    bool isUnitless_m = false;
    // UnitState_t stateOfLastBoundP_m;

    /// holds the actual time of the integration
    double t_m;

    /// the position along design trajectory
    double spos_m;

    /*
       flags to tell if we are a DC-beam
     */
    bool dcBeam_m;
    double periodLength_m;

    /// Temporary E field container used to store temporary E field during binned solver
    std::shared_ptr<VField_t<T, Dim>> Etmp_m;

public:

    /**
     * @brief Construct a PartBunch with given macro charge/mass and configuration.
     *
     * @param qi              Charge per macroparticle [C].
     * @param mi              Mass per macroparticle [GeV/c^2].
     * @param totalP          Total number of macroparticles.
     * @param lbt             Load-balancer timescale.
     * @param integration_method Name of the integrator (e.g. "LF2").
     * @param OPALFieldSolver Field solver command providing mesh and binning configuration.
     * @param dataSink        Shared pointer to the global DataSink used for diagnostics.
     */
    PartBunch(double qi,
              double mi,
              size_t totalP,
              double lbt,
              std::string integration_method,
              std::shared_ptr<FieldSolverCmd>& OPALFieldSolver,
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
  
    ~PartBunch() {
        *gmsg << level2 << "* PartBunch Destructor: Finished time step: " << this->it_m << " time: " << this->time_m << endl;
    }

    std::shared_ptr<ParticleContainer_t> getParticleContainer() {
        return this->pcontainer_m;
    }

    void setSolver();

    void setBins();

    void pre_run() override ;

    void performBunchSanityChecks() const;

public:
    std::shared_ptr<VField_t<T, Dim>> getTempEField() { return this->Etmp_m; }
    void setTempEField(std::shared_ptr<VField_t<T, Dim>> Etmp) { this->Etmp_m = Etmp; }

    std::shared_ptr<AdaptBins_t> getBins() { return bins_m; } // TODO: Binning
    
    void setBins(std::shared_ptr<AdaptBins_t> bins) { bins_m = bins; } // TODO: Binning

    void setBCHandler(std::shared_ptr<BCHandler_t> bcHandler) { bcHandler_m = bcHandler; }
    std::shared_ptr<BCHandler_t> getBCHandler() const { return bcHandler_m; }

    void updateMoments(){
        this->pcontainer_m->updateMoments();
    }

    size_t getTotalNum() const {
        return this->pcontainer_m->getTotalNum();
    }

    size_t getLocalNum() const {
        return this->pcontainer_m->getLocalNum();
    }

    Vector_t<double, Dim> R(size_t /*i*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }

    Vector_t<double, Dim> P(size_t /*i*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }

    Vector_t<double, Dim> Ef(size_t /*i*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }

    Vector_t<double, Dim> Bf(size_t /*i*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }

    Vector_t<double, Dim> dt(size_t /*i*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return Vector_t<double, Dim>(0.0);
    }

    void advance() override {
        // \todo needs to go
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }

    void par2grid() override {
        scatterCIC();
    }

    void grid2par() override {
        gatherCIC();
    }

    void gatherCIC();

    void scatterCIC() {
        scatterCICPerBin(-1);
    } 

    void scatterCICPerBin(binIndex_t binIndex);

    /*
      Up to here it is like the opaltest
    */

    T getCouplingConstant() const;

    void calcBeamParameters();

    void do_binaryRepart();

    void setCharge() {
        this->getParticleContainer()->Q = qi_m;
    }
    
    void setMass() {
        this->getParticleContainer()->M = mi_m;
    }

    double getCharge() const {
        return qi_m*this->getTotalNum();
    }

    double getChargePerParticle() const {
        return qi_m;
    }
    double getMassPerParticle() const {
        return mi_m;
    }

    double getQ() const {
        return this->getCharge();
    }
    double getM() const {
        return  mi_m*this->getTotalNum();
    }

    double getdE() const;

    double getGamma(int /*i*/) const {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }
    double getBeta(int /*i*/) const {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }

    void actT() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }
  
    const PartData* getReference() const {
        return reference_m;
    }


  /// \todo constructor could set this
    void setReference (const PartData* ref) {
        reference_m = ref;
    }
  
    double getEmissionDeltaT() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 1.0;
    }

    void gatherLoadBalanceStatistics();

    size_t getLoadBalance(int p) {
        return globalPartPerNode_m[p];
    }

    void resizeMesh() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }

    bool isGridFixed() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return false;
    }

    void boundp() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }

    size_t boundp_destroyT() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 1;
    }

    void setBCAllOpen() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }
    void setBCForDCBeam() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }
    void setupBCs() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }
    void setBCAllPeriodic() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }

    void resetInterpolationCache(bool /*clearCache = false*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }
    void swap(unsigned int /*i*/, unsigned int /*j*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }
    double getRho(int /*x*/, int /*y*/, int /*z*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }
    void gatherStatistics(unsigned int /*totalP*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
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
        if (isUnitless_m) {
            throw OpalException("PartBunch::switchToUnitlessPositions",
                                "PartBunch is already in unitless positions!");
        }

        // Divide by c*dt
        double unitless_factor = 1.0 / (Physics::c * this->getdT());
        auto Rview  = this->getParticleContainer()->R.getView();
        auto dtview = this->getParticleContainer()->dt.getView();
        Kokkos::parallel_for(
                             "switchToUnitlessPositions", ippl::getRangePolicy(Rview),
                             KOKKOS_LAMBDA(const size_t i) {
                                double fac = use_dt_per_particle ? (1.0 / (Physics::c * dtview(i))) 
                                                                 : unitless_factor;
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
        Kokkos::parallel_for(
                             "switchOffUnitlessPositions", ippl::getRangePolicy(Rview),
                             KOKKOS_LAMBDA(const size_t i) {
                                double fac = use_dt_per_particle ? (Physics::c * dtview(i)) 
                                                                 : unitless_factor;
                                Rview(i) *= fac;
                             });
        isUnitless_m = false;

        /// \todo remove later
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

    void computeSelfFields();
    void dumpBinConfig(bool preMerge);

    Inform& print(Inform& os);

    bool hasFieldSolver() const {
        return this->fsolver_m != nullptr;
    }

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

    bool getFieldSolverType() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return false;
    }

    bool getIfBeamEmitting() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return false;
    }
    int getLastEmittedEnergyBin() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0;
    }
    size_t getNumberOfEmissionSteps() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0;
    }
    int getNumberOfEnergyBins() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0;
    }

    void Rebin() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }

    void setEnergyBins(int /*numberOfEnergyBins*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }
    bool weHaveEnergyBins() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return false;
    }
    void setTEmission(double /*t*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }
    double getTEmission() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }
    bool weHaveBins() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return false;
    }
    // void setPBins(PartBins* pbin) {}
    size_t emitParticles(double /*eZ*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0;
    }
    void updateNumTotal() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }
    void rebin() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }
    int getLastemittedBin() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0;
    }
    void setLocalBinCount(size_t /*num*/, int /*bin*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }
    void calcGammas() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }
    double getBinGamma(int /*bin*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }
    bool hasBinning() const {
        return this->bins_m != nullptr;
    }

    /** Number of field solver bins (1 if no binning, else getBins()->getCurrentBinCount()). */
    int getCurrentNBins() const {
        return hasBinning() ? static_cast<int>(bins_m->getCurrentBinCount()) : 1;
    }
    double calcMeanPhi() {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }
    double getPx(int /*i*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }
    double getPy(int /*i*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }
    double getPz(int /*i*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }
    double getPx0(int /*i*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }
    double getPy0(int /*i*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }
    double getX(int /*i*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }
    double getY(int /*i*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }
    double getZ(int /*i*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }
    double getX0(int /*i*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }
    double getY0(int /*i*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }

    void setZ(int /*i*/, double /*zcoo*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }

    void get_bounds(Vector_t<double, Dim>& rmin, Vector_t<double, Dim>& rmax) {
        rmin = rmin_m;
        rmax = rmax_m;
    }

    void getLocalBounds(Vector_t<double, Dim>& /*rmin*/, Vector_t<double, Dim>& /*rmax*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
    }

    void get_PBounds(Vector_t<double, Dim>& /*min*/, Vector_t<double, Dim>& /*max*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
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
    Vector_t<double, Dim> get_68Percentile() const {
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
    }
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

    /*
      Some quantities related to integrations/tracking
     */

    void setStepsPerTurn(int n) {
        stepsPerTurn_m = n;
    }

    int getStepsPerTurn() const {
        return stepsPerTurn_m;
    }

    /// step in multiple TRACK commands
    void setGlobalTrackStep(long long n) {
        globalTrackStep_m = n;
    }

    long long getGlobalTrackStep() const {
        return globalTrackStep_m;
    }

    /// step in a TRACK command
    void setLocalTrackStep(long long n) {
        localTrackStep_m = n;
    }

    void incTrackSteps() {
        globalTrackStep_m++;
        localTrackStep_m++;
    }

    long long getLocalTrackStep() const {
        return localTrackStep_m;
    }

    void setNumBunch(short n) {
        numBunch_m = n;
        bunchTotalNum_m.resize(n);
        bunchLocalNum_m.resize(n);
    }

    short getNumBunch() const {
        return numBunch_m;
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

    void setSteptoLastInj(int n) {
        SteptoLastInj_m = n;
    }

    int getSteptoLastInj() const {
        return SteptoLastInj_m;
    }

    double calculateAngle(double /*x*/, double /*y*/) {
        *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
        return 0.0;
    }

    // Sanity check functions
    void spaceChargeEFieldCheck(Vector_t<double, 3> efScale);

};

// Explicit instantiations
extern template class PartBunch<double, 3>;

#endif
