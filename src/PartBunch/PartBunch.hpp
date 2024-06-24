#ifndef PARTBUNCH_H
#define PARTBUNCH_H

/*
  Notes:


Additional Functions
--------------------

   setDistribution
   getDistribution


Main loop
---------


   for integration
     if not injectionDome
       injected = distr_m->create()
       injectrd = distr_m->inject()
       update()

     drift()
     update()
     kick()
     drift()

     diagnostics()
     updateExternalFields()
     update()

updateExternalFields(): check if bunch has access to the fields of eternal elements. Maybe phase
out some elements and read in new elements


diagnostics(): calculate statistics and maybe write tp h5 and stat files

*/

#include <memory>

#include "Algorithms/BoostMatrix.h"
#include "Algorithms/CoordinateSystemTrafo.h"
//#include "Algorithms/DistributionMoments.h"
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

#include "Structure/FieldSolverCmd.h"

#include "Algorithms/PartData.h"

using view_type = typename ippl::detail::ViewType<ippl::Vector<double, 3>, 1>::view_type;


template <typename T, unsigned Dim>
class PartBunch
    : public ippl::PicManager<
          T, Dim, ParticleContainer<T, Dim>, FieldContainer<T, Dim>, LoadBalancer<T, Dim>> {
public:
    using ParticleContainer_t = ParticleContainer<T, Dim>;
    using FieldContainer_t    = FieldContainer<T, Dim>;
    using FieldSolver_t       = FieldSolver<T, Dim>;
    using LoadBalancer_t      = LoadBalancer<T, Dim>;
    using Base                = ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>;

    double time_m;

    size_type totalP_m;

    int nt_m;

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

    /// step in a TRACK command
    long long localTrackStep_m;

    /// if multiple TRACK commands
    long long globalTrackStep_m;


    std::shared_ptr<Distribution> OPALdist_m;

    std::shared_ptr<FieldSolverCmd> OPALFieldSolver_m;

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
    std::unique_ptr<size_t[]> binemitted_m;

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

    PartData* reference_m;

    double couplingConstant_m;

    // unit state of PartBunch
    // UnitState_t unit_state_m;
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

public:
    PartBunch(
              double qi, double mi, size_t totalP, int nt, double lbt, std::string integration_method,
        std::shared_ptr<Distribution> &OPALdistribution,
        std::shared_ptr<FieldSolverCmd> &OPALFieldSolver)
        : ippl::PicManager<
            T, Dim, ParticleContainer<T, Dim>, FieldContainer<T, Dim>, LoadBalancer<T, Dim>>(),
          time_m(0.0),
          totalP_m(totalP),
          nt_m(nt),
          lbt_m(lbt),
          dt_m(0),
          it_m(0),
          integration_method_m(integration_method),
          solver_m(""),
          isFirstRepartition_m(true),
          qi_m(qi),
          mi_m(mi),
          rmsDensity_m(0.0),
          localTrackStep_m(0),
          globalTrackStep_m(0),
          OPALdist_m(OPALdistribution),
          OPALFieldSolver_m(OPALFieldSolver){

        Inform m("PartBunch() ");

        static IpplTimings::TimerRef gatherInfoPartBunch = IpplTimings::getTimer("gatherInfoPartBunch");
        IpplTimings::startTimer(gatherInfoPartBunch);

        /*
          get the needed information from OPAL FieldSolver command
        */

        nr_m = Vector_t<int, Dim>(
            OPALFieldSolver_m->getNX(), OPALFieldSolver_m->getNY(), OPALFieldSolver_m->getNZ());

        const Vector_t<bool, 3> domainDecomposition = OPALFieldSolver_m->getDomDec();

        for (unsigned i = 0; i < Dim; i++) {
            this->domain_m[i] = ippl::Index(nr_m[i]);
            this->decomp_m[i] = domainDecomposition[i];
        }

        bool isAllPeriodic = true;  // \fixme need to get BCs from OPAL Fieldsolver

        /*
          set stuff for pre_run i.e. warmup
          this will be reset when the correct computational
          domain is set
        */

        Vector_t<double, Dim> length (6.0);
        this->hr_m = length / this->nr_m;
        this->origin_m = -3.0;
        this->dt_m = 0.5 / this->nr_m[2];

        rmin_m = origin_m;
        rmax_m = origin_m + length;

        this->setFieldContainer( std::make_shared<FieldContainer_t>(hr_m, rmin_m, rmax_m, decomp_m, domain_m, origin_m, isAllPeriodic) );

        this->setParticleContainer(std::make_shared<ParticleContainer_t>(
            this->fcontainer_m->getMesh(), this->fcontainer_m->getFL()));

        IpplTimings::stopTimer(gatherInfoPartBunch);

        static IpplTimings::TimerRef setSolverT = IpplTimings::getTimer("setSolver");
        IpplTimings::startTimer(setSolverT);
        setSolver(OPALFieldSolver_m->getType());
        IpplTimings::stopTimer(setSolverT);

        static IpplTimings::TimerRef prerun = IpplTimings::getTimer("prerun");
        IpplTimings::startTimer(prerun);
        pre_run();
        IpplTimings::stopTimer(prerun);

        /*
          end wrmup 
        */

    }

    void bunchUpdate();

    ~PartBunch() {
        Inform m("PartBunch Destructor ");
        m << "Finished time step: " << this->it_m << " time: " << this->time_m << endl;
    }

    std::shared_ptr<ParticleContainer_t> getParticleContainer() {
        return this->pcontainer_m;
    }

    void setSolver(std::string solver) {
        Inform m("setSolver  ");

        if (this->solver_m != "")
            m << "Warning solver already initiated but overwrite ..." << endl;

        this->solver_m = solver;

        this->fcontainer_m->initializeFields(this->solver_m);

        this->setFieldSolver(std::make_shared<FieldSolver_t>(
            this->solver_m, &this->fcontainer_m->getRho(), &this->fcontainer_m->getE(),
            &this->fcontainer_m->getPhi()));

        this->fsolver_m->initSolver();

        /// ADA we need to be able to set a load balancer when not having a field solver
        this->setLoadBalancer(std::make_shared<LoadBalancer_t>(
            this->lbt_m, this->fcontainer_m, this->pcontainer_m, this->fsolver_m));
    }

    void pre_run() override {
        static IpplTimings::TimerRef DummySolveTimer = IpplTimings::getTimer("solveWarmup");
        IpplTimings::startTimer(DummySolveTimer);

        this->fcontainer_m->getRho() = 0.0;

        this->fsolver_m->runSolver();

        IpplTimings::stopTimer(DummySolveTimer);

    }

public:
    void updateMoments(){
        this->pcontainer_m->updateMoments();
    }

    size_t getTotalNum() const {
        return this->pcontainer_m->getTotalNum();
    }

    size_t getLocalNum() const {
        return this->pcontainer_m->getLocalNum();
    }

    Vector_t<double, Dim> R(size_t i) {
        return Vector_t<double, Dim>(0.0);
    }

    Vector_t<double, Dim> P(size_t i) {
        return Vector_t<double, Dim>(0.0);
    }

    Vector_t<double, Dim> Ef(size_t i) {
        return Vector_t<double, Dim>(0.0);
    }

    Vector_t<double, Dim> Bf(size_t i) {
        return Vector_t<double, Dim>(0.0);
    }

    Vector_t<double, Dim> dt(size_t i) {
        return Vector_t<double, Dim>(0.0);
    }

    void advance() override {
        // \todo needs to go
    }

    void par2grid() override {
        scatterCIC();
    }

    void grid2par() override {
        gatherCIC();
    }

    void gatherCIC() {
/*
        using Base = ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>;
        typename Base::particle_position_type* Ep = &this->pcontainer_m->E;
        typename Base::particle_position_type* R = &this->pcontainer_m->R;
        VField_t<T, Dim>* Ef = &this->fcontainer_m->getE();
        gather(*Ep, *Ef, *R);
*/
    }

    void scatterCIC(); 

    /*
      Up to here it is like the opaltest
    */

    double getCouplingConstant() const {
        return 1.0;
    }

    void setCouplingConstant(double c) {
    }

    void calcBeamParameters();

    void do_binaryRepart()
    {
        using FieldContainer_t    = FieldContainer<T, Dim>;
        std::shared_ptr<FieldContainer_t> fc    = this->fcontainer_m; 

        size_type totalP = this->getTotalNum() ;
        int it           = this->it_m;
        if (this->loadbalancer_m->balance(totalP, it + 1)) {
            auto* mesh = &fc->getRho().get_mesh();
            auto* FL   = &fc->getFL();
            this->loadbalancer_m->repartition(FL, mesh, isFirstRepartition_m);
        }
    }

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

    double getdE() const {
        return this->pcontainer_m->getStdKineticEnergy();
    }

    double getGamma(int i) const {
        return 0.0;
    }
    double getBeta(int i) const {
        return 0.0;
    }

    void actT() {
    }

    PartData* getReference() {
        return reference_m;
    }

    double getEmissionDeltaT() {
        return 1.0;
    }

    void gatherLoadBalanceStatistics() {
    }

    size_t getLoadBalance(int p) {
        return 0;
    }

    void resizeMesh() {
    }

    /*

    Mesh_t<Dim>& getMesh() {
    }

    FieldLayout_t<Dim>& getFieldLayout() {
        return nullptr;
    }
    */



    bool isGridFixed() {
        return false;
    }

    void boundp() {
    }

    size_t boundp_destroyT() {
        return 1;
    }

    /*
    void setTotalNum(size_t newTotalNum) {
    }

    void set_meshEnlargement(double dh) {
    }
    */

    void setBCAllOpen() {
    }
    void setBCForDCBeam() {
    }
    void setupBCs() {
    }
    void setBCAllPeriodic() {
    }

    void resetInterpolationCache(bool clearCache = false) {
    }
    void swap(unsigned int i, unsigned int j) {
    }
    double getRho(int x, int y, int z) {
    }
    void gatherStatistics(unsigned int totalP) {
    }
    void switchToUnitlessPositions(bool use_dt_per_particle = false) {
    }
    void switchOffUnitlessPositions(bool use_dt_per_particle = false) {
    }

    size_t calcNumPartsOutside(Vector_t<double, Dim> x) {
        return 0;
    }

    void calcLineDensity(
        unsigned int nBins, std::vector<double>& lineDensity, std::pair<double, double>& meshInfo) {
    }

    void setBeamFrequency(double v) {
    }

    Vector_t<double, Dim> getEExtrema() {
       return Vector_t<double, Dim>(0);
    }

    void computeSelfFields();

    Inform& print(Inform& os);


    bool hasFieldSolver() {
        return true;
    }

    bool getFieldSolverType() {
        return false;
    }

    bool getIfBeamEmitting() {
        return false;
    }
    int getLastEmittedEnergyBin() {
        return 0;
    }
    size_t getNumberOfEmissionSteps() {
        return 0;
    }
    int getNumberOfEnergyBins() {
        return 0;
    }

    void Rebin() {
    }

    void setEnergyBins(int numberOfEnergyBins) {
    }
    bool weHaveEnergyBins() {
        return false;
    }
    void setTEmission(double t) {
    }
    double getTEmission() {
         return 0.0;
    }
    bool weHaveBins() {
          return false;
    }
    // void setPBins(PartBins* pbin) {}
    size_t emitParticles(double eZ) {
           return 0;
    }
    void updateNumTotal() {
    }
    void rebin() {
    }
    int getLastemittedBin() {
        return 0;
    }
    void setLocalBinCount(size_t num, int bin) {
    }
    void calcGammas() {
    }
    double getBinGamma(int bin) {
        return 0.0;
    }
    bool hasBinning() {
        return false;
    }
    void setBinCharge(int bin, double q) {
    }
    void setBinCharge(int bin) {
    }
    double calcMeanPhi() {
        return 0.0;
    }
    bool resetPartBinID2(const double eta) {
        return false;
    }
    bool resetPartBinBunch() {
        return false;
    }
    double getPx(int i) {
        return 0.0;
    }
    double getPy(int i) {
        return 0.0;
    }
    double getPz(int i) {
        return 0.0;
    }
    double getPx0(int i) {
        return 0.0;
    }
    double getPy0(int i) {
        return 0.0;
    }
    double getX(int i) {
        return 0.0;
    }
    double getY(int i) {
        return 0.0;
    }
    double getZ(int i) {
        return 0.0;
    }
    double getX0(int i) {
        return 0.0;
    }
    double getY0(int i) {
        return 0.0;
    }

    void setZ(int i, double zcoo) {
    }

    void get_bounds(Vector_t<double, Dim>& rmin, Vector_t<double, Dim>& rmax) {
        rmin = rmin_m;
        rmax = rmax_m;
    }

    void getLocalBounds(Vector_t<double, Dim>& rmin, Vector_t<double, Dim>& rmax) {
    }

    void get_PBounds(Vector_t<double, Dim>& min, Vector_t<double, Dim>& max) {
    }

    /*

      Misc Bunch related quantities

    */

    /// get 2nd order beam matrix
    // matrix_t getSigmaMatrix() const;

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

    /**
     * get the spos of the primary beam
     *
     * @param none
     *
     */

    void set_sPos(double s) {
        spos_m = s;
    }

    double get_sPos() const {
        return spos_m;
    }

    double get_gamma() const {
        return 1.00;
    }

    double get_meanKineticEnergy() {
        return 0.0;
    }

    Vector_t<double, Dim> get_origin() const {
        return rmin_m;
    }
    Vector_t<double, Dim> get_maxExtent() const {
        return rmax_m;
    }

    // in opal, MeanPosition is return for get_centroid, which I think is wrong. We already have get_rmean()
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
    Vector_t<double, Dim> get_pmean_Distribution() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_emit() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_norm_emit() const {
        return this->pcontainer_m->getNormEmit();
    }
    Vector_t<double, Dim> get_halo() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_68Percentile() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_95Percentile() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_99Percentile() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_99_99Percentile() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_normalizedEps_68Percentile() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_normalizedEps_95Percentile() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_normalizedEps_99Percentile() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_normalizedEps_99_99Percentile() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_hr() const {
        return Vector_t<double, Dim>(0.0);
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

    double calculateAngle(double x, double y) {
        return 0.0;
    }
};

/* \todo
Inform& operator<<(Inform& os, PartBunch<double, 3>& p) {
    return p.print(os);
}
*/

#endif
