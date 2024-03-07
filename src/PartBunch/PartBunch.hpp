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

updateExternalFields(): check if bunch has access to the fields of eternal elements. Maybee phase
out some elements and read in new elements


diagnostics(): calculate statistics and maybee write tp h5 and stat files

*/

#include <memory>

#include "Algorithms/CoordinateSystemTrafo.h"
#include "Attributes/Attributes.h"
#include "Manager/BaseManager.h"
#include "Manager/PicManager.h"
#include "PartBunch/FieldContainer.hpp"
#include "PartBunch/FieldSolver.hpp"
#include "PartBunch/LoadBalancer.hpp"
#include "PartBunch/ParticleContainer.hpp"
#include "Random/Distribution.h"
#include "Random/InverseTransformSampling.h"
#include "Random/NormalDistribution.h"
#include "Random/Randn.h"

#include "Distribution/Distribution.h"

// #include "Structure/FieldSolverCmd.h"

#include "Algorithms/PartData.h"

using view_type = typename ippl::detail::ViewType<ippl::Vector<double, 3>, 1>::view_type;

// define functions used in sampling particles
struct CustomDistributionFunctions {
    struct CDF {
        KOKKOS_INLINE_FUNCTION double operator()(
            double x, unsigned int d, const double* params_p) const {
            return x
                   + (params_p[d * 2 + 0] / params_p[d * 2 + 1])
                         * Kokkos::sin(params_p[d * 2 + 1] * x);
        }
    };

    struct PDF {
        KOKKOS_INLINE_FUNCTION double operator()(
            double x, unsigned int d, double const* params_p) const {
            return 1.0 + params_p[d * 2 + 0] * Kokkos::cos(params_p[d * 2 + 1] * x);
        }
    };

    struct Estimate {
        KOKKOS_INLINE_FUNCTION double operator()(
            double u, unsigned int d, double const* params_p) const {
            return u + params_p[d] * 0.;
        }
    };
};

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

    Vector_t<int, Dim> nr_m;

    size_type totalP_m;

    int nt_m;
    double lbt_m;
    double dt_m;
    int it_m;

    std::string integration_method_m;
    std::string solver_m;

    Vector_t<double, Dim> origin_m;
    Vector_t<double, Dim> rmin_m;
    Vector_t<double, Dim> rmax_m;

    /// mesh size [m]
    Vector_t<double, Dim> hr_m;

    /// total charge in the bunch [Cb]
    double totalQ_m;

    // Landau damping specific
    double Bext_m;
    double alpha_m;
    double DrInv_m;

    ippl::NDIndex<Dim> domain_m;
    std::array<bool, Dim> decomp_m;
    bool isFirstRepartition_m;

    /*
      Up to here it is like the opaltest
    */

    /**
      Reference particle structures
     */

    Vector_t<double, Dim> RefPartR_m;
    Vector_t<double, Dim> RefPartP_m;

private:
    CoordinateSystemTrafo toLabTrafo_m;

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

    /// step in a TRACK command
    long long localTrackStep_m;

    /// if multiple TRACK commands
    long long globalTrackStep_m;

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
    /// Mesh enlargement
    double dh_m;  /// relative enlargement of the mesh

    PartData* reference_m;

    double couplingConstant_m;
    double qi_m;
    double massPerParticle_m;

    // unit state of PartBunch
    // UnitState_t unit_state_m;
    // UnitState_t stateOfLastBoundP_m;

    /// holds the centroid of the beam
    double centroid_m[2 * Dim];

    /// 6x6 matrix of the moments of the beam
    matrix_t moments_m;

    /// holds the actual time of the integration
    double t_m;

    /// the position along design trajectory
    double spos_m;

    /*
       flags to tell if we are a DC-beam
     */
    bool dcBeam_m;
    double periodLength_m;

    Distribution* OPALdist_m;

public:
    PartBunch(
        double totalCharge, Vector_t<int, Dim> nr, size_t totalP, int nt, double lbt,
        std::string integration_method, Distribution* OPALdistribution)
        : ippl::PicManager<
            T, Dim, ParticleContainer<T, Dim>, FieldContainer<T, Dim>, LoadBalancer<T, Dim>>(),
          time_m(0.0),
          nr_m(nr),
          totalP_m(totalP),
          nt_m(nt),
          lbt_m(lbt),
          dt_m(0),
          it_m(0),
          integration_method_m(integration_method),
          solver_m(""),
          totalQ_m(totalCharge),
          isFirstRepartition_m(true),
          OPALdist_m(OPALdistribution) {
        Inform m("PartBunch() ");

        for (unsigned i = 0; i < Dim; i++) {
            this->domain_m[i] = ippl::Index(nr[i]);
        }

        /*
          via OPALdist_m we get all information about
          the distribution

          - read fromFile
          - emmit from cathod
          - inject
         */

        Vector_t<double, 3> sigmaR = OPALdist_m->getSigmaR();
        Vector_t<double, 3> sigmaP = OPALdist_m->getSigmaP();
        this->decomp_m.fill(true);

        // some fake setup to get a distribution going

        this->rmin_m = 0.0;
        this->rmax_m = sigmaR;

        Vector_t<double, Dim> length = this->rmax_m - this->rmin_m;
        this->hr_m                   = length / this->nr_m;

        this->Bext_m   = 5.0;
        this->origin_m = this->rmin_m;

        this->dt_m = 0.5 / this->nr_m[2];

        this->alpha_m = -0.5 * this->dt_m;
        this->DrInv_m = 1.0 / (1 + (std::pow((this->alpha_m * this->Bext_m), 2)));

        using ParticleContainer_t = ParticleContainer<T, Dim>;
        using FieldContainer_t    = FieldContainer<T, Dim>;

        bool isAllPeriodic = true;

        this->setFieldContainer(std::make_shared<FieldContainer_t>(
            this->hr_m, this->rmin_m, this->rmax_m, this->decomp_m, this->domain_m, this->origin_m,
            isAllPeriodic));

        this->setParticleContainer(std::make_shared<ParticleContainer_t>(
            this->fcontainer_m->getMesh(), this->fcontainer_m->getFL()));
    }

    ~PartBunch() {
        Inform m("PartBunch Destructor ");
        m << "Finished time step: " << this->it_m << " time: " << this->time_m << endl;
    }

    std::shared_ptr<ParticleContainer_t> getParticleContainer() {
        return this->pcontainer_m;
    }

    void setSolver(std::string solver) {
        Inform m("PartBunch pre_run ");

        if (this->solver_m != "")
            m << "Warning solver already initiated" << endl;
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
        initializeParticles();

        static IpplTimings::TimerRef DummySolveTimer = IpplTimings::getTimer("solveWarmup");
        IpplTimings::startTimer(DummySolveTimer);

        this->fcontainer_m->getRho() = 0.0;

        this->fsolver_m->runSolver();

        IpplTimings::stopTimer(DummySolveTimer);

        this->par2grid();

        static IpplTimings::TimerRef SolveTimer = IpplTimings::getTimer("solve");
        IpplTimings::startTimer(SolveTimer);

        this->fsolver_m->runSolver();

        IpplTimings::stopTimer(SolveTimer);

        this->grid2par();
    }

public:
    size_t getTotalNum() const {
        return 0;
    }

    size_t getLocalNum() const {
        return 0;
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

    void initializeParticles() {
        Inform m("Initialize Particles");

        auto* mesh                   = &this->fcontainer_m->getMesh();
        auto* FL                     = &this->fcontainer_m->getFL();
        Vector_t<double, Dim> length = this->rmax_m - this->rmin_m;

        Vector_t<double, Dim> mu, sd;
        for (unsigned d = 0; d < Dim; d++) {
            mu[d] = 0.5 * length[d] + this->origin_m[d];
        }
        sd[0] = 0.15 * length[0];
        sd[1] = 0.05 * length[1];
        sd[2] = 0.20 * length[2];

        using DistR_t =
            ippl::random::Distribution<double, Dim, 2 * Dim, CustomDistributionFunctions>;

        double* parR = new double[2 * Dim];
        parR[0]      = mu[0];
        parR[1]      = sd[0];
        parR[2]      = mu[1];
        parR[3]      = sd[1];
        parR[4]      = mu[2];
        parR[5]      = sd[2];
        DistR_t distR(parR);

        Vector_t<double, Dim> origin = this->origin_m;
        if ((this->lbt_m != 1.0) && (ippl::Comm->size() > 1)) {
            m << "Starting first repartition" << endl;
            const ippl::NDIndex<Dim>& lDom = FL->getLocalNDIndex();
            const int nghost               = this->fcontainer_m->getRho().getNghost();
            auto rhoview                   = this->fcontainer_m->getRho().getView();

            using index_array_type = typename ippl::RangePolicy<Dim>::index_array_type;

            const auto hr = this->hr_m;
            ippl::parallel_for(
                "Assign initial rho based on PDF",
                this->fcontainer_m->getRho().getFieldRangePolicy(),
                KOKKOS_LAMBDA(const index_array_type& args) {
                    // local to global index conversion
                    Vector_t<double, Dim> xvec = (args + lDom.first() - nghost + 0.5) * hr + origin;

                    // ippl::apply accesses the view at the given indices and obtains a
                    // reference; see src/Expression/IpplOperations.h
                    ippl::apply(rhoview, args) = distR.getFullPdf(xvec);
                });

            Kokkos::fence();

            this->loadbalancer_m->initializeORB(FL, mesh);
            this->loadbalancer_m->repartition(FL, mesh, this->isFirstRepartition_m);
            isFirstRepartition_m = false;
        }

        // Sample particle positions:
        ippl::detail::RegionLayout<double, Dim, Mesh_t<Dim>> rlayout;
        rlayout = ippl::detail::RegionLayout<double, Dim, Mesh_t<Dim>>(*FL, *mesh);

        int seed        = 42;
        using size_type = ippl::detail::size_type;
        Kokkos::Random_XorShift64_Pool<> rand_pool64((size_type)(seed + 100 * ippl::Comm->rank()));

        using samplingR_t = ippl::random::InverseTransformSampling<
            double, Dim, Kokkos::DefaultExecutionSpace, DistR_t>;
        samplingR_t samplingR(distR, this->rmax_m, this->rmin_m, rlayout, this->totalP_m);

        size_type nlocal = samplingR.getLocalSamplesNum();

        this->pcontainer_m->create(nlocal);

        view_type* R_m = &this->pcontainer_m->R.getView();
        samplingR.generate(*R_m, rand_pool64);

        view_type* P_m = &this->pcontainer_m->P.getView();

        double muP[Dim] = {0.0, 0.0, 0.0};
        double sdP[Dim] = {1.0, 1.0, 1.0};
        Kokkos::parallel_for(nlocal, ippl::random::randn<double, Dim>(*P_m, rand_pool64, muP, sdP));

        Kokkos::fence();
        ippl::Comm->barrier();

        this->pcontainer_m->Q = this->totalQ_m / this->totalP_m;
        m << "particles created and initial conditions assigned " << endl;
    }

    void advance() override {
        if (this->integration_method_m == "LeapFrog") {
            LeapFrogStep();
        }
    }
    void LeapFrogStep() {
        // LeapFrog time stepping https://en.wikipedia.org/wiki/Leapfrog_integration
        // Here, we assume a constant charge-to-mass ratio of -1 for
        // all the particles hence eliminating the need to store mass as
        // an attribute
        Inform m("LeapFrog");
        Vector_t<double, Dim> len = this->rmax_m - this->rmin_m;
        const double alpha        = this->alpha_m;
        const double Bext         = this->Bext_m;
        const double DrInv        = this->DrInv_m;
        const double V0           = 30 * len[2];

        const auto ori = this->origin_m;

        std::shared_ptr<ParticleContainer_t> pc = this->pcontainer_m;
        std::shared_ptr<FieldContainer_t> fc    = this->fcontainer_m;

        auto Rview = pc->R.getView();
        auto Pview = pc->P.getView();
        auto Eview = pc->E.getView();

        Kokkos::parallel_for(
            "Kick1", pc->getLocalNum(), KOKKOS_LAMBDA(const size_t j) {
                double Eext_x =
                    -(Rview(j)[0] - ori[0] - 0.5 * len[0]) * (V0 / (2 * Kokkos::pow(len[2], 2)));
                double Eext_y =
                    -(Rview(j)[1] - ori[1] - 0.5 * len[1]) * (V0 / (2 * Kokkos::pow(len[2], 2)));
                double Eext_z =
                    (Rview(j)[2] - ori[2] - 0.5 * len[2]) * (V0 / (Kokkos::pow(len[2], 2)));

                Eext_x += Eview(j)[0];
                Eext_y += Eview(j)[1];
                Eext_z += Eview(j)[2];

                Pview(j)[0] += alpha * (Eext_x + Pview(j)[1] * Bext);
                Pview(j)[1] += alpha * (Eext_y - Pview(j)[0] * Bext);
                Pview(j)[2] += alpha * Eext_z;
            });
        Kokkos::fence();
        ippl::Comm->barrier();

        // drift
        pc->R = pc->R + dt_m * pc->P;

        // Since the particles have moved spatially update them to correct processors
        pc->update();

        size_type totalP = this->totalP_m;
        int it           = this->it_m;

        if (this->loadbalancer_m->balance(totalP, it + 1)) {
            auto* mesh = &fc->getRho().get_mesh();
            auto* FL   = &fc->getFL();
            this->loadbalancer_m->repartition(FL, mesh, isFirstRepartition_m);
        }

        // scatter the charge onto the underlying grid
        this->par2grid();

        // Field solve
        this->fsolver_m->runSolver();

        // gather E field
        this->grid2par();

        auto R2view = pc->R.getView();
        auto P2view = pc->P.getView();
        auto E2view = pc->E.getView();

        Kokkos::parallel_for(
            "Kick2", pc->getLocalNum(), KOKKOS_LAMBDA(const size_t j) {
                double Eext_x =
                    -(R2view(j)[0] - ori[0] - 0.5 * len[0]) * (V0 / (2 * Kokkos::pow(len[2], 2)));
                double Eext_y =
                    -(R2view(j)[1] - ori[1] - 0.5 * len[1]) * (V0 / (2 * Kokkos::pow(len[2], 2)));
                double Eext_z =
                    (R2view(j)[2] - ori[2] - 0.5 * len[2]) * (V0 / (Kokkos::pow(len[2], 2)));

                Eext_x += E2view(j)[0];
                Eext_y += E2view(j)[1];
                Eext_z += E2view(j)[2];

                P2view(j)[0] = DrInv
                               * (P2view(j)[0]
                                  + alpha * (Eext_x + P2view(j)[1] * Bext + alpha * Bext * Eext_y));
                P2view(j)[1] = DrInv
                               * (P2view(j)[1]
                                  + alpha * (Eext_y - P2view(j)[0] * Bext - alpha * Bext * Eext_x));
                P2view(j)[2] += alpha * Eext_z;
            });
        Kokkos::fence();
        ippl::Comm->barrier();
    }

    void par2grid() override {
        scatterCIC();
    }

    void grid2par() override {
        gatherCIC();
    }

    void gatherCIC() {
        gather(this->pcontainer_m->E, this->fcontainer_m->getE(), this->pcontainer_m->R);
    }

    void scatterCIC() {
        Inform m("scatter ");
        this->fcontainer_m->getRho() = 0.0;

        ippl::ParticleAttrib<double>* q          = &this->pcontainer_m->Q;
        typename Base::particle_position_type* R = &this->pcontainer_m->R;
        Field_t<Dim>* rho                        = &this->fcontainer_m->getRho();
        double Q                                 = this->totalQ_m;
        Vector_t<double, Dim> rmin               = rmin_m;
        Vector_t<double, Dim> rmax               = rmax_m;
        Vector_t<double, Dim> hr                 = hr_m;

        scatter(*q, *rho, *R);
        double relError = std::fabs((Q - (*rho).sum()) / Q);

        m << relError << endl;

        size_type TotalParticles = 0;
        size_type localParticles = this->pcontainer_m->getLocalNum();

        ippl::Comm->reduce(localParticles, TotalParticles, 1, std::plus<size_type>());

        if (ippl::Comm->rank() == 0) {
            if (TotalParticles != totalP_m || relError > 1e-10) {
                m << "Time step: " << it_m << endl;
                m << "Total particles in the sim. " << totalP_m << " "
                  << "after update: " << TotalParticles << endl;
                m << "Rel. error in charge conservation: " << relError << endl;
                ippl::Comm->abort();
            }
        }

        double cellVolume = std::reduce(hr.begin(), hr.end(), 1., std::multiplies<double>());
        (*rho)            = (*rho) / cellVolume;

        //        double rhoNorm = norm(*rho);
        // rho = rho_e - rho_i (only if periodic BCs)
        if (this->fsolver_m->getStype() != "OPEN") {
            double size = 1;
            for (unsigned d = 0; d < Dim; d++) {
                size *= rmax[d] - rmin[d];
            }
            *rho = *rho - (Q / size);
        }
    }

    Inform& print(Inform& os) {
        // if (this->getLocalNum() != 0) {  // to suppress Nans
        Inform::FmtFlags_t ff = os.flags();

        os << std::scientific;
        os << level1 << "\n";
        os << "* ************** B U N C H "
              "********************************************************* \n";
        os << "* PARTICLES       = " << this->pcontainer_m->getLocalNum() << "\n";
        os << "* CORES           = " << ippl::Comm->size() << "\n";
        os << "* FIELD SOLVER    = " << solver_m << "\n";
        os << "* INTEGRATOR      = " << integration_method_m << "\n";
        os << "* MESH SPACIN     = " << this->fcontainer_m->getMesh().getMeshSpacing() << "\n";
        os << "* FIELD LAYOUT    = " << this->fcontainer_m->getFL() << "\n";
        os << "* "
              "********************************************************************************"
              "** "
           << endl;
        os.flags(ff);
        // }
        return os;
    }

    /*
      Up to here it is like the opaltest
    */

    double getCouplingConstant() const {
        return 1.0;
    }
    void setCouplingConstant(double c) {
    }
    void calcBeamParameters() {
    }

    void setCharge(double q) {
    }
    void setMass(double mass) {
    }

    double getCharge() const {
        return 1.0;
    }
    double getChargePerParticle() const {
        return 1.0;
    }
    double getMassPerParticle() const {
        return 1.0;
    }

    double getQ() const {
        return 1.0;
    }
    double getM() const {
        return 1.0;
    }

    void setPType(const std::string& type) {
    }

    double getdE() const {
        return 1.0;
    }

    double getGamma(int i) const {
    }
    double getBeta(int i) const {
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

    Mesh_t<Dim>& getMesh() {
    }

    FieldLayout_t<Dim>& getFieldLayout() {
        return nullptr;
    }

    bool isGridFixed() {
    }

    void boundp() {
    }

    size_t boundp_destroyT() {
        return 1;
    }

    void setTotalNum(size_t newTotalNum) {
    }

    void set_meshEnlargement(double dh) {
    }
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
    }

    void calcLineDensity(
        unsigned int nBins, std::vector<double>& lineDensity, std::pair<double, double>& meshInfo) {
    }

    void setBeamFrequency(double v) {
    }

    Vector_t<double, Dim> getEExtrema() {
    }

    void computeSelfFields() {
    }
    void computeSelfFields(int b) {
    }

    bool hasFieldSolver() {
        return true;
    }

    bool getFieldSolverType() {
    }

    bool getIfBeamEmitting() {
    }
    int getLastEmittedEnergyBin() {
    }
    size_t getNumberOfEmissionSteps() {
    }
    int getNumberOfEnergyBins() {
    }

    void Rebin() {
    }

    void setEnergyBins(int numberOfEnergyBins) {
    }
    bool weHaveEnergyBins() {
    }
    void setTEmission(double t) {
    }
    double getTEmission() {
    }
    bool weHaveBins() {
    }
    // void setPBins(PartBins* pbin) {}
    size_t emitParticles(double eZ) {
    }
    void updateNumTotal() {
    }
    void rebin() {
    }
    int getLastemittedBin() {
    }
    void setLocalBinCount(size_t num, int bin) {
    }
    void calcGammas() {
    }
    double getBinGamma(int bin) {
    }
    bool hasBinning() {
    }
    void setBinCharge(int bin, double q) {
    }
    void setBinCharge(int bin) {
    }
    double calcMeanPhi() {
    }
    bool resetPartBinID2(const double eta) {
    }
    bool resetPartBinBunch() {
    }
    double getPx(int i) {
    }
    double getPy(int i) {
    }
    double getPz(int i) {
    }
    double getPx0(int i) {
    }
    double getPy0(int i) {
    }
    double getX(int i) {
    }
    double getY(int i) {
    }
    double getZ(int i) {
    }
    double getX0(int i) {
    }
    double getY0(int i) {
    }

    void setZ(int i, double zcoo) {
    }

    void get_bounds(Vector_t<double, Dim>& rmin, Vector_t<double, Dim>& rmax) {
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

    Vector_t<double, Dim> get_centroid() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_rrms() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_rprms() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_rmean() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_prms() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_pmean() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_pmean_Distribution() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_emit() const {
        return Vector_t<double, Dim>(0.0);
    }
    Vector_t<double, Dim> get_norm_emit() const {
        return Vector_t<double, Dim>(0.0);
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
        return 0.0;
    }
    double get_Dy() const {
        return 0.0;
    }
    double get_DDx() const {
        return 0.0;
    }
    double get_DDy() const {
        return 0.0;
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
