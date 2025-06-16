#ifndef IPPL_LANDAU_DAMPING_MANAGER_H
#define IPPL_LANDAU_DAMPING_MANAGER_H

#include <memory>

#include "FieldContainer.hpp"
#include "FieldSolver.hpp"
#include "LoadBalancer.hpp"
#include "AlpineManager.h"
#include "Manager/BaseManager.h"
#include "ParticleContainer.hpp"
#include "Random/Distribution.h"
#include "Random/InverseTransformSampling.h"
#include "Random/NormalDistribution.h"
#include "Random/Randn.h"

#include "AdaptBins.h" // TODO: Binning

using view_type = typename ippl::detail::ViewType<ippl::Vector<double, Dim>, 1>::view_type;

// define functions used in sampling particles
struct CustomDistributionFunctions {
  struct CDF{
       KOKKOS_INLINE_FUNCTION double operator()(double x, unsigned int d, const double *params_p) const {
           return x + (params_p[d * 2 + 0] / params_p[d * 2 + 1]) * Kokkos::sin(params_p[d * 2 + 1] * x);
       }
  };

  struct PDF{
       KOKKOS_INLINE_FUNCTION double operator()(double x, unsigned int d, double const *params_p) const {
           return 1.0 + params_p[d * 2 + 0] * Kokkos::cos(params_p[d * 2 + 1] * x);
       }
  };

  struct Estimate{
        KOKKOS_INLINE_FUNCTION double operator()(double u, unsigned int d, double const *params_p) const {
            return u + params_p[d] * 0.;
	}
  };
};

template <typename T, unsigned Dim>
class LandauDampingManager : public AlpineManager<T, Dim> {
public:
    using ParticleContainer_t = ParticleContainer<T, Dim>;
    using FieldContainer_t = FieldContainer<T, Dim>;
    using FieldSolver_t= FieldSolver<T, Dim>;
    using LoadBalancer_t= LoadBalancer<T, Dim>;

    // TODO: Binning
    using BinningSelector_t = typename ParticleBinning::CoordinateSelector<ParticleContainer_t>;
    using AdaptBins_t       = typename ParticleBinning::AdaptBins<ParticleContainer_t, BinningSelector_t>;

    VField_t<double, 3> E_tmp; // temporary field for adding up the lorentz transformed E field
    std::string samplingMethod_m;

    LandauDampingManager(size_type totalP_, int nt_, Vector_t<int, Dim> &nr_,
                       double lbt_, std::string& solver_, std::string& stepMethod_, std::string& samplingMethod_)
        : AlpineManager<T, Dim>(totalP_, nt_, nr_, lbt_, solver_, stepMethod_),
          samplingMethod_m(samplingMethod_) {}

    ~LandauDampingManager(){}

    void pre_run() override {
        Inform m("Pre Run");

        if (this->solver_m == "OPEN") {
            throw IpplException("LandauDamping", "Open boundaries solver incompatible with this simulation!");
        }

        for (unsigned i = 0; i < Dim; i++) {
            this->domain_m[i] = ippl::Index(this->nr_m[i]);
        }

        this->decomp_m.fill(true);
        this->kw_m    = 0.5;
        this->alpha_m = 0.05;
        this->rmin_m  = 0.0;
        this->rmax_m  = 2 * pi / this->kw_m;

        this->hr_m = this->rmax_m / this->nr_m;
        // Q = -\int\int f dx dv
        this->Q_m = std::reduce(this->rmax_m.begin(), this->rmax_m.end(), -1., std::multiplies<double>());
        this->origin_m = this->rmin_m;
        this->dt_m     = std::min(.05, 0.5 * *std::min_element(this->hr_m.begin(), this->hr_m.end()));
        this->it_m     = 0;
        this->time_m   = 0.0;

        m << "Discretization:" << endl
          << "nt " << this->nt_m << " Np= " << this->totalP_m << " grid = " << this->nr_m << endl;

        this->isAllPeriodic_m = true;

        this->setFieldContainer( std::make_shared<FieldContainer_t>( this->hr_m, this->rmin_m, this->rmax_m, this->decomp_m, this->domain_m, this->origin_m, this->isAllPeriodic_m) );

        this->setParticleContainer( std::make_shared<ParticleContainer_t>( this->fcontainer_m->getMesh(), this->fcontainer_m->getFL()) );

        this->fcontainer_m->initializeFields(this->solver_m);

        this->setFieldSolver( std::make_shared<FieldSolver_t>( this->solver_m, &this->fcontainer_m->getRho(), &this->fcontainer_m->getE(), &this->fcontainer_m->getPhi()) );

        this->fsolver_m->initSolver();

        this->setLoadBalancer( std::make_shared<LoadBalancer_t>( this->lbt_m, this->fcontainer_m, this->pcontainer_m, this->fsolver_m) );

        if (samplingMethod_m == "Landau") {
            initializeParticles();
        } else if (samplingMethod_m == "Flattop") {
            initializeParticlesFlattop();
        } else {
            throw IpplException("LandauDamping", "Unknown sampling method!");
        }

        //TODO: Binning - Create the bins object
        this->setBins(std::make_shared<AdaptBins_t>(
            this->getParticleContainer(), 
            BinningSelector_t(2), // no need to be a pointer, is only used inside the AdaptBins class
            256,
            1.0, 1.5, 0.05 // cost function parameters
        ));
        this->bins_m->debug();

        // TODO: Binning - After initializing the particles, create the limits
        //this->bins_m->initLimits();
        //this->bins_m->doFullRebin(10, true, ParticleBinning::HistoReductionMode::TeamBased); // test with 10 bins
        this->bins_m->doFullRebin(10);
        this->bins_m->print(); // TODO For debugging...
        this->E_tmp.initialize(this->fcontainer_m->getMesh(), this->fcontainer_m->getFL()); // initialize temporary field 


        static IpplTimings::TimerRef DummySolveTimer  = IpplTimings::getTimer("solveWarmup");
        IpplTimings::startTimer(DummySolveTimer);

        this->fcontainer_m->getRho() = 0.0;

        this->fsolver_m->runSolver();

        IpplTimings::stopTimer(DummySolveTimer);

        if (this->pcontainer_m->getLocalNum() > 0) this->par2grid();

        static IpplTimings::TimerRef SolveTimer = IpplTimings::getTimer("solve");
        IpplTimings::startTimer(SolveTimer);

        this->fsolver_m->runSolver();

        IpplTimings::stopTimer(SolveTimer);

        this->grid2par();

        this->dump();

        m << "Done";
    }

    void initializeParticles(){
        Inform m("Initialize Particles");

        auto *mesh = &this->fcontainer_m->getMesh();
        auto *FL = &this->fcontainer_m->getFL();
        using DistR_t = ippl::random::Distribution<double, Dim, 2 * Dim, CustomDistributionFunctions>;
        double parR[2 * Dim];
        for(unsigned int i=0; i<Dim; i++){
            parR[i * 2   ]  = this->alpha_m;
            parR[i * 2 + 1] = this->kw_m[i];
        }
        DistR_t distR(parR);

        Vector_t<double, Dim> kw     = this->kw_m;
        Vector_t<double, Dim> hr     = this->hr_m;
        Vector_t<double, Dim> origin = this->origin_m;
        static IpplTimings::TimerRef domainDecomposition = IpplTimings::getTimer("loadBalance");
        if ((this->lbt_m != 1.0) && (ippl::Comm->size() > 1)) {
            m << "Starting first repartition" << endl;
            IpplTimings::startTimer(domainDecomposition);
            this->isFirstRepartition_m           = true;
            const ippl::NDIndex<Dim>& lDom = FL->getLocalNDIndex();
            const int nghost               = this->fcontainer_m->getRho().getNghost();
            auto rhoview                   = this->fcontainer_m->getRho().getView();

            using index_array_type = typename ippl::RangePolicy<Dim>::index_array_type;
            ippl::parallel_for(
                "Assign initial rho based on PDF", this->fcontainer_m->getRho().getFieldRangePolicy(),
                KOKKOS_LAMBDA (const index_array_type& args) {
                    // local to global index conversion
                    Vector_t<double, Dim> xvec =
                        (args + lDom.first() - nghost + 0.5) * hr + origin;

                    // ippl::apply accesses the view at the given indices and obtains a
                    // reference; see src/Expression/IpplOperations.h
                    ippl::apply(rhoview, args) = distR.getFullPdf(xvec);
                });

            Kokkos::fence();

            this->loadbalancer_m->initializeORB(FL, mesh);
            this->loadbalancer_m->repartition(FL, mesh, this->isFirstRepartition_m);
            IpplTimings::stopTimer(domainDecomposition);
        }

        static IpplTimings::TimerRef particleCreation = IpplTimings::getTimer("particlesCreation");
        IpplTimings::startTimer(particleCreation);

        // Sample particle positions:
        ippl::detail::RegionLayout<double, Dim, Mesh_t<Dim>> rlayout;
        rlayout = ippl::detail::RegionLayout<double, Dim, Mesh_t<Dim>>( *FL, *mesh );

        // unsigned int
        size_type totalP = this->totalP_m;
        int seed           = 42;
        Kokkos::Random_XorShift64_Pool<> rand_pool64((size_type)(seed + 100 * ippl::Comm->rank()));

        using samplingR_t =
            ippl::random::InverseTransformSampling<double, Dim, Kokkos::DefaultExecutionSpace,
                                                   DistR_t>;
        Vector_t<double, Dim> rmin = this->rmin_m;
        Vector_t<double, Dim> rmax = this->rmax_m;
        samplingR_t samplingR(distR, rmax, rmin, rlayout, totalP);
        size_type nlocal = samplingR.getLocalSamplesNum();

        this->pcontainer_m->create(nlocal);

        view_type* R = &(this->pcontainer_m->R.getView());
        samplingR.generate(*R, rand_pool64);
        // this->pcontainer_m->R = (this->pcontainer_m->R * this->pcontainer_m->R) / 13; // change distribution a bit for binning tests

        view_type* P = &(this->pcontainer_m->P.getView());

        double mu[Dim];
        double sd[Dim];
        for(unsigned int i=0; i<Dim; i++) {
            mu[i] = 0.0;
            sd[i] = 1.0;
        }
        Kokkos::parallel_for(nlocal, ippl::random::randn<double, Dim>(*P, rand_pool64, mu, sd));
        Kokkos::fence();
        ippl::Comm->barrier();

        IpplTimings::stopTimer(particleCreation);

        this->pcontainer_m->q = this->Q_m/totalP;
        m << "particles created and initial conditions assigned " << endl;
    }

    void initializeParticlesFlattop() {
        // Reserve a few particles!
        double overalloc = ippl::Comm->getDefaultOverallocation();
        size_type totalP = this->totalP_m;
        size_type nlocal = (totalP / ippl::Comm->size() + 1) * overalloc;
        this->pcontainer_m->create(nlocal);

        Kokkos::View<bool*> tmp_invalid("tmp_invalid", 0);
        this->pcontainer_m->destroy(tmp_invalid, nlocal);
    }

    double flatTopProfile(double t, double riseTime, double sigmaTRise,
                          double flattopTime, double sigmaTFall, double fallTime) {
        if (t < riseTime) {
            // Rising tail; center is at riseTime.
            return std::exp( - std::pow((t - riseTime) / sigmaTRise, 2) / 2.0 );
        } else if (t >= riseTime && t < riseTime + flattopTime) {
            return 1.0;
        } else if (t >= riseTime + flattopTime && t < riseTime + flattopTime + fallTime) {
            // Falling tail; center is at (riseTime + flattopTime + fallTime)
            return std::exp( - std::pow((t - (riseTime + flattopTime + fallTime)) / sigmaTFall, 2) / 2.0 );
        } else {
            return 0.0;
        }
    }
    double integrateTrapezoidal(double x1, double x2, double y1, double y2) {
        return 0.5 * (y1 + y2) * std::fabs(x2 - x1);
    }
    size_type countEnteringParticles(double t0, double tf, double riseTime, double sigmaTRise,
                                     double flattopTime, double sigmaTFall, double fallTime, size_type totalN,
                                     double distArea) {
        double y1 = flatTopProfile(t0, riseTime, sigmaTRise, flattopTime, sigmaTFall, fallTime);
        double y2 = flatTopProfile(tf, riseTime, sigmaTRise, flattopTime, sigmaTFall, fallTime);
        double tArea = integrateTrapezoidal(t0, tf, y1, y2);

        std::size_t totalNew = static_cast<size_type>(std::floor(totalN * tArea / distArea));
        return totalNew;
    }
    size_type computeLocalParticleCount(size_type totalNew) {
        int rank = ippl::Comm->rank();
        int numRanks = ippl::Comm->size();
        std::size_t nlocal = totalNew / numRanks;
        std::size_t remainder = totalNew % numRanks;
        if (rank == 0) {
            nlocal += remainder;
        }
        return nlocal;
    }


    void emitFlattop(double t, double dt) {
        Inform m("Emit Flattop");

        size_type totalNew = countEnteringParticles(t, t + dt, 2.0, 1.0, 5.0, 1.0, 2.0, this->totalP_m, 1.0);
        size_type nNew = computeLocalParticleCount(totalNew);
        
        if (nNew <= 0) {
            m << "No new particles to emit at time " << t << " (dt=" << dt << ")." << endl;
            return;
        }
        
        m << "Emitting " << nNew << " new particles at time " << t << endl;
        
        size_type oldNum = this->pcontainer_m->getLocalNum();
        
        this->pcontainer_m->create(nNew);

        auto *mesh = &this->fcontainer_m->getMesh();
        auto *FL   = &this->fcontainer_m->getFL();
        
        ippl::detail::RegionLayout<double, Dim, Mesh_t<Dim>> rlayout(*FL, *mesh);
        
        using DistR_t = ippl::random::Distribution<double, Dim, 2 * Dim, CustomDistributionFunctions>;
        double parR[2 * Dim];
        for (unsigned int i = 0; i < Dim; i++) {
            parR[i * 2]     = this->alpha_m;
            parR[i * 2 + 1] = this->kw_m[i];
        }
        DistR_t distR(parR);

        using samplingR_t = ippl::random::InverseTransformSampling<double, Dim, Kokkos::DefaultExecutionSpace, DistR_t>;
        Vector_t<double, Dim> rmin = this->rmin_m;
        Vector_t<double, Dim> rmax = this->rmax_m;
        samplingR_t samplingR(distR, rmax, rmin, rlayout, nNew);

        // Get the view for positions. We only want to fill in the newly created indices.
        view_type fullR = this->pcontainer_m->R.getView();
        auto new_R = Kokkos::subview(fullR, std::make_pair(oldNum, oldNum + nNew));
        
        Kokkos::Random_XorShift64_Pool<> rand_pool64((size_type)(42 + 100 * ippl::Comm->rank() + this->getNt()));

        samplingR.generate(new_R, rand_pool64);

        view_type fullP = this->pcontainer_m->P.getView();
        auto new_P = Kokkos::subview(fullP, std::make_pair(oldNum, oldNum + nNew));
        
        double mu[Dim];
        double sd[Dim];
        for (unsigned int i = 0; i < Dim; i++) {
            mu[i] = 0.0;
            sd[i] = 1.0;
        }
        
        Kokkos::parallel_for("SampleVelocities", new_P.extent(0),
            ippl::random::randn<double, Dim>(new_P, rand_pool64, mu, sd));
        Kokkos::fence();

        this->pcontainer_m->q = this->Q_m/this->totalP_m;

        m << "Particles emitted and new positions/velocities sampled." << endl;
    }

    void advance() override {
        if (this->stepMethod_m == "LeapFrog") {
            LeapFrogStep();
            //this->bins_m->print(); // just some debug output TODO: binning, remove later
        } else {
            throw IpplException(TestName, "Step method is not set/recognized!");
        }
    }

    void LeapFrogStep(){
        // LeapFrog time stepping https://en.wikipedia.org/wiki/Leapfrog_integration
        // Here, we assume a constant charge-to-mass ratio of -1 for
        // all the particles hence eliminating the need to store mass as
        // an attribute

        // Create 5 new particles every timestep with random position values and 0 velocity
        /*{
            size_type nnew = 1;
            size_type nlocal = this->pcontainer_m->getLocalNum();
            this->pcontainer_m->create(nnew);
            view_type R = this->pcontainer_m->R.getView();
            view_type P = this->pcontainer_m->P.getView();
            Kokkos::deep_copy(Kokkos::subview(P, Kokkos::make_pair(nlocal, nlocal + nnew)), 0.0); // set new particles' velocity to 0

            Kokkos::Random_XorShift64_Pool<> rand_pool64((size_type)(42 + 100 * ippl::Comm->rank()));
            Kokkos::parallel_for(Kokkos::RangePolicy(nlocal, nlocal+nnew), ippl::random::randn<double, Dim>(R, rand_pool64));
            
            ippl::Comm->reduce(this->pcontainer_m->getLocalNum(), this->totalP_m, 1, std::plus<size_type>());
            this->pcontainer_m->q = this->Q_m/this->totalP_m;
        }*/

        static IpplTimings::TimerRef PTimer           = IpplTimings::getTimer("pushVelocity");
        static IpplTimings::TimerRef RTimer           = IpplTimings::getTimer("pushPosition");
        static IpplTimings::TimerRef updateTimer      = IpplTimings::getTimer("update");
        static IpplTimings::TimerRef domainDecomposition = IpplTimings::getTimer("loadBalance");

        static IpplTimings::TimerRef runBinnedSolverT = IpplTimings::getTimer("runBinnedSolver");
        //static IpplTimings::TimerRef GenAdaptiveHistogram = IpplTimings::getTimer("genAdaptiveHistogram");
        //static IpplTimings::TimerRef FullRebin128 = IpplTimings::getTimer("FullRebin128");
        static IpplTimings::TimerRef TotalBinningTimer = IpplTimings::getTimer("TotalBinning");

        double dt                               = this->dt_m;
        std::shared_ptr<ParticleContainer_t> pc = this->pcontainer_m;
        std::shared_ptr<FieldContainer_t> fc    = this->fcontainer_m;

        IpplTimings::startTimer(PTimer);
        pc->P = pc->P - 0.5 * dt * pc->E;
        IpplTimings::stopTimer(PTimer);

        // drift
        IpplTimings::startTimer(RTimer);
        pc->R = pc->R + dt * pc->P;
        IpplTimings::stopTimer(RTimer);
        
        // Since the particles have moved spatially update them to correct processors
        IpplTimings::startTimer(updateTimer);
        if (samplingMethod_m == "Flattop") emitFlattop(this->it_m * 1.0, 1.0);
        pc->update();
        IpplTimings::stopTimer(updateTimer);
        
        size_type totalP        = this->totalP_m;
        int it                  = this->it_m;
        bool isFirstRepartition = false;
        if (this->loadbalancer_m->balance(totalP, it + 1)) {
                IpplTimings::startTimer(domainDecomposition);
                auto* mesh = &fc->getRho().get_mesh();
                auto* FL = &fc->getFL();
                this->loadbalancer_m->repartition(FL, mesh, isFirstRepartition);
                IpplTimings::stopTimer(domainDecomposition);
        }
        
        IpplTimings::startTimer(TotalBinningTimer);

        // TODO: binning
        /*std::cout << "=============== Starting Reduction Experiment parallel_reduce ===============" << std::endl;
        for (size_t i = 2; i < 129; i++) {
            this->bins_m->initLimits();
            this->bins_m->setCurrentBinCount(i);
            this->bins_m->assignBinsToParticles();
            this->bins_m->initLocalHisto(ParticleBinning::HistoReductionMode::ParallelReduce);
            // this->bins_m->doFullRebin(i, true, HistoReductionMode::Standard);
        }

        std::cout << "=============== Starting Reduction Experiment team_for ===============" << std::endl;
        for (size_t i = 2; i < 256; i++) {
            this->bins_m->initLimits();
            this->bins_m->setCurrentBinCount(i);
            this->bins_m->assignBinsToParticles();
            this->bins_m->initLocalHisto(ParticleBinning::HistoReductionMode::TeamBased);
            // this->bins_m->doFullRebin(i, true, HistoReductionMode::Standard);
        }

        std::cout << "=============== Starting Reduction Experiment Atomics ===============" << std::endl;
        for (size_t i = 2; i < 256; i++) {
            this->bins_m->initLimits();
            this->bins_m->setCurrentBinCount(i);
            this->bins_m->assignBinsToParticles();

            // Now do the reduction. init kokkos view
            auto bins = this->pcontainer_m->Bin.getView();
            Kokkos::View<size_type*> to_reduce("to_reduce", i);
            auto start = std::chrono::high_resolution_clock::now(); // TODO: remove
            Kokkos::parallel_for("atomic", this->pcontainer_m->getLocalNum(), KOKKOS_LAMBDA(const size_t j) {
                Kokkos::atomic_increment(&to_reduce(bins(j)));
            });
            auto end = std::chrono::high_resolution_clock::now(); // TODO: remove
            long long duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); // TODO: remove
            std::cout << "executeAtomicReduction;" << this->pcontainer_m->getLocalNum() << ";" << i << ";" << duration << std::endl; // TODO: remove
        }*/

        /*for (size_t i = 2; i < 257; i++) {
            //this->bins_m->initLimits();
            //this->bins_m->assignBinsToParticles();
            //this->bins_m->initLocalHisto(ParticleBinning::HistoReductionMode::TeamBased);
            //this->bins_m->doFullRebin(10);
            this->bins_m->doFullRebin(i);
            this->bins_m->sortContainerByBin();
        }*/

        //IpplTimings::startTimer(FullRebin128);
        this->bins_m->doFullRebin(128);
        this->bins_m->print(); // For debugging...

        // this->bins_m->doFullRebin(10, true, HistoReductionMode::Standard);
        //IpplTimings::stopTimer(FullRebin128);
        //this->bins_m->print(); // for debugging...
        
        this->bins_m->sortContainerByBin(); // sort particles after creating bins for scatter() operation inside LeapFrogStep 
        // this->bins_m->initLocalHisto(HistoReductionMode::Standard);

        //IpplTimings::startTimer(GenAdaptiveHistogram);
        this->bins_m->genAdaptiveHistogram(); // merge bins with width/N_part ratio of 1.0
        IpplTimings::stopTimer(TotalBinningTimer);
        //IpplTimings::stopTimer(GenAdaptiveHistogram);

        this->bins_m->print(); // For debugging...


        IpplTimings::startTimer(runBinnedSolverT);
        E_tmp = 0.0; // reset temporary field
        //runBinnedSolver();
        IpplTimings::stopTimer(runBinnedSolverT);


        // kick
        IpplTimings::startTimer(PTimer);
        pc->P = pc->P - 0.5 * dt * pc->E;
        IpplTimings::stopTimer(PTimer);
    }

    void runBinnedSolver() {
        static IpplTimings::TimerRef PerBinSolver = IpplTimings::getTimer("PerBinSolver");
        static IpplTimings::TimerRef CombineEFieldTimer = IpplTimings::getTimer("CombineEField");
        /*
         * Strategy:
         * Initialize E field to 0.
         * for every bin:
         *      1. par2grid() for all charges(bin) onto rho_m
         *      2. Solve Poisson equation, only for Base:SOL
         *      3. Add grad(phi_m)*\gamma to E field
         * grid2par()
         */
        Inform msg("runBinnedSolver");
        //static IpplTimings::TimerRef SolveTimer = IpplTimings::getTimer("solve");
        using binIndex_t       = typename ParticleContainer_t::bin_index_type;
        using binIndexView_t   = typename ippl::ParticleAttrib<binIndex_t>::view_type;

        //this->bins_m->print();

        // Defines used views
        std::shared_ptr<ParticleContainer_t> pc = this->pcontainer_m;
        std::shared_ptr<FieldContainer_t> fc    = this->fcontainer_m;
        view_type viewP                         = pc->P.getView();
        binIndexView_t bin                      = pc->Bin.getView();

        for (binIndex_t i = 0; i < this->bins_m->getCurrentBinCount(); ++i) {
            if (this->bins_m->getNPartInBin(i) == 0) { continue; }
            // Scatter only for current bin index
            this->par2gridPerBin(i);

            // Run solver: obtains phi_m only for what was scattered in the previous step
            IpplTimings::startTimer(PerBinSolver);
            this->fsolver_m->runSolver();
            IpplTimings::stopTimer(PerBinSolver);

            IpplTimings::startTimer(CombineEFieldTimer);
            E_tmp = E_tmp + this->bins_m->LTrans(fc->getE(), i);
            IpplTimings::stopTimer(CombineEFieldTimer);
        }

        // TODO: remove. A little debug output:
        /*{
            this->par2grid();
            this->fsolver_m->runSolver();
            msg << "E assigned from " << this->bins_m->getCurrentBinCount() << " bins. Norm = " << ParticleBinning::vnorm(E_tmp) << endl;
            msg << "              Single bin norm = " << ParticleBinning::vnorm(fc->getE()) << endl;
        }*/

        // gather E field from locally built up E_m
        gather(pc->E, E_tmp, this->pcontainer_m->R); 
        /*
        Alternative: don't use a temporary field, but directly gather the field for every particles inside the loop.
        Problem: gather routine is way more expensive than simple additions on a temporary field instance. 
        However: if the field is big and memory is crucial, one might consider this approach and remove E_tmp.
        */
        msg << "Field gathered" << endl;
    }

    void dump() override {
        static IpplTimings::TimerRef dumpDataTimer = IpplTimings::getTimer("dumpData");
        IpplTimings::startTimer(dumpDataTimer);
        dumpLandau(this->fcontainer_m->getE().getView());
        IpplTimings::stopTimer(dumpDataTimer);
    }

    template <typename View>
    void dumpLandau(const View& Eview) {
        const int nghostE = this->fcontainer_m->getE().getNghost();

        using index_array_type = typename ippl::RangePolicy<Dim>::index_array_type;
        double localEx2 = 0, localExNorm = 0;
        ippl::parallel_reduce(
            "Ex stats", ippl::getRangePolicy(Eview, nghostE),
            KOKKOS_LAMBDA(const index_array_type& args, double& E2, double& ENorm) {
                // ippl::apply<unsigned> accesses the view at the given indices and obtains a
                // reference; see src/Expression/IpplOperations.h
                double val = ippl::apply(Eview, args)[0];
                double e2  = Kokkos::pow(val, 2);
                E2 += e2;

                double norm = Kokkos::fabs(ippl::apply(Eview, args)[0]);
                if (norm > ENorm) {
                    ENorm = norm;
                }
            },
            Kokkos::Sum<double>(localEx2), Kokkos::Max<double>(localExNorm));

        double globaltemp = 0.0;
        ippl::Comm->reduce(localEx2, globaltemp, 1, std::plus<double>());

        double fieldEnergy =
            std::reduce(this->fcontainer_m->getHr().begin(), this->fcontainer_m->getHr().end(), globaltemp, std::multiplies<double>());

        double ExAmp = 0.0;
        ippl::Comm->reduce(localExNorm, ExAmp, 1, std::greater<double>());

        if (ippl::Comm->rank() == 0) {
            std::stringstream fname;
            fname << "data/FieldLandau_";
            fname << ippl::Comm->size();
            fname << "_manager";
            fname << ".csv";
            Inform csvout(NULL, fname.str().c_str(), Inform::APPEND);
            csvout.precision(16);
            csvout.setf(std::ios::scientific, std::ios::floatfield);
            if ( std::fabs(this->time_m) < 1e-14 ) {
                csvout << "time, Ex_field_energy, Ex_max_norm" << endl;
            }
            csvout << this->time_m << " " << fieldEnergy << " " << ExAmp << endl;
        }
        ippl::Comm->barrier();
    }

};
#endif
