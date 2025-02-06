#ifndef IPPL_ALPINE_MANAGER_H
#define IPPL_ALPINE_MANAGER_H

#include <memory>

#include "FieldContainer.hpp"
#include "FieldSolver.hpp"
#include "LoadBalancer.hpp"
#include "Manager/BaseManager.h"
#include "Manager/PicManager.h"
#include "ParticleContainer.hpp"
#include "Random/Distribution.h"
#include "Random/InverseTransformSampling.h"
#include "Random/NormalDistribution.h"
#include "Random/Randn.h"

#include "AdaptBins.h" // TODO: Binning

using view_type = typename ippl::detail::ViewType<ippl::Vector<double, Dim>, 1>::view_type;

template <typename T, unsigned Dim>
class AlpineManager
    : public ippl::PicManager<T, Dim, ParticleContainer<T, Dim>, FieldContainer<T, Dim>,
                              LoadBalancer<T, Dim>> {
public:
    using ParticleContainer_t = ParticleContainer<T, Dim>;
    using FieldContainer_t = FieldContainer<T, Dim>;
    using FieldSolver_t= FieldSolver<T, Dim>;
    using LoadBalancer_t= LoadBalancer<T, Dim>;
    using Base= ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>;

    // TODO: Binning
    using BinningSelector_t = typename ParticleBinning::CoordinateSelector<ParticleContainer_t>;
    using AdaptBins_t       = typename ParticleBinning::AdaptBins<ParticleContainer_t, BinningSelector_t>;
    
    using binIndex_t = typename ParticleContainer_t::bin_index_type;
    using binIndexView_t = typename ippl::ParticleAttrib<binIndex_t>::view_type;
    using charge_view_t = typename ippl::ParticleAttrib<double>::view_type;

protected:
    size_type totalP_m;
    int nt_m;
    Vector_t<int, Dim> nr_m;
    double lbt_m;
    std::string solver_m;
    std::string stepMethod_m;
    std::shared_ptr<AdaptBins_t> bins_m; // TODO: Binning

public:
    AlpineManager(size_type totalP_, int nt_, Vector_t<int, Dim>& nr_, double lbt_, std::string& solver_, std::string& stepMethod_)
        : ippl::PicManager<T, Dim, ParticleContainer<T, Dim>, FieldContainer<T, Dim>, LoadBalancer<T, Dim>>()
        , totalP_m(totalP_)
        , nt_m(nt_)
        , nr_m(nr_)
        , lbt_m(lbt_)
        , solver_m(solver_)
        , stepMethod_m(stepMethod_)
        , bins_m(nullptr) {} // TODO: Binning
    ~AlpineManager(){}

protected:
    double time_m;
    double dt_m;
    int it_m;
    Vector_t<double, Dim> kw_m;
    double alpha_m;
    Vector_t<double, Dim> rmin_m;
    Vector_t<double, Dim> rmax_m;
    Vector_t<double, Dim> hr_m;
    double Q_m;
    Vector_t<double, Dim> origin_m;
    bool isAllPeriodic_m;
    bool isFirstRepartition_m;
    ippl::NDIndex<Dim> domain_m;
    std::array<bool, Dim> decomp_m;
    double rhoNorm_m;

public:
    std::shared_ptr<AdaptBins_t> getBins() { return bins_m; } // TODO: Binning
    void setBins(std::shared_ptr<AdaptBins_t> bins) { bins_m = bins; } // TODO: Binning

    size_type getTotalP() const { return totalP_m; }

    void setTotalP(size_type totalP_) { totalP_m = totalP_; }

    int getNt() const { return nt_m; }

    void setNt(int nt_) { nt_m = nt_; }

    const std::string& getSolver() const { return solver_m; }

    void setSolver(const std::string& solver_) { solver_m = solver_; }

    double getLoadBalanceThreshold() const { return lbt_m; }

    void setLoadBalanceThreshold(double lbt_) { lbt_m = lbt_; }

    const std::string& getStepMethod() const { return stepMethod_m; }

    void setStepMethod(const std::string& stepMethod_) { stepMethod_m = stepMethod_; }

    const Vector_t<int, Dim>& getNr() const { return nr_m; }

    void setNr(const Vector_t<int, Dim>& nr_) { nr_m = nr_; }

    double getTime() const { return time_m; }

    void setTime(double time_) { time_m = time_; }

    virtual void dump() { /* default does nothing */ };

    void pre_step() override {
        Inform m("Pre-step");
        m << "Done" << endl;
    }

    void post_step() override {
        // Update time
        this->time_m += this->dt_m;
        this->it_m++;
        // wrtie solution to output file
        this->dump();

        Inform m("Post-step:");
        m << "Finished time step: " << this->it_m << " time: " << this->time_m << endl;
    }

    void grid2par() override { gatherCIC(); }

    void gatherCIC() {
        gather(this->pcontainer_m->E, this->fcontainer_m->getE(), this->pcontainer_m->R);
    }

    void par2grid() override { scatterCIC(); }

    void scatterCIC() {
        Inform m("scatter ");
        this->fcontainer_m->getRho() = 0.0;

        ippl::ParticleAttrib<double> *q = &this->pcontainer_m->q;
        typename Base::particle_position_type *R = &this->pcontainer_m->R;
        Field_t<Dim> *rho               = &this->fcontainer_m->getRho();
        double Q                        = Q_m;
        Vector_t<double, Dim> rmin	= rmin_m;
        Vector_t<double, Dim> rmax	= rmax_m;
        Vector_t<double, Dim> hr        = hr_m;

        scatter(*q, *rho, *R);
        double relError = std::fabs((Q-(*rho).sum())/Q);

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
        (*rho)          = (*rho) / cellVolume;

        rhoNorm_m = norm(*rho);

        // rho = rho_e - rho_i (only if periodic BCs)
        if (this->fsolver_m->getStype() != "OPEN") {
            double size = 1;
            for (unsigned d = 0; d < Dim; d++) {
                size *= rmax[d] - rmin[d];
            }
            *rho = *rho - (Q / size);
        }
    }

    
    /*
     * TODO: binning
     * The following contains extra scatter functions for the scatter-per-bin operations!
     */
    void par2gridPerBin(binIndex_t binIndex) {
        Inform m("scatterPerBin");
        this->fcontainer_m->getRho() = 0.0;
        
        ippl::ParticleAttrib<double> *q          = &this->pcontainer_m->q;
        //charge_view_t viewQ                      = q->getView();
        typename Base::particle_position_type *R = &this->pcontainer_m->R;
        Field_t<Dim> *rho                        = &this->fcontainer_m->getRho();
        Vector_t<double, Dim> rmin	             = rmin_m;
        Vector_t<double, Dim> rmax	             = rmax_m;
        Vector_t<double, Dim> hr                 = hr_m;
        binIndexView_t bin                       = this->pcontainer_m->Bin.getView();
        size_type localParticles                 = this->pcontainer_m->getLocalNum();
        double Q                                 = Q_m * this->bins_m->getNPartInBin(binIndex, true)/totalP_m; // Q_m;

        static IpplTimings::TimerRef binSortingAndScatterT = IpplTimings::getTimer("perBinScatter");
        IpplTimings::startTimer(binSortingAndScatterT);
        scatter(*q, *rho, *R, this->bins_m->getBinIterationPolicy(binIndex), this->bins_m->getHashArray());
        IpplTimings::stopTimer(binSortingAndScatterT);
        double relError = std::fabs((Q-(*rho).sum())/Q);
        
        /*
         * Didn't change anything after here...
         */
        m << relError << endl;

        size_type TotalParticles = 0;
        size_type totalP_tmp = this->pcontainer_m->getTotalNum(); // use this, since not all particles might be emitted

        ippl::Comm->reduce(localParticles, TotalParticles, 1, std::plus<size_type>());

        if (ippl::Comm->rank() == 0) {
            if (TotalParticles != totalP_tmp || relError > 1e-10) {
                m << "Time step: " << it_m << endl;
                m << "Total particles in the sim. " << totalP_tmp << " "
                  << "after update: " << TotalParticles << endl;
                m << "Rel. error in charge conservation: " << relError << endl;
                ippl::Comm->abort();
            }
	    }

	    double cellVolume = std::reduce(hr.begin(), hr.end(), 1., std::multiplies<double>());
        (*rho)            = (*rho) / cellVolume;

        rhoNorm_m = norm(*rho);

        // rho = rho_e - rho_i (only if periodic BCs)
        if (this->fsolver_m->getStype() != "OPEN") {
            double size = 1;
            for (unsigned d = 0; d < Dim; d++) {
                size *= rmax[d] - rmin[d];
            }
            *rho = *rho - (Q / size);
        }
   }
};
#endif
