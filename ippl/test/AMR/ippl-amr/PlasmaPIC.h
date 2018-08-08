#ifndef PLASMA_PIC_H
#define PLASMA_PIC_H

#include <memory>

#include <boost/filesystem.hpp>

#include "Ippl.h"

#include <AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "../AmrOpal.h"

#include "../AbstractSolver.h"

#include "PhaseDist.h"

using amrex::Vector;
using amrex::Geometry;
using amrex::MultiFab;
using amrex::BoxArray;
using amrex::DistributionMapping;

class PlasmaPIC
{
public:
    
    typedef Vektor<double, AMREX_SPACEDIM> Vector_t;
    typedef Vektor<int, AMREX_SPACEDIM> iVector_t;
    typedef AmrOpal::amrplayout_t amrplayout_t;
    typedef AmrOpal::amrbase_t amrbase_t;
    typedef AmrOpal::amrbunch_t amrbunch_t;
    typedef std::unique_ptr<amrbunch_t> amrbunch_p;
    typedef AmrOpal amropal_t;
    typedef std::unique_ptr<amropal_t> amropal_p;
    typedef AbstractSolver* solver_p;
    typedef Vector<std::unique_ptr<MultiFab> > container_t;
    
public:
    PlasmaPIC();
    
    ~PlasmaPIC();
    
    void execute(Inform& msg);
    
private:
    
    void parseParticleInfo_m();
    void parseBoxInfo_m();
    
    void initAmr_m();
    
    void initBunch_m();
    
    void initDistribution_m();
    
    void initSolver_m();
    
    void initAmrexFMG_m();
    
#ifdef HAVE_AMR_MG_SOLVER
    void initAmrMG();
#endif
    
    void depositCharge_m();
    
    void solve_m();
    
    void integrate_m();
    
    void applyPeriodicity_m();
    
    void electricField_m(double& field_energy, double& amplitude);
    
    void dump_m();
    
    void volWeightedSum_m(double& sum, double& vol);
    
    void dumpFields_m();
    
    void deposit2D_m(int step);
    
private:
    int maxgrid_m;
    int blocking_factor_m;
    double dt_m;
    double tstop_m;
    double tcurrent_m;
    int nlevel_m;
    
    double threshold_m;
    std::string test_m;
    
    boost::filesystem::path dir_m;
    
    Vector_t vmax_m;
    Vector_t vmin_m;
    
    double alpha_m;
    
    iVector_t pNx_m;
    iVector_t pNv_m;
    
    double wavenumber_m;
    double boxlength_m;
    Vector_t left_m;
    Vector_t right_m;
    Vector<int> bNx_m;
    
    amropal_p amropal_m;
    amrbunch_p bunch_m;
    solver_p solver_mp;
    
    
    container_t rho_m;
    container_t phi_m;
    container_t efield_m;
    
    
    PhaseDist pd_m;
};

#endif
