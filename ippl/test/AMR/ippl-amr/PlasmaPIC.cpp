#include "PlasmaPIC.h"

#include <boost/filesystem.hpp>

#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include "Physics/Physics.h"

#include "../Solver.h"

#include "../Distribution.h"

#include "../helper_functions.h"

#ifdef HAVE_AMR_MG_SOLVER
    #include "../trilinos/AmrMultiGrid.h"
#endif


PlasmaPIC::PlasmaPIC() : tcurrent_m(0.0) {
    
    amrex::ParmParse pp;
    
    pp.get("dt", dt_m);
    pp.get("tstop", tstop_m);
    
    this->parseParticleInfo_m();
    
    this->parseBoxInfo_m();
    
    this->initializeAmr_m();
    
    this->initializeBunch_m();
    
    this->initializeSolver_m();
}


PlasmaPIC::~PlasmaPIC() {
    if ( solver_mp )
        delete solver_mp;
}


void PlasmaPIC::execute() {
    
    amropal_m->setBunch( bunch_m.get() );
    
    for (int i = 0; i <= amropal_m->finestLevel() &&
                    i < amropal_m->maxLevel(); ++i)
    {
        amropal_m->regrid(i /*lbase*/, 0.0 /*time*/);
    }
    
    while ( tcurrent_m < tstop_m ) {
        
        if ( tcurrent_m + dt_m > tstop_m )
            dt_m = tstop_m - tcurrent_m;
        
        this->integrate_m();
    }
}


void PlasmaPIC::parseParticleInfo_m() {
    amrex::ParmParse pp("particle");
    
    pp.get("test", test_m);
    
    dir_m = test_m + "-data-grid-" +
            AMREX_D_TERM(std::to_string(bNx_m[0]),
                         + "-" + std::to_string(bNx_m[1]),
                         + "-" + std::to_string(bNx_m[2]));
    
    boost::filesystem::path dir(dir_m);
    if ( Ippl::myNode() == 0 )
        boost::filesystem::create_directory(dir);
    dir_m = dir.string();
    
    pp.get("alpha", alpha_m);
    pp.get("threshold", threshold_m);
    
    int nx = 0;
    pp.get("nx", nx);
    pNx_m = iVector_t(D_DECL(nx, nx, nx));
    
    int nv = 0;
    pp.get("nv", nv);
    pNv_m = iVector_t(D_DECL(nv, nv, nv));
    
    double vmin = 0.0;
    pp.get("vmin", vmin);
    vmin_m = Vector_t(D_DECL(vmin, vmin, vmin));
    
    double vmax = 0.0;
    pp.get("vmax", vmax);
    vmax_m = Vector_t(D_DECL(vmax, vmax, vmax));
}


void PlasmaPIC::parseBoxInfo_m() {
    amrex::ParmParse pp("box");
    
    pp.get("maxgrid", maxgrid_m);
    pp.get("blocking_factor", blocking_factor_m);
    pp.get("nlevel", nlevel_m);
    pp.get("boxlength", boxlength_m);
    
    pp.get("wavenumber", wavenumber_m);
    
    left_m = Vector_t(D_DECL(0.0, 0.0, 0.0));
    
    using Physics::two_pi;
    double length = two_pi * boxlength_m / wavenumber_m;
    right_m = Vector_t(D_DECL(length, length, length));
    
    int nx = 0;
    pp.get("nx", nx);
    bNx_m = Vector<int>(D_DECL(nx, nx, nx));
}


void PlasmaPIC::initializeAmr_m() {
    
    amrex::ParmParse pp("amr");
    pp.add("max_grid_size", maxgrid_m);
    
    amrex::Vector<int> bf(nlevel_m, blocking_factor_m);
    pp.addarr("blocking_factor", bf);
    
    amrex::Vector<int> error_buf(nlevel_m, 0);
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
    amrex::ParmParse geom("geometry");
    amrex::Vector<int> is_per = { D_DECL(1, 1, 1) };
    geom.addarr("is_periodic", is_per);
    
    
    std::vector<int> rr(nlevel_m);
    for (int i = 0; i < nlevel_m; ++i) {
        rr[i] = 2;
    }
    
    amrex::RealBox domain;
    
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        domain.setLo(i, left_m[i]);
        domain.setHi(i, right_m[i]);
    }
    
    amropal_m.reset( new amropal_t(&domain,
                                   nlevel_m - 1,
                                   bNx_m,
                                   0 /* cartesian */,
                                   rr) );
    
    amropal_m->setTagging(AmrOpal::kChargeDensity);
    amropal_m->setCharge(1.0e-14);
}


void PlasmaPIC::initializeBunch_m() {
    
    const auto& geom = amropal_m->Geom();
    const auto& dmap = amropal_m->DistributionMap();
    const auto& ba   = amropal_m->boxArray();
    
    Vector<int> rr(nlevel_m);
    for (int i = 0; i < nlevel_m; ++i) {
        rr[i] = 2;
    }
    
    amrplayout_t* playout = new amrplayout_t(geom, dmap, ba, rr);
    
    bunch_m.reset( new amrbunch_t() );
    bunch_m->initialize(playout);
    bunch_m->initializeAmr(); // add attributes: level, grid
    
    bunch_m->setAllowParticlesNearBoundary(true);
    
    Distribution dist;
    
//     dist.plasma(amropal_m,
//                 bunch_m,
//                 vmin_m, vmax_m,
//                 pNx_m, pNv_m,
//                 test_m);
    
    dist.injectBeam( *(bunch_m.get()) );
    
    bunch_m->update();
}


void PlasmaPIC::initializeSolver_m() {
    amrex::ParmParse pp("solver");
    
    std::string type = "AMReX_FMG";
    
    pp.get("type", type);
    
    if ( type == "AMReX_FMG" )
        this->initializeAmrexFMG_m();
#ifdef HAVE_AMR_MG_SOLVER
    else if ( type == "AMR_MG" )
        this->initializeAmrMG();
#endif
}


void PlasmaPIC::initializeAmrexFMG_m() {
    int maxiter = 100;
    int maxiter_b = 100;
    int verbose = 0;
    bool usecg = true;
    double bottom_solver_eps = 1.0e-4;
    int max_nlevel = 1024;
    
    /* MG_SMOOTHER_GS_RB  = 1
     * MG_SMOOTHER_JACOBI = 2
     * MG_SMOOTHER_MINION_CROSS = 5
     * MG_SMOOTHER_MINION_FULL = 6
     * MG_SMOOTHER_EFF_RB = 7
     */
    int smoother = 1;
    
    // #smoothings at each level on the way DOWN the V-cycle
    int nu_1 = 2;
    
    // #smoothings at each level on the way UP the V-cycle
    int nu_2 = 2;
    
    // #smoothings before and after the bottom solver
    int nu_b = 0;
    
    // #smoothings
    int nu_f = 8;

    /* MG_FCycle = 1  (full multigrid)
     * MG_WCycle = 2
     * MG_VCycle = 3
     * MG_FVCycle = 4
     */
    int cycle = 1;
    
    bool cg_solver = true;
    
    /* if cg_solver == true:
     * - BiCG --> 1
     * - CG --> 2
     * - CABiCG --> 3
     * 
     * else if cg_solver == false
     * - CABiCG is taken
     */
    int bottom_solver = 1;
    
    amrex::ParmParse pp("mg");

    pp.add("maxiter", maxiter);
    pp.add("maxiter_b", maxiter_b);
    pp.add("nu_1", nu_1);
    pp.add("nu_2", nu_2);
    pp.add("nu_b", nu_b);
    pp.add("nu_f", nu_f);
    pp.add("v"   , verbose);
    pp.add("usecg", usecg);
    pp.add("cg_solver", cg_solver);

    pp.add("rtol_b", bottom_solver_eps);
    pp.add("numLevelsMAX", max_nlevel);
    pp.add("smoother", smoother);
    pp.add("cycle_type", cycle); // 1 -> F, 2 -> W, 3 -> V, 4 -> F+V
    //
    // The C++ code usually sets CG solver type using cg.cg_solver.
    // We'll allow people to also use mg.cg_solver but pick up the former as well.
    //
    if (!pp.query("cg_solver", cg_solver))
    {
        amrex::ParmParse pp("cg");

        pp.add("cg_solver", cg_solver);
    }

    pp.add("bottom_solver", bottom_solver);
    
    solver_mp = new Solver();
}


#ifdef HAVE_AMR_MG_SOLVER
void PlasmaPIC::initializeAmrMG() {
    std::string interp = "PC";
    std::string norm = "LINF";
    std::string bs = "SA";
    std::string smoother = "GS";
    std::string prec = "NONE";
    bool rebalance = true;
    std::size_t nsweeps = 12;
    solver_mp = new AmrMultiGrid(*(amropal_m.get()),
                                 bs, prec, rebalance,
                                 "periodic", "periodic",
                                 "periodic", smoother,
                                 nsweeps, interp, norm,
                                 0, 0);
}
#endif

void PlasmaPIC::depositCharge_m() {
    int base_level   = 0;
    int finest_level = amropal_m->finestLevel();
    
    container_t partMF(nlevel_m);
    for (int lev = 0; lev < nlevel_m; lev++) {
        const BoxArray& ba = amropal_m->boxArray()[lev];
        const DistributionMapping& dmap = amropal_m->DistributionMap(lev);
        partMF[lev].reset(new MultiFab(ba, dmap, 1, 2));
        partMF[lev]->setVal(0.0, 2);
   }
    
    bunch_m->AssignDensityFort(bunch_m->qm, partMF, base_level, 1, finest_level);
    
    for (int lev = 0; lev < nlevel_m; ++lev) {
        amrex::MultiFab::Copy(*rho_m[lev], *partMF[lev], 0, 0, 1, 0);
        
        if ( rho_m[lev]->contains_nan() || rho_m[lev]->contains_inf() )
            throw std::runtime_error("\033[1;31mError: NANs or INFs on charge grid.\033[0m");
    }
    
    
    
    const Vector<Geometry>& geom = amropal_m->Geom();
    // Check charge conservation
    double totCharge = totalCharge(rho_m, finest_level, geom);
    double totCharge_composite = totalCharge_composite(rho_m, finest_level, geom);
    
    if ( Ippl::myNode() == 0 ) {
        std::cout << "Total Charge (computed):  " << totCharge << " C" << std::endl
                  << "Total Charge (composite): " << totCharge_composite << " C" << std::endl;
    }
    
    //subtract the background charge of the ions
    for (int i = 0; i <= finest_level; ++i) {
        rho_m[i]->plus(1.0, 0, 1);
    }
    
    // \Delta\phi = -\rho
    for (int i = 0; i <= finest_level; ++i) {
        rho_m[i]->mult(-1.0, 0, 1);
    }
}


void PlasmaPIC::solve_m() {
    static IpplTimings::TimerRef solvTimer = IpplTimings::getTimer("solve");
    
    /*
     * reset grid data
     */
    rho_m.resize(nlevel_m);
    phi_m.resize(nlevel_m);
    efield_m.resize(nlevel_m);
    
    // Define the density on level 0 from all particles at all levels                                                                                                                                            
    int base_level   = 0;
    int finest_level = amropal_m->finestLevel();
    
    for (int lev = 0; lev < nlevel_m; ++lev) {
        initGridData(rho_m, phi_m, efield_m,
                     amropal_m->boxArray()[lev],
                     amropal_m->DistributionMap(lev), lev);
    }
    
    
    /*
     *  deposit charge
     */
    this->depositCharge_m();
    
    // normalize each level
    double l0norm = rho_m[finest_level]->norm0(0);
    for (int i = 0; i <= finest_level; ++i) {
        rho_m[i]->mult(1.0 / l0norm, 0, 1);
    }
    
    IpplTimings::startTimer(solvTimer);
    
    solver_mp->solve(amropal_m,
                     rho_m, phi_m, efield_m,
                     base_level, finest_level, false);
    
    IpplTimings::stopTimer(solvTimer);
    
    // undo normalization
    for (int i = 0; i <= finest_level; ++i) {
        phi_m[i]->mult(l0norm, 0, 1);
    }
    
    bunch_m->phi = 0;
    
    bunch_m->InterpolateFort(bunch_m->phi, phi_m, base_level, finest_level);
    
    // undo scale
    for (int i = 0; i <= finest_level; ++i)
        efield_m[i]->mult(l0norm, 0, AMREX_SPACEDIM);
    
    bunch_m->E = 0;
    
    bunch_m->InterpolateFort(bunch_m->E, efield_m, base_level, finest_level);
}


void PlasmaPIC::integrate_m() {
    
    
}
