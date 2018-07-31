#include "PlasmaPIC.h"

#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include "Physics/Physics.h"

#include <string>
#include <fstream>

#include "../Solver.h"

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
    
    this->initAmr_m();
    
    this->initBunch_m();
    
    this->initDistribution_m();
    
    this->initSolver_m();
}


PlasmaPIC::~PlasmaPIC() {
    if ( solver_mp )
        delete solver_mp;
}


void PlasmaPIC::execute(Inform& msg) {
    
    amropal_m->setBunch( bunch_m.get() );
    
    for (int i = 0; i <= amropal_m->finestLevel() &&
                    i < amropal_m->maxLevel(); ++i)
    {
        amropal_m->regrid(i /*lbase*/,tcurrent_m);
    }
    
    while ( tcurrent_m < tstop_m ) {
        
        if ( tcurrent_m + dt_m > tstop_m )
            dt_m = tstop_m - tcurrent_m;
        
        msg << "At time " << tcurrent_m
            << " with timestep " << dt_m << ". ";
        
        this->integrate_m();
        
        this->dump_m();
        
        tcurrent_m += dt_m;
        
        msg << "Done." << endl;
    }
}


void PlasmaPIC::parseParticleInfo_m() {
    amrex::ParmParse pp("particle");
    
    pp.get("test", test_m);
    
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
    
    pp.get("wavenumber", wavenumber_m);
    
    left_m = Vector_t(D_DECL(0.0, 0.0, 0.0));
    
    using Physics::two_pi;
    double length = two_pi / wavenumber_m;
    right_m = Vector_t(D_DECL(length, length, length));
    
    int nx = 0;
    pp.get("nx", nx);
    bNx_m = Vector<int>(D_DECL(nx, nx, nx));
    
    std::string dir = test_m + "-data-grid-" +
                      AMREX_D_TERM(std::to_string(bNx_m[0]),
                                   + "-" + std::to_string(bNx_m[1]),
                                   + "-" + std::to_string(bNx_m[2]));
    
    boost::filesystem::path dir_m(dir);
    if ( Ippl::myNode() == 0 )
        boost::filesystem::create_directory(dir_m);
}


void PlasmaPIC::initAmr_m() {
    
    amrex::ParmParse pp("amr");
    pp.add("max_grid_size", maxgrid_m);
    
    pp.add("max_level", nlevel_m);
    
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


void PlasmaPIC::initBunch_m() {
    
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
}


void PlasmaPIC::initDistribution_m() {
    
    const BoxArray ba = amropal_m->boxArray(0);
    const DistributionMapping& dmap = amropal_m->DistributionMap(0);
    
    Vector_t length = right_m - left_m;
    Vector_t hx = length / pNx_m;
    Vector_t hv = length / pNv_m;
    
    std::size_t ip = 0;
    
    for (amrex::MFIter mfi(ba, dmap); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        
        for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
            for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
                amrex::IntVect iv(D_DECL(i, j, 0));
                
                double x = (iv[0] + 0.5) * hx[0] + left_m[0];
                double y = (iv[1] + 0.5) * hx[1] + left_m[1];
                
                for (int ii = 0; ii < pNv_m[0]; ii++) {
                    for (int jj = 0; jj < pNv_m[1]; jj++) {
                        double vx = (ii + 0.5)*hv[0] + vmin_m[0];
                        double vy = (jj + 0.5)*hv[1] + vmin_m[1];
                        double v2 = vx * vx + vy * vy;
                        double f = (1.0 / (2.0 * Physics::pi)) * std::exp(-0.5 * v2) *
                                    ( 1.0 + alpha_m *
                                      std::cos(wavenumber_m * x) * std::cos(wavenumber_m * y)
                                    );
                        double m = hx[0] * hv[0] * hx[1] * hv[1] * f;
                        double q = -m;
                        
                        if ( m > threshold_m ) {
                            
                            bunch_m->create(1);
                            
                            bunch_m->R[ip](0) = x;
                            bunch_m->R[ip](1) = y;                            
                            bunch_m->P[ip](0) = vx;
                            bunch_m->P[ip](1) = vy;
                            bunch_m->qm[ip] = q;
                            bunch_m->mass[ip] = m;
                            
                            ++ip;
                        }
                    }
                }
            }
        }
    }
    bunch_m->update();
}


void PlasmaPIC::initSolver_m() {
    amrex::ParmParse pp("solver");
    
    std::string type = "AMReX_FMG";
    
    pp.get("type", type);
    
    if ( type == "AMReX_FMG" )
        this->initAmrexFMG_m();
#ifdef HAVE_AMR_MG_SOLVER
    else if ( type == "AMR_MG" )
        this->initAmrMG();
#endif
}


void PlasmaPIC::initAmrexFMG_m() {
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
void PlasmaPIC::initAmrMG() {
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
    // RK-4
    this->solve_m();
    
    bunch_m->k1 = bunch_m->E;
    
    bunch_m->R += 0.5 * bunch_m->P * dt_m + 0.125 * bunch_m->k1 * dt_m * dt_m;
    
    this->applyPeriodicity_m();
    
    this->solve_m();
    
    bunch_m->k2 = bunch_m->E;
    
    bunch_m->R += 0.5 * bunch_m->P * dt_m + (0.5 * bunch_m->k2 - 0.125 * bunch_m->k1) * dt_m * dt_m;
    
    this->applyPeriodicity_m();
    
    this->solve_m();
    
    bunch_m->R = bunch_m->R - bunch_m->P * dt_m - 0.5 * bunch_m->k2 * dt_m * dt_m;
    
    bunch_m->R += bunch_m->P * dt_m + (1.0 / 6.0 ) * (bunch_m->k1 + 2.0 * bunch_m->k2) * dt_m * dt_m;
    bunch_m->P += (1.0 / 6.0) * (bunch_m->k1 + 4.0 * bunch_m->k2 + bunch_m->E) * dt_m;
    
    this->applyPeriodicity_m();
    
//     // RK-2
//     
//     this->solve_m();
//     
//     bunch_m->k1 = bunch_m->E;
//     
//     bunch_m->R += bunch_m->P * dt_m;
//     
//     this->applyPeriodicity_m();
//     
//     this->solve_m();
//     
//     bunch_m->R += 0.5 * bunch_m->k1 * dt_m * dt_m;
//     bunch_m->P += 0.5 * (bunch_m->k1 + bunch_m->E) * dt_m;
//     
//     this->applyPeriodicity_m();
    
    
    
    
    
    
    
//     /* Leap-Frog
//      * 
//      * v_{i+1/2} = v_i + a_i * dt / 2    (a_i = F_i = q * E_i )
//      * x_{i+1} = x_i + v_{i+1/2} * dt
//      * v_{i+1} = v_{i+1/2} + a_{i+1} * dt / 2
//      */
//     
//     this->solve_m();
//     
//     bunch_m->P = bunch_m->P + 0.5 * dt_m * bunch_m->qm / bunch_m->mass * bunch_m->E;
//     bunch_m->R = bunch_m->R + dt_m * bunch_m->P;
//     
//     this->applyPeriodicity_m();
//     
//     this->solve_m();
//     
//     bunch_m->P = bunch_m->P + 0.5 * dt_m * bunch_m->qm / bunch_m->mass * bunch_m->E;
}


void PlasmaPIC::applyPeriodicity_m()
{
    for (std::size_t j = 0; j < bunch_m->getLocalNum(); ++j) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if ( bunch_m->R[j](d) > right_m[d] )
                bunch_m->R[j](d) = bunch_m->R[j](d) - right_m[d];
            else if ( bunch_m->R[j](d) < 0.0 )
                bunch_m->R[j](d) = bunch_m->R[j](d) + right_m[d];
        }
    }
    
    if ( amropal_m->maxLevel() > 0 )
        for (int k = 0; k <= amropal_m->finestLevel() && k < amropal_m->maxLevel(); ++k)
            amropal_m->regrid(k /*lbase*/, tcurrent_m);
    else
        bunch_m->update();
}


void PlasmaPIC::electricField_m(double& field_energy, double& amplitude)
{
    std::vector<amrex::MultiFab> cp_efield(efield_m.size());
    
    for (uint lev = 0; lev < efield_m.size(); ++lev) {
        
        cp_efield[lev].define(efield_m[lev]->boxArray(),
                              efield_m[lev]->DistributionMap(), AMREX_SPACEDIM, 1);
        
        cp_efield[lev].setVal(0.0, 1);
        
        amrex::MultiFab::Copy(cp_efield[lev], *efield_m[lev], 0, 0, AMREX_SPACEDIM, 1);
    }
    
    field_energy = 0.0;
    amplitude = -1.0;
    
    for (uint lev = 0; lev < cp_efield.size() - 1; ++lev) {
        // get boxarray with refined cells
        amrex::BoxArray ba = cp_efield[lev].boxArray();
        ba.refine(2);
        ba = amrex::intersect(cp_efield[lev+1].boxArray(), ba);
        ba.coarsen(2);
        for ( uint b = 0; b < ba.size(); ++b)
            for (int c = 0; c < AMREX_SPACEDIM; ++c)
                cp_efield[lev].mult(0.0, ba[b], c, AMREX_SPACEDIM, 1);
    }
    
    const auto& geom = amropal_m->Geom();
    
    for (uint lev = 0; lev < efield_m.size(); ++lev) {
            amrex::MultiFab::AddProduct(cp_efield[lev],
                                        cp_efield[lev], 0,
                                        cp_efield[lev], 0,
                                        0, AMREX_SPACEDIM, 0);
        for (int c = 0; c < AMREX_SPACEDIM; ++c) {
            field_energy += cp_efield[lev].sum(c);
            amplitude = std::max(amplitude, std::sqrt(cp_efield[lev].max(c)));
        }
        
        field_energy *= 0.5 * AMREX_D_TERM(  geom[lev].CellSize(0),
                                           * geom[lev].CellSize(1),
                                           * geom[lev].CellSize(2) );
            
    }
}


void PlasmaPIC::dump_m() {
    // kinetic energy
    double ekin = 0.5 * sum( dot(bunch_m->P, bunch_m->P) );
    
    double epot = 0;
    for (unsigned i=0; i<bunch_m->getLocalNum(); ++i) {
        epot += 0.5 * (bunch_m->qm[i] *  bunch_m->phi[i]);
    }
    
    double glob_epot = 0.0;
    MPI_Reduce(&epot, &glob_epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    double ampl = 0.0;
    double field_energy = 0.0;
    electricField_m(field_energy, ampl);
    
    if(Ippl::myNode()==0) {
        std::ofstream csvout;
        csvout.precision(10);
        csvout.setf(std::ios::scientific, std::ios::floatfield);

        std::stringstream fname;
        fname << (dir_m / "energy.csv").string();

        // open a new data file for this iteration
        // and start with header
        csvout.open(fname.str().c_str(), std::ios::out | std::ofstream::app);
        
        if (tcurrent_m < dt_m) {
            csvout << "it, field energy, kinetic energy, "
                   << "total energy, potential energy" << std::endl;
        }
        
        csvout << tcurrent_m << ", "
               << field_energy << ","
               << ekin << ","
               << glob_epot + ekin << "," 
               << epot << std::endl;
        
        csvout.close();
        
        
        fname.str("");
        fname << (dir_m / "amplitude.csv").string();
        
        csvout.open(fname.str().c_str(), std::ios::out | std::ofstream::app);
        if (tcurrent_m < dt_m){
                csvout << "it,max(|Ez|)" << std::endl;
        }
        csvout << tcurrent_m << ", "
                << ampl << std::endl;
        csvout.close();
    }
}