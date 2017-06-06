/*!
 * @file testUnifSphere.cpp
 * @author Matthias Frey
 * @date 20. Dec. 2016
 * @brief Solve the electrostatic potential for a cube
 *        with Dirichlet boundary condition.
 * @details In this example we initialize 1e6 particles
 *          within a sphere of radius R = 0.005 [m].
 *          We then solve the Poisson equation.\n
 *          Domain: [-0.05 (m), 0.05 (m)]^3
 */
#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <sstream>

#include <AMReX_ParmParse.H>

#include "../MGTSolver.h"
#include "../AmrOpal.h"

#include "../helper_functions.h"

// #include "writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"
#include <random>

typedef AmrOpal::amrplayout_t amrplayout_t;
typedef AmrOpal::amrbase_t amrbase_t;
typedef AmrOpal::amrbunch_t amrbunch_t;

typedef Vektor<double, BL_SPACEDIM> Vector_t;

double domainMapping(amrbase_t& PData,
                     const Vector_t& lold,
                     const Vector_t& uold,
                     const Vector_t& lnew,
                     const Vector_t& unew)
{
    // [lold, uold] --> [lnew, unew]
    Vector_t invdiff = 1.0 / ( uold - lold );
    Vector_t slope = ( unew - lnew ) * invdiff;
    Vector_t intercept = ( uold * lnew - lold * unew ) * invdiff;
    
    for (unsigned int i = 0; i < PData.getLocalNum(); ++i)
        PData.R[i] = slope * PData.R[i] + intercept;
    
    Vector_t dnew = unew - lnew;
    Vector_t dold = uold - lold;
    
    return  std::sqrt( dot(dnew, dnew) / dot(dold, dold) ); //( unew - lnew ) / ( uold - lold );
}

void initSphere(double r, amrbunch_t* bunch, int nParticles) {
    bunch->create(nParticles / Ippl::getNodes());
    
    std::mt19937_64 eng;
    
    if ( Ippl::myNode() )
        eng.seed(42 + Ippl::myNode() );
    
    std::uniform_real_distribution<> ph(-1.0, 1.0);
    std::uniform_real_distribution<> th(0.0, 2.0 * Physics::pi);
    std::uniform_real_distribution<> u(0.0, 1.0);
    
    
    std::string outfile = "amr-particles-level-" + std::to_string(0);
    std::ofstream out(outfile);
    
    long double qi = 4.0 * Physics::pi * Physics::epsilon_0 * r * r / double(nParticles);
    
    for (int i = 0; i < nParticles; ++i) {
        // 17. Dec. 2016,
        // http://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability
        // http://stackoverflow.com/questions/5408276/sampling-uniformly-distributed-random-points-inside-a-spherical-volume
        double phi = std::acos( ph(eng) );
        double theta = th(eng);
        double radius = r * std::cbrt( u(eng) );
        
        double x = radius * std::cos( theta ) * std::sin( phi );
        double y = radius * std::sin( theta ) * std::sin( phi );
        double z = radius * std::cos( phi );
        
        out << x << " " << y << " " << z << std::endl;
        
        bunch->R[i] = Vector_t( x, y, z );    // m
        bunch->qm[i] = qi;   // C
    }
    out.close();
}

void doSolve(AmrOpal& myAmrOpal, amrbunch_t* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& efield,
             const Array<Geometry>& geom,
             const Array<int>& rr, int nLevels,
             Inform& msg)
{
    static IpplTimings::TimerRef solvTimer = IpplTimings::getTimer("solv");
    // =======================================================================                                                                                                                                   
    // 4. prepare for multi-level solve                                                                                                                                                                          
    // =======================================================================

    rhs.resize(nLevels);
    phi.resize(nLevels);
    efield.resize(nLevels);
    

    // Define the density on level 0 from all particles at all levels                                                                                                                                            
    int base_level   = 0;
    int finest_level = myAmrOpal.finestLevel();
    
    for (int lev = 0; lev < nLevels; ++lev) {
        initGridData(rhs, phi, efield,
                     myAmrOpal.boxArray()[lev], myAmrOpal.DistributionMap(lev), lev);
    }
    
   container_t partMF(nLevels);
   for (int lev = 0; lev < nLevels; lev++) {
        const amrex::BoxArray& ba = myAmrOpal.boxArray()[lev];
        const amrex::DistributionMapping& dmap = myAmrOpal.DistributionMap(lev);
        partMF[lev].reset(new amrex::MultiFab(ba, dmap, 1, 2));
        partMF[lev]->setVal(0.0, 2);
   }
    
    bunch->AssignDensityFort(bunch->qm, partMF, base_level, 1, finest_level);
    
    for (int lev = 0; lev < nLevels; ++lev) {
        amrex::MultiFab::Copy(*rhs[lev], *partMF[lev], 0, 0, 1, 0);
    }
    
//     writeScalarField(rhs, *(geom[0].CellSize()), -0.05, "amr-rho_scalar-level-");
    
    // Check charge conservation
    double totCharge = totalCharge(rhs, finest_level, geom);
    
    msg << "Total Charge (computed): " << totCharge << " C" << endl
        << "Vacuum permittivity: " << Physics::epsilon_0 << " F/m (= C/(m V)" << endl;
    
    Real vol = (*(geom[0].CellSize()) * *(geom[0].CellSize()) * *(geom[0].CellSize()) );
    msg << "Cell volume: " << *(geom[0].CellSize()) << "^3 = " << vol << " m^3" << endl;
    
    // eps in C / (V * m)
    double constant = -1.0 / Physics::epsilon_0;  // in [V m / C]
    for (int i = 0; i <=finest_level; ++i) {
        rhs[i]->mult(constant, 0, 1);       // in [V m]
    }
    
    // **************************************************************************                                                                                                                                
    // Compute the total charge of all particles in order to compute the offset                                                                                                                                  
    //     to make the Poisson equations solvable                                                                                                                                                                
    // **************************************************************************                                                                                                                                

    Real offset = 0.;
    
//     for (int lev = base_level; lev <= finest_level; lev++) {
//         phi[lev]->setBndry(0.0);
// //         rhs[lev]->FillBoundary();
//         phi[lev]->FillBoundary();
//     }

    // solve                                                                                                                                                                                                     
    MGTSolver sol;
    IpplTimings::startTimer(solvTimer);
    sol.solve(rhs,            // [V m]
              phi,            // [V m^3]
              efield,       // [V m^2]
              geom);
    
//     for (int i = 0; i <=finest_level; ++i) {
//         rhs[i]->mult(1.0 / constant, 0, 1);
//     }
    
    IpplTimings::stopTimer(solvTimer);
}

void doAMReX(const Vektor<size_t, 3>& nr,
              int nLevels,
              size_t maxBoxSize,
              Inform& msg)
{
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    std::array<double, BL_SPACEDIM> lower = {{-0.02, -0.02, -0.02}}; // m
    std::array<double, BL_SPACEDIM> upper = {{ 0.02,  0.02,  0.02}}; // m
    
    RealBox domain;
    
    // in helper_functions.h
    init(domain, nr, lower, upper);
    
    
    /*
     * create an Amr object
     */
    ParmParse pp("amr");
    pp.add("max_grid_size", int(maxBoxSize));
    
    Array<int> error_buf(nLevels, 0);
    
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
    ParmParse pgeom("geometry");
    Array<int> is_per = { 0, 0, 0};
    pgeom.addarr("is_periodic", is_per);
    
    Array<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = nr[i];
    
    
    std::vector<int> rr(nLevels);
    Array<int> rrr(nLevels);
    for (int i = 0; i < nLevels; ++i) {
        rr[i] = 2;
        rrr[i] = 2;
    }
    
    AmrOpal myAmrOpal(&domain, nLevels - 1, nCells, 0 /* cartesian */, rr);
    
    // ========================================================================
    // 2. initialize all particles (just single-level)
    // ========================================================================
    
    const Array<BoxArray>& ba = myAmrOpal.boxArray();
    const Array<DistributionMapping>& dmap = myAmrOpal.DistributionMap();
    const Array<Geometry>& geom = myAmrOpal.Geom();
    
    
    amrplayout_t* playout = new amrplayout_t(geom, dmap, ba, rrr);
    
    std::unique_ptr<amrbunch_t> bunch( new amrbunch_t() );
    bunch->initialize(playout);
    bunch->initializeAmr(); // add attributes: level, grid
    
    bunch->setAllowParticlesNearBoundary(true);
    
    // initialize a particle distribution
    double R = 0.005; // radius of sphere [m]
    int nParticles = 1e6;
    initSphere(R, bunch.get(), nParticles);
    
    msg << "Bunch radius: " << R << " m" << endl
        << "#Particles: " << nParticles << endl
        << "Charge per particle: " << bunch->qm[0] << " C" << endl
        << "Total charge: " << nParticles * bunch->qm[0] << " C" << endl
        << "#Cells per dim for bunch: " << 2.0 * R / *(geom[0].CellSize()) << endl;
    
    // redistribute on single-level
    bunch->update();
    
    myAmrOpal.setBunch(bunch.get());
    
    bunch->gatherStatistics();
    
    // ========================================================================
    // 3. multi-level redistribute
    // ========================================================================
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    
    container_t rhs;
    container_t phi;
    container_t efield;
    
    std::string plotsolve = amrex::Concatenate("plt", 0, 4);
    
    
    // map particles
    Vector_t rmin, rmax;
    bounds(bunch->R, rmin, rmax);
    msg << "rmin = " << rmin << endl << "rmax = " << rmax << endl;
    
    Vector_t low(-1.0, -1.0, -1.0);
    Vector_t up = - low;
    
    msg << endl << "Transformed positions" << endl << endl;
    
//     double factor = domainMapping(*bunch, 1.05 * rmin, 1.05 * rmax, low , up);
    
    bunch->update();
    
    Vector_t tmp_min, tmp_max;
    bounds(bunch->R, tmp_min, tmp_max);
    msg << "rmin = " << tmp_min << endl << "rmax = " << tmp_max << endl;
    
    doSolve(myAmrOpal, bunch.get(), rhs, phi, efield, geom, rrr, nLevels, msg);
    
    msg << endl << "Back to normal positions" << endl << endl;
    
//     domainMapping(*bunch, low , up, 1.05 * rmin, 1.05 * rmax);
    
    bunch->update();
    
    bounds(bunch->R, tmp_min, tmp_max);
    
    msg << "rmin = " << tmp_min << endl << "rmax = " << tmp_max << endl;
    
    
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
        if ( rhs[i]->contains_nan(false) )
            msg << "rho: Nan" << endl;
    }
    
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
        if ( efield[i]->contains_nan(false) )
            msg << "efield: Nan" << endl;
    }
    
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
        
        msg << i << ": " << phi[i]->nGrow() << " "
                  << phi[i]->min(0, 1) << " " << phi[i]->max(0, 1) << endl;
        
        if ( phi[i]->contains_nan() )
            msg << "phi: Nan" << endl;
    }
    
    
//     // scale solution
//     double factor2 = factor * factor;
//     for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
//         phi[i]->mult(factor, 0, 1);
//         
//         efield[i]->mult(factor2, 0, 3, 1);
//     }
    
//     for (int lev = phi.size() - 2; lev >= 0; lev--) {
//         std::cout << "average down: " << lev + 1 << " --> " << lev << std::endl;
//         amrex::average_down(*phi[lev+1], *phi[lev], 0, 1, rr[lev]);
//     }
    
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
        msg << "Max. potential level " << i << ": "<< phi[i]->max(0) << endl
            << "Min. potential level " << i << ": " << phi[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(0) << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(0) << endl;
    }
    
//     writeScalarField(phi, *(geom[0].CellSize()), lower[0], "amr-phi_scalar-level-");
//     writeVectorField(efield, *(geom[0].CellSize()), lower[0]);
    
//     writePlotFile(plotsolve, rhs, phi, efield, rr, geom, 0);
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg("Solver");
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);

    std::stringstream call;
    call << "Call: mpirun -np [#procs] " << argv[0]
         << " [#gridpoints x] [#gridpoints y] [#gridpoints z] [#levels] [max. box size]";
    
    if ( argc < 6 ) {
        msg << call.str() << endl;
        return -1;
    }
    
    // number of grid points in each direction
    Vektor<size_t, 3> nr(std::atoi(argv[1]),
                         std::atoi(argv[2]),
                         std::atoi(argv[3]));
    
    
    msg << "Particle test running with" << endl
        << "- grid       = " << nr << endl;
        
    amrex::Initialize(argc,argv, false);
    size_t nLevels = std::atoi(argv[4]) + 1; // i.e. nLevels = 0 --> only single level
    size_t maxBoxSize = std::atoi(argv[5]);
    doAMReX(nr, nLevels, maxBoxSize, msg);
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();

    return 0;
}
