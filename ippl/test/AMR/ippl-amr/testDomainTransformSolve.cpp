/*!
 * @file testDomainTransformSolve.cpp
 * @author Matthias Frey
 * @date May 2017
 * @details Compute \f$\Delta\phi = -\rho\f$ where the charge
 * distribution is a uniformly charged sphere.
 * 
 * You can perform either a simulation using scaled or no-scaled
 * particle coordinates. One should get the same result.
 * It additionally writes the electric field at the particle location
 * (original locations, i.e. without scaling) to a file in order to
 * compare the two approaches. The domain in both cases is chosen such that
 * the distance of the boundary of the box to the particles is the same.
 * 
 * BC:      Dirichlet boundary conditions\n
 * #Particles 1000
 * Sphere radius: 0.005 [m]
 * 
 * Call:\n
 *  mpirun -np [#cores] testNewTracker [#gridpoints x] [#gridpoints y] [#gridpoints z]
 *                                     [#particles] [#levels] [max. box size] [w or w/o scaling]
 * 
 * @brief Computes \f$\Delta\phi = -\rho\f$
 */

#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <sstream>
#include <iomanip>
#include <fstream>


#include <AMReX_ParmParse.H>


#include "../Distribution.h"
#include "../Solver.h"
#include "../AmrOpal.h"

#include "../helper_functions.h"

#include <cmath>

#include "Physics/Physics.h"

typedef AmrOpal::amrplayout_t amrplayout_t;
typedef AmrOpal::amrbase_t amrbase_t;
typedef AmrOpal::amrbunch_t amrbunch_t;

typedef Vektor<double, BL_SPACEDIM> Vector_t;
typedef std::array<double, BL_SPACEDIM> bc_t;


void initSphere(double r, amrbunch_t* bunch, int nParticles) {
    bunch->create(nParticles / Ippl::getNodes());
    
    std::mt19937_64 eng;
    
    if ( Ippl::myNode() )
        eng.seed(42 + Ippl::myNode() );
    
    std::uniform_real_distribution<> ph(-1.0, 1.0);
    std::uniform_real_distribution<> th(0.0, 2.0 * Physics::pi);
    std::uniform_real_distribution<> u(0.0, 1.0);
    
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
        
        bunch->R[i] = Vector_t( x, y, z );    // m
        bunch->qm[i] = qi;   // C
    }
}


void setup(AmrOpal* &myAmrOpal, std::unique_ptr<amrbunch_t>& bunch,
           const bc_t& lower, const bc_t& upper,
           const Vektor<size_t, 3>& nr, size_t nParticles,
           int nLevels, size_t maxBoxSize, Inform& msg)
{
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
    
    myAmrOpal = new AmrOpal(&domain, nLevels - 1, nCells, 0 /* cartesian */);
    
    
    // ========================================================================
    // 2. initialize all particles (just single-level)
    // ========================================================================
    
    const Array<BoxArray>& ba = myAmrOpal->boxArray();
    const Array<DistributionMapping>& dmap = myAmrOpal->DistributionMap();
    const Array<Geometry>& geom = myAmrOpal->Geom();
    
    Array<int> rr(nLevels);
    for (int i = 0; i < nLevels; ++i)
        rr[i] = 2;
    
    
    amrplayout_t* playout = new amrplayout_t(geom, dmap, ba, rr);
    
    bunch.reset( new amrbunch_t() );
    bunch->initialize(playout);
    bunch->initializeAmr(); // add attributes: level, grid
    
    bunch->setAllowParticlesNearBoundary(true);
    
    // initialize a particle distribution
    double R = 0.005; // radius of sphere [m]
    initSphere(R, bunch.get(), nParticles);
    
    msg << "Bunch radius: " << R << " m" << endl
        << "#Particles: " << nParticles << endl
        << "Charge per particle: " << bunch->qm[0] << " C" << endl
        << "Total charge: " << nParticles * bunch->qm[0] << " C" << endl
        << "#Cells per dim for bunch: " << 2.0 * R / *(geom[0].CellSize()) << endl;
    
    // redistribute on single-level
    bunch->update();
    
    msg << "Single-level statistics" << endl;
    bunch->gatherStatistics();
    
    msg << "#Particles: " << nParticles << endl
        << "Charge per particle: " << bunch->qm[0] << " C" << endl
        << "Total charge: " << nParticles * bunch->qm[0] << " C" << endl;
    
    
    myAmrOpal->setBunch(bunch.get());
    
    Vector_t rmin , rmax;
    bounds(bunch->R, rmin, rmax);
    
    std::cout << rmin << " " << rmax << std::endl;
}


double domainMapping(amrbase_t& PData, const double& scale, bool inverse = false)
{
    Vector_t rmin, rmax;
    bounds(PData.R, rmin, rmax);
    
    double absmax = scale;
    
    if ( !inverse ) {
        Vector_t tmp = Vector_t(std::max( std::abs(rmin[0]), std::abs(rmax[0]) ),
                                std::max( std::abs(rmin[1]), std::abs(rmax[1]) ),
                                std::max( std::abs(rmin[2]), std::abs(rmax[2]) )
                               );
        
        absmax = std::max( tmp[0], tmp[1] );
        absmax = std::max( absmax, tmp[2] );
    }
    
    Vector_t vscale = Vector_t(absmax, absmax, absmax);
    
    for (unsigned int i = 0; i < PData.getLocalNum(); ++i) {
        PData.R[i] /= vscale;
    }
    
    return 1.0 / absmax;
}


void doSolve(AmrOpal* myAmrOpal, amrbunch_t* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& efield,
             int nLevels,
             Inform& msg,
             const double& scale)
{
    // =======================================================================                                                                                                                                   
    // 4. prepare for multi-level solve                                                                                                                                                                          
    // =======================================================================
    
    rhs.resize(nLevels);
    phi.resize(nLevels);
    efield.resize(nLevels);
    
    for (int lev = 0; lev < nLevels; ++lev) {
        initGridData(rhs, phi, efield,
                     myAmrOpal->boxArray()[lev], myAmrOpal->DistributionMap(lev), lev);
    }

    // Define the density on level 0 from all particles at all levels                                                                                                                                            
    int base_level   = 0;
    int finest_level = myAmrOpal->finestLevel();
    
   container_t partMF(nLevels);
   for (int lev = 0; lev < nLevels; lev++) {
        const amrex::BoxArray& ba = myAmrOpal->boxArray()[lev];
        const amrex::DistributionMapping& dmap = myAmrOpal->DistributionMap(lev);
        partMF[lev].reset(new amrex::MultiFab(ba, dmap, 1, 2));
        partMF[lev]->setVal(0.0, 2);
   }
    
    bunch->AssignDensityFort(bunch->qm, partMF, base_level, 1, finest_level);
    
    for (int lev = 0; lev < nLevels; ++lev) {
        amrex::MultiFab::Copy(*rhs[lev], *partMF[lev], 0, 0, 1, 0);
    }
    
    // eps in C / (V * m)
    double constant = -1.0 / Physics::epsilon_0 * scale;  // in [V m / C]
    for (int i = 0; i <=finest_level; ++i) {
#ifdef UNIQUE_PTR
        rhs[i]->mult(constant, 0, 1);       // in [V m]
#else
        rhs[i].mult(constant, 0, 1);
#endif
    }
    
    // **************************************************************************                                                                                                                                
    // Compute the total charge of all particles in order to compute the offset                                                                                                                                  
    //     to make the Poisson equations solvable                                                                                                                                                                
    // **************************************************************************                                                                                                                                

    Real offset = 0.;
    
    // solve                                                                                                                                                                                                     
    Solver sol;
    sol.solve_for_accel(rhs,
                        phi,
                        efield,
                        myAmrOpal->Geom(),
                        base_level,
                        finest_level,
                        offset,
                        false);
    
    // for plotting unnormalize
    for (int i = 0; i <=finest_level; ++i) {
#ifdef UNIQUE_PTR
        rhs[i]->mult(1.0 / constant, 0, 1);       // in [V m]
#else
        rhs[i].mult(1.0 / constant, 0, 1);
#endif
    }
    
    
    bunch->InterpolateFort(bunch->E, efield, base_level, finest_level);
    
    Vector_t vscale = Vector_t(scale, scale, scale);
    
    bunch->E *= vscale;
}

void doWithScaling(const Vektor<size_t, 3>& nr, size_t nParticles,
                   int nLevels, size_t maxBoxSize, Inform& msg)
{
    bc_t lower = {{-1.025, -1.025, -1.025}}; // m
    bc_t upper = {{ 1.025,  1.025,  1.025}}; // m
    
    
    AmrOpal* myAmrOpal = 0;
    std::unique_ptr<amrbunch_t> bunch;
    
    setup(myAmrOpal, bunch, lower, upper, nr,
          nParticles, nLevels, maxBoxSize, msg);
    
    
    container_t rhs;
    container_t phi;
    container_t efield;
    
//     bunch->python_format(0);
    
    // map particles
    double scale = 1.0;
    
    scale = domainMapping(*bunch, scale);
    
    msg << endl << "Transformed positions" << endl << endl;
    
    bunch->update();
    
    for (int i = 0; i <= myAmrOpal->finestLevel() && i < myAmrOpal->maxLevel(); ++i)
        myAmrOpal->regrid(i /*lbase*/, 0.0 /*time*/);
    
    msg << "Multi-level statistics" << endl;
    bunch->gatherStatistics();
    doSolve(myAmrOpal, bunch.get(), rhs, phi, efield, nLevels, msg, scale);
    
    msg << endl << "Back to normal positions" << endl << endl;
    
    domainMapping(*bunch, scale, true);
    
    bunch->update();
    
    for (int i = 0; i <= myAmrOpal->finestLevel(); ++i) {
        if ( efield[i]->contains_nan(false) )
            msg << "Efield: Nan" << endl;
    }
    
    for (int i = 0; i <= myAmrOpal->finestLevel(); ++i) {
        msg << "Max. potential level " << i << ": "<< phi[i]->max(0) << endl
            << "Min. potential level " << i << ": " << phi[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(0) * scale << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(0) * scale << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(1) * scale << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(1) * scale << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(2) * scale << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(2) * scale << endl;
    }
    
    
    std::ofstream out;
    for (int i = 0; i < Ippl::getNodes(); ++i) {
        
        if ( i == Ippl::myNode() ) {
            
            if ( i == 0 )
                out.open("scaling.dat", std::ios::out);
            else
                out.open("scaling.dat", std::ios::app);
            
            for (std::size_t i = 0; i < bunch->getLocalNum(); ++i)
                out << bunch->E[i](0) << " "
                    << bunch->E[i](1) << " "
                    << bunch->E[i](2) << std::endl;
            out.close();
        }
        Ippl::Comm->barrier();
    }
    
    delete myAmrOpal;
}

void doWithoutScaling(const Vektor<size_t, 3>& nr, size_t nParticles,
                      int nLevels, size_t maxBoxSize, Inform& msg)
{
    double max = 1.025 * 0.004843681885;
    bc_t lower = {{-max, -max, -max}}; // m
    bc_t upper = {{ max,  max,  max}}; // m
    
    
    AmrOpal* myAmrOpal = 0;
    std::unique_ptr<amrbunch_t> bunch;
    
    setup(myAmrOpal, bunch, lower, upper, nr,
          nParticles, nLevels, maxBoxSize, msg);
    
    container_t rhs;
    container_t phi;
    container_t efield;
    
//     bunch->python_format(0);
    
    // map particles
    double scale = 1.0;
    
    for (int i = 0; i <= myAmrOpal->finestLevel() && i < myAmrOpal->maxLevel(); ++i)
        myAmrOpal->regrid(i /*lbase*/, 0.0 /*time*/);
    
    bunch->update();
    
    msg << "Multi-level statistics" << endl;
    bunch->gatherStatistics();
    
    
    doSolve(myAmrOpal, bunch.get(), rhs, phi, efield, nLevels, msg, scale);
    
    for (int i = 0; i <= myAmrOpal->finestLevel(); ++i) {
        if ( efield[i]->contains_nan(false) )
            msg << "Efield: Nan" << endl;
    }
    
    for (int i = 0; i <= myAmrOpal->finestLevel(); ++i) {
        msg << "Max. potential level " << i << ": "<< phi[i]->max(0) << endl
            << "Min. potential level " << i << ": " << phi[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(0) << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(1) << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(1) << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(2) << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(2) << endl;
    }
    
    std::ofstream out;
    for (int i = 0; i < Ippl::getNodes(); ++i) {
        
        if ( i == Ippl::myNode() ) {
            
            if ( i == 0 )
                out.open("no_scaling.dat", std::ios::out);
            else
                out.open("no_scaling.dat", std::ios::app);
            
            for (std::size_t i = 0; i < bunch->getLocalNum(); ++i)
                out << bunch->E[i](0) << " "
                    << bunch->E[i](1) << " "
                    << bunch->E[i](2) << std::endl;
            out.close();
        }
        Ippl::Comm->barrier();
    }
    
    delete myAmrOpal;
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg("Solver");

    std::stringstream call;
    call << "Call: mpirun -np [#procs] " << argv[0]
         << " [#gridpoints x] [#gridpoints y] [#gridpoints z] "
         << "[#levels] [max. box size] [w or w/o scaling]";
    
    if ( argc < 7 ) {
        msg << call.str() << endl;
        return -1;
    }
    
    // number of grid points in each direction
    Vektor<size_t, 3> nr(std::atoi(argv[1]),
                         std::atoi(argv[2]),
                         std::atoi(argv[3]));
    
    
    size_t nParticles = 1e3;
    
    msg << "Particle test running with" << endl
        << "- #particles = " << nParticles << endl
        << "- grid       = " << nr << endl;
        
    amrex::Initialize(argc,argv, false);
    size_t nLevels = std::atoi(argv[4]) + 1; // i.e. nLevels = 0 --> only single level
    size_t maxBoxSize = std::atoi(argv[5]);
    
    if ( std::atoi( argv[6] ) ) {
        msg << "With Scaling" << endl;
        doWithScaling(nr, nParticles, nLevels, maxBoxSize, msg);
    } else {
        msg << "Without Scaling" << endl;
        doWithoutScaling(nr, nParticles, nLevels, maxBoxSize, msg);
    }
    
    return 0;
}
