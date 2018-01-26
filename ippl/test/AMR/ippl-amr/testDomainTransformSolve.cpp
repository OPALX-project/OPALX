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
 * #Particles 1e6
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
#include <AMReX_PlotFileUtil.H>


#include "../Distribution.h"
#include "../Solver.h"
#include "../AmrOpal.h"

#include "../helper_functions.h"

#include "../writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"

#include <getopt.h>

using namespace amrex;


typedef AmrOpal::amrplayout_t amrplayout_t;
typedef AmrOpal::amrbase_t amrbase_t;
typedef AmrOpal::amrbunch_t amrbunch_t;

typedef Vektor<double, AMREX_SPACEDIM> Vector_t;
typedef std::array<double, AMREX_SPACEDIM> bc_t;

// #include "../AmrWriter.h"

struct param_t {
    Vektor<size_t, 3> nr;
    size_t nLevels;
    size_t maxBoxSize;
    bool scaling;
    bool isHelp;
    bool isWriteYt;
};


bool parseProgOptions(int argc, char* argv[], param_t& params, Inform& msg) {
    /* Parsing Command Line Arguments
     * 
     * 26. June 2017
     * https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example
     */
    
    params.isHelp = false;
    
    int c = 0;
    
    int cnt = 0;
    
    int required = 5;
    
    params.scaling = false;
    params.isWriteYt = false;
    
    while ( true ) {
        static struct option long_options[] = {
            { "gridx",          required_argument, 0, 'x' },
            { "gridy",          required_argument, 0, 'y' },
            { "gridz",          required_argument, 0, 'z' },
            { "level",          required_argument, 0, 'l' },
            { "maxgrid",        required_argument, 0, 'm' },
            { "help",           no_argument,       0, 'h' },
            { "scaling",        no_argument,       0, 's' },
            { "writeYt",        no_argument,       0, 'w' },
            { 0,                0,                 0,  0  }
        };
        
        int option_index = 0;
        
        c = getopt_long(argc, argv, "x:y:z:l:m:hs", long_options, &option_index);
        
        if ( c == -1 )
            break;
        
        switch ( c ) {
            case 'x':
                params.nr[0] = std::atoi(optarg); ++cnt; break;
            case 'y':
                params.nr[1] = std::atoi(optarg); ++cnt; break;
            case 'z':
                params.nr[2] = std::atoi(optarg); ++cnt; break;
            case 'l':
                params.nLevels = std::atoi(optarg) + 1; ++cnt; break;
            case 'm':
                params.maxBoxSize = std::atoi(optarg); ++cnt; break;
            case 's':
                params.scaling = true;
                break;
            case 'w':
                params.isWriteYt = true;
                break;
            case 'h':
                msg << "Usage: " << argv[0]
                    << endl
                    << "--gridx [#gridpoints in x]" << endl
                    << "--gridy [#gridpoints in y]" << endl
                    << "--gridz [#gridpoints in z]" << endl
                    << "--level [#levels]" << endl
                    << "--maxgrid [max. grid]" << endl
                    << "--scaling (optional)" << endl
                    << "--writeYt (optional)" << endl;
                params.isHelp = true;
                break;
            case '?':
                break;
            
            default:
                break;
            
        }
    }
    
    return ( cnt == required );
}


void writeYt(container_t& rho,
             const container_t& phi,
             const container_t& efield,
             const Array<Geometry>& geom,
             const Array<int>& rr,
             const double& scalefactor,
             std::string dir)
{
    double time = 0.0;
    
    for (unsigned int i = 0; i < rho.size(); ++i)
        rho[i]->mult(- Physics::epsilon_0 / scalefactor, 0, 1);
    
    writePlotFile(dir, rho, phi, efield, rr, geom, time, scalefactor);
}


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
    
    msg << "Bounds: rmin = " << rmin << " rmax = " << rmax << endl;
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
    
    // undo scale
    for (int i = 0; i <= finest_level; ++i)
        efield[i]->mult(scale, 0, 3);
    
    bunch->InterpolateFort(bunch->E, efield, base_level, finest_level);
    
//     Vector_t vscale = Vector_t(scale, scale, scale);
    
//     bunch->E *= vscale;
}

void doWithScaling(const param_t& params, size_t nParticles, Inform& msg)
{
    bc_t lower = {{-1.025, -1.025, -1.025}}; // m
    bc_t upper = {{ 1.025,  1.025,  1.025}}; // m
    
    
    AmrOpal* myAmrOpal = 0;
    std::unique_ptr<amrbunch_t> bunch;
    
    setup(myAmrOpal, bunch, lower, upper, params.nr,
          nParticles, params.nLevels, params.maxBoxSize, msg);
    
    
    container_t rhs;
    container_t phi;
    container_t efield;
    
//     bunch->python_format(0);
    
    // map particles
    double scale = 1.0;
    
    scale = domainMapping(*bunch, scale);
    
    msg << "Scale: " << scale << endl;
    
    msg << endl << "Transformed positions" << endl << endl;
    
    bunch->update();
    
    for (int i = 0; i <= myAmrOpal->finestLevel() && i < myAmrOpal->maxLevel(); ++i)
        myAmrOpal->regrid(i /*lbase*/, 0.0 /*time*/);
    
    msg << "Multi-level statistics" << endl;
    bunch->gatherStatistics();
    doSolve(myAmrOpal, bunch.get(), rhs, phi, efield, params.nLevels, msg, scale);
    
    msg << endl << "Back to normal positions" << endl << endl;
    
    domainMapping(*bunch, scale, true);
    
    bunch->update();
    
    for (int i = 0; i <= myAmrOpal->finestLevel(); ++i) {
        if ( efield[i]->contains_nan(false) )
            msg << "Efield: Nan" << endl;
    }
    
    
    
//     AmrWriter writer;
//     std::string dir = "grid-data";
//     if ( !writer.save(dir, amrex::GetArrOfPtrs(phi), myAmrOpal->Geom(), rr, 0.0) )
//         msg << "Couldn't write potential" << endl;
    
    for (int i = 0; i <= myAmrOpal->finestLevel(); ++i) {
        msg << "Max. potential level " << i << ": "<< phi[i]->max(0) << endl
            << "Min. potential level " << i << ": " << phi[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(0) /* * scale */ << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(0) /* * scale */ << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(1) /* * scale */ << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(1) /* * scale */ << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(2) /* * scale */ << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(2) /* * scale */ << endl;
    }
    
    
//     std::ofstream out;
//     for (int i = 0; i < Ippl::getNodes(); ++i) {
//         
//         if ( i == Ippl::myNode() ) {
//             
//             if ( i == 0 )
//                 out.open("scaling.dat", std::ios::out);
//             else
//                 out.open("scaling.dat", std::ios::app);
//             
//             for (std::size_t i = 0; i < bunch->getLocalNum(); ++i)
//                 out << bunch->E[i](0) << " "
//                     << bunch->E[i](1) << " "
//                     << bunch->E[i](2) << std::endl;
//             out.close();
//         }
//         Ippl::Comm->barrier();
//     }
    
    Array<int> rr(params.nLevels);
    for (size_t i = 0; i < params.nLevels; ++i)
        rr[i] = 2;
    
    if ( params.isWriteYt ) {
        const Array<Geometry>& geom = myAmrOpal->Geom();
        writeYt(rhs, phi, efield, geom, rr, scale, "yt-scaling");
    }
    
    delete myAmrOpal;
}

void doWithoutScaling(const param_t& params, size_t nParticles, Inform& msg)
{
//     double max = 1.025 * 0.004843681885; // 1e3 particles
//     double max = 1.025 * 0.004996545409; // 1e6 particles
    double max = 1.025 * 0.004998988423;   // 1e7 particles
    bc_t lower = {{-max, -max, -max}}; // m
    bc_t upper = {{ max,  max,  max}}; // m
    
    AmrOpal* myAmrOpal = 0;
    std::unique_ptr<amrbunch_t> bunch;
    
    setup(myAmrOpal, bunch, lower, upper, params.nr,
          nParticles, params.nLevels, params.maxBoxSize, msg);
    
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
    
    
    doSolve(myAmrOpal, bunch.get(), rhs, phi, efield, params.nLevels, msg, scale);
    
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
    
//     std::ofstream out;
//     for (int i = 0; i < Ippl::getNodes(); ++i) {
//         
//         if ( i == Ippl::myNode() ) {
//             
//             if ( i == 0 )
//                 out.open("no_scaling.dat", std::ios::out);
//             else
//                 out.open("no_scaling.dat", std::ios::app);
//             
//             for (std::size_t i = 0; i < bunch->getLocalNum(); ++i)
//                 out << bunch->E[i](0) << " "
//                     << bunch->E[i](1) << " "
//                     << bunch->E[i](2) << std::endl;
//             out.close();
//         }
//         Ippl::Comm->barrier();
//     }
    
    
    Array<int> rr(params.nLevels);
    for (size_t i = 0; i < params.nLevels; ++i)
        rr[i] = 2;
    
    if ( params.isWriteYt ) {
        const Array<Geometry>& geom = myAmrOpal->Geom();
        writeYt(rhs, phi, efield, geom, rr, scale, "yt-no-scaling");
    }
    
    delete myAmrOpal;
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg("Solver");
    
    param_t params;

    amrex::Initialize(argc,argv, false);
    
    try {
        if ( !parseProgOptions(argc, argv, params, msg) && !params.isHelp )
            throw std::runtime_error("\033[1;31mError: Check the program options.\033[0m");
        else if ( params.isHelp )
            return 0;
        
        
        size_t nParticles = 1e7;
        
        std::string scale = "no";
        
        if ( params.scaling )
            scale = "yes";
    
        msg << "Particle test running with" << endl
            << "- #particles = " << nParticles << endl
            << "- grid       = " << params.nr << endl
            << "- scaling    = " << scale << endl;
        
        if ( params.scaling )
            doWithScaling(params, nParticles, msg);
        else
            doWithoutScaling(params, nParticles, msg);
        
    } catch(const std::exception& ex) {
        msg << ex.what() << endl;
    }
    
    return 0;
}
