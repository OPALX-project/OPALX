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

#include "../Solver.h"
#include "../MGTSolver.h"
#include "../AmrOpal.h"

#include "../helper_functions.h"

#include "../writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"
#include <random>

#include <getopt.h>

typedef AmrOpal::amrplayout_t amrplayout_t;
typedef AmrOpal::amrbase_t amrbase_t;
typedef AmrOpal::amrbunch_t amrbunch_t;

typedef Vektor<double, BL_SPACEDIM> Vector_t;

struct param_t {
    Vektor<size_t, 3> nr;
    size_t nLevels;
    size_t maxBoxSize;
    double radius;
    double length;
    size_t nParticles;
    double pCharge;
    bool isFixedCharge;
    bool isWriteYt;
    bool isWriteCSV;
    bool isWriteParticles;
    bool isHelp;
    bool useMgtSolver;
};


bool parseProgOptions(int argc, char* argv[], param_t& params, Inform& msg) {
    /* Parsing Command Line Arguments
     * 
     * 26. June 2017
     * https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example
     */
    
    params.isWriteYt = false;
    params.isFixedCharge = false;
    params.isWriteCSV = false;
    params.isWriteParticles = false;
    params.isHelp = false;
    params.useMgtSolver = false;
    
    int c = 0;
    
    int cnt = 0;
    
    int required = 8;
    
    while ( true ) {
        static struct option long_options[] = {
            { "gridx",          required_argument, 0, 'x' },
            { "gridy",          required_argument, 0, 'y' },
            { "gridz",          required_argument, 0, 'z' },
            { "level",          required_argument, 0, 'l' },
            { "maxgrid",        required_argument, 0, 'm' },
            { "radius",         required_argument, 0, 'r' },
            { "boxlength",      required_argument, 0, 'b' },
            { "nparticles",     required_argument, 0, 'n' },
            { "writeYt",        no_argument,       0, 'w' },
            { "help",           no_argument,       0, 'h' },
	    { "pcharge",        required_argument, 0, 'c' },
            { "writeCSV",       no_argument,       0, 'v' },
            { "writeParticles", no_argument,       0, 'p' },
            { "use-mgt-solver", no_argument,       0, 's' },
            { 0,                0,                 0,  0  }
        };
        
        int option_index = 0;
        
        c = getopt_long(argc, argv, "x:y:z:l:m:r:b:n:w", long_options, &option_index);
        
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
            case 'r':
                params.radius = std::atof(optarg); ++cnt; break;
            case 'b':
                params.length = std::atof(optarg); ++cnt; break;
            case 'n':
                params.nParticles = std::atoi(optarg); ++cnt; break;
            case 'c':
                params.pCharge = std::atof(optarg);
                params.isFixedCharge = true;
                break;
            case 'w':
                params.isWriteYt = true;
                break;
            case 'v':
                params.isWriteCSV = true;
                break;
            case 'p':
                params.isWriteParticles = true;
                break;
            case 's':
                params.useMgtSolver = true;
                break;
            case 'h':
                msg << "Usage: " << argv[0]
                    << endl
                    << "--gridx [#gridpoints in x]" << endl
                    << "--gridy [#gridpoints in y]" << endl
                    << "--gridz [#gridpoints in z]" << endl
                    << "--level [#levels]" << endl
                    << "--maxgrid [max. grid]" << endl
                    << "--radius [sphere radius]" << endl
                    << "--boxlength [cube side length]" << endl
                    << "--nparticles [#particles]" << endl
                    << "--pcharge [charge per particle] (optional)" << endl
                    << "--writeYt (optional)" << endl
                    << "--writeCSV (optional)" << endl
                    << "--writeParticles (optional)" << endl
                    << "--use-mgt-solver (optional)" << endl;
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


void writeCSV(const container_t& phi,
              const container_t& efield,
              double lower, double dx)
{
    // Immediate debug output:
    // Output potential and e-field along axis
    std::string outfile = "potential.grid";
    std::ofstream out;
    for (MFIter mfi(*phi[0]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const FArrayBox& lhs = (*phi[0])[mfi];
        
        for (int proc = 0; proc < amrex::ParallelDescriptor::NProcs(); ++proc) {
            if ( proc == amrex::ParallelDescriptor::MyProc() ) {
                
                if ( proc == 0 ) {
                    out.open(outfile, std::ios::out);
                    out << "$x$ [m], $\\Phi$ [V]" << std::endl;
                } else
                    out.open(outfile, std::ios::app);
                
                int j = 0.5 * (bx.hiVect()[1] - bx.loVect()[1]);
                int k = 0.5 * (bx.hiVect()[2] - bx.loVect()[2]);
                
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
                    IntVect ivec(i, j, k);
                    // add one in order to have same convention as PartBunch::computeSelfField()
                    out << lower + i * dx << ", " << lhs(ivec, 0)  << std::endl;
                }
                out.close();
            }
            amrex::ParallelDescriptor::Barrier();
        }
    }
    
    outfile = "efield.grid";
    for (MFIter mfi(*efield[0]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const FArrayBox& lhs = (*efield[0])[mfi];
        
        for (int proc = 0; proc < amrex::ParallelDescriptor::NProcs(); ++proc) {
            if ( proc == amrex::ParallelDescriptor::MyProc() ) {
                
                if ( proc == 0 ) {
                    out.open(outfile, std::ios::out);
                    out << "$x$ [m], $E_x$ [V/m], $E_x$ [V/m], $E_x$ [V/m]" << std::endl;
                } else
                    out.open(outfile, std::ios::app);
                
                int j = 0.5 * (bx.hiVect()[1] - bx.loVect()[1]);
                int k = 0.5 * (bx.hiVect()[2] - bx.loVect()[2]);
                
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
                    IntVect ivec(i, j, k);
                    // add one in order to have same convention as PartBunch::computeSelfField()
                    out << lower + i * dx << ", " << lhs(ivec, 0) << ", "
                        << lhs(ivec, 1) << ", " << lhs(ivec, 2) << std::endl;
                }
                out.close();
            }
            amrex::ParallelDescriptor::Barrier();
        }
    }
    
}


void writeYt(container_t& rho,
             const container_t& phi,
             const container_t& efield,
             const Array<Geometry>& geom,
             const Array<int>& rr,
             const double& scalefactor)
{
    std::string dir = "yt-testUnifSphere";
    
    double time = 0.0;
    
    for (unsigned int i = 0; i < rho.size(); ++i)
        rho[i]->mult(- Physics::epsilon_0 / scalefactor, 0, 1);
    
    writePlotFile(dir, rho, phi, efield, rr, geom, time, scalefactor);
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

void initSphere(double r,
                std::unique_ptr<amrbunch_t>& bunch,
                int nParticles,
                double charge,
                bool isFixedCharge,
                bool isWriteParticles)
{
    
    int nLocParticles = nParticles / Ippl::getNodes();
    
    bunch->create(nLocParticles);
    
    std::mt19937_64 eng[3];
    
    
    for (int i = 0; i < 3; ++i) {
        eng[i].seed(42 + 3 * i);
        eng[i].discard( nLocParticles * Ippl::myNode());
    }
    
    std::uniform_real_distribution<> ph(-1.0, 1.0);
    std::uniform_real_distribution<> th(0.0, 2.0 * Physics::pi);
    std::uniform_real_distribution<> u(0.0, 1.0);
    
    
    std::string outfile = "amr-particles-level-" + std::to_string(0);
    std::ofstream out;
    
    if ( isWriteParticles && Ippl::getNodes() == 1 )
        out.open(outfile, std::ios::out);

    long double qi = 0.0;

    if ( isFixedCharge )
        qi = charge;
    else
        qi = 4.0 * Physics::pi * Physics::epsilon_0 * r * r / double(nParticles);
    
    for (uint i = 0; i < bunch->getLocalNum(); ++i) {
        // 17. Dec. 2016,
        // http://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability
        // http://stackoverflow.com/questions/5408276/sampling-uniformly-distributed-random-points-inside-a-spherical-volume
        double phi = std::acos( ph(eng[0]) );
        double theta = th(eng[1]);
        double radius = r * std::cbrt( u(eng[2]) );
        
        double x = radius * std::cos( theta ) * std::sin( phi );
        double y = radius * std::sin( theta ) * std::sin( phi );
        double z = radius * std::cos( phi );
        
        if ( isWriteParticles && Ippl::getNodes() == 1 )
            out << x << " " << y << " " << z << std::endl;
        
        bunch->R[i] = Vector_t( x, y, z );    // m
        bunch->qm[i] = qi;   // C
    }
    
    if ( isWriteParticles && Ippl::getNodes() == 1 )
        out.close();
}

void doSolve(AmrOpal& myAmrOpal, amrbunch_t* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& efield,
             const Array<int>& rr,
             Inform& msg,
             const double& scale, const param_t& params)
{
    static IpplTimings::TimerRef solvTimer = IpplTimings::getTimer("solve");
    // =======================================================================                                                                                                                                   
    // 4. prepare for multi-level solve                                                                                                                                                                          
    // =======================================================================

    rhs.resize(params.nLevels);
    phi.resize(params.nLevels);
    efield.resize(params.nLevels);
    

    // Define the density on level 0 from all particles at all levels                                                                                                                                            
    int base_level   = 0;
    int finest_level = myAmrOpal.finestLevel();
    
    for (uint lev = 0; lev < params.nLevels; ++lev) {
        initGridData(rhs, phi, efield,
                     myAmrOpal.boxArray()[lev], myAmrOpal.DistributionMap(lev), lev);
    }
    
   container_t partMF(params.nLevels);
   for (uint lev = 0; lev < params.nLevels; lev++) {
        const amrex::BoxArray& ba = myAmrOpal.boxArray()[lev];
        const amrex::DistributionMapping& dmap = myAmrOpal.DistributionMap(lev);
        partMF[lev].reset(new amrex::MultiFab(ba, dmap, 1, 2));
        partMF[lev]->setVal(0.0, 2);
   }
    
    bunch->AssignDensityFort(bunch->qm, partMF, base_level, 1, finest_level);
    
    for (uint lev = 0; lev < params.nLevels; ++lev) {
        amrex::MultiFab::Copy(*rhs[lev], *partMF[lev], 0, 0, 1, 0);
    }
    
    const Array<Geometry>& geom = myAmrOpal.Geom();
    
    for (uint i = 0; i < rhs.size(); ++i)
        if ( rhs[i]->contains_nan() || rhs[i]->contains_inf() )
            throw std::runtime_error("\033[1;31mError: NANs or INFs on charge grid.\033[0m");
    
    
    // Check charge conservation
    double totCharge = totalCharge(rhs, finest_level, geom);
    
    msg << "Total Charge (computed): " << totCharge << " C" << endl
        << "Vacuum permittivity: " << Physics::epsilon_0 << " F/m (= C/(m V)" << endl;
    
    Real vol = (*(geom[0].CellSize()) * *(geom[0].CellSize()) * *(geom[0].CellSize()) );
    msg << "Cell volume: " << *(geom[0].CellSize()) << "^3 = " << vol << " m^3" << endl;
    
    // eps in C / (V * m)
    double constant = -1.0 / Physics::epsilon_0 * scale;  // in [V m / C]
    for (int i = 0; i <=finest_level; ++i) {
        rhs[i]->mult(constant, 0, 1);       // in [V m]
    }
    
    // **************************************************************************                                                                                                                                
    // Compute the total charge of all particles in order to compute the offset                                                                                                                                  
    //     to make the Poisson equations solvable                                                                                                                                                                
    // **************************************************************************                                                                                                                                

    Real offset = 0.;

    // solve
    if ( params.useMgtSolver ) {
        MGTSolver sol;
        IpplTimings::startTimer(solvTimer);
        sol.solve(rhs,            // [V m]
                  phi,            // [V m^3]
                  efield,       // [V m^2]
                  geom);
    } else {
        Solver sol;
        IpplTimings::startTimer(solvTimer);
        sol.solve_for_accel(rhs,            // [V m]
                            phi,            // [V m^3]
                            efield,       // [V m^2]
                            geom,
                            base_level,
                            finest_level,
                            offset);
    }
    
    // undo scale
    for (int i = 0; i <= finest_level; ++i)
        efield[i]->mult(scale, 0, 3);
    
    IpplTimings::stopTimer(solvTimer);
}

void doAMReX(const param_t& params, Inform& msg)
{
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    double halflength = 0.5 * params.length;
    
    std::array<double, BL_SPACEDIM> lower = {{-halflength, -halflength, -halflength}}; // m
    std::array<double, BL_SPACEDIM> upper = {{ halflength,  halflength,  halflength}}; // m
    
    RealBox domain;
    
    // in helper_functions.h
    init(domain, params.nr, lower, upper);
    
    msg << "Domain: " << domain << endl;
    
    /*
     * create an Amr object
     */
    ParmParse pp("amr");
    pp.add("max_grid_size", int(params.maxBoxSize));
    
    Array<int> error_buf(params.nLevels, 0);
    
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
    ParmParse pgeom("geometry");
    Array<int> is_per = { 0, 0, 0};
    pgeom.addarr("is_periodic", is_per);
    
    Array<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = params.nr[i];
    
    
    std::vector<int> rr(params.nLevels);
    Array<int> rrr(params.nLevels);
    for (unsigned int i = 0; i < params.nLevels; ++i) {
        rr[i] = 2;
        rrr[i] = 2;
    }
    
    AmrOpal myAmrOpal(&domain, params.nLevels - 1, nCells, 0 /* cartesian */, rr);
    
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
    initSphere(params.radius,
               bunch,
               params.nParticles,
               params.pCharge,
               params.isFixedCharge,
               params.isWriteParticles);
    
    msg << "Bunch radius: " << params.radius << " m" << endl
        << "#Particles: " << bunch->getTotalNum() << endl
        << "Charge per particle: " << bunch->qm[0] << " C" << endl
        << "Total charge: " << params.nParticles * bunch->qm[0] << " C" << endl
        << "#Cells per dim for bunch: " << 2.0 * params.radius / *(geom[0].CellSize()) << endl;
    
    // map particles
    double scale = 1.0;
    
    scale = domainMapping(*bunch, scale);
    // redistribute on single-level
    bunch->update();
    
    myAmrOpal.setBunch(bunch.get());
    
    bunch->gatherStatistics();
    
    // ========================================================================
    // 3. multi-level redistribute
    // ========================================================================
    
    container_t rhs;
    container_t phi;
    container_t efield;
    
    msg << endl << "Transformed positions" << endl << endl;
    
    msg << "Bunch radius: " << params.radius * scale << " m" << endl
        << "#Particles: " << params.nParticles << endl
        << "Charge per particle: " << bunch->qm[0] << " C" << endl
        << "Total charge: " << params.nParticles * bunch->qm[0] << " C" << endl
        << "#Cells per dim for bunch: " << 2.0 * params.radius * scale / *(geom[0].CellSize()) << endl;
    
    bunch->update();
    
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, scale/*0.0*/ /*time*/);
    
    doSolve(myAmrOpal, bunch.get(), rhs, phi, efield, rrr, msg, scale, params);
    
    msg << endl << "Back to normal positions" << endl << endl;
    
    domainMapping(*bunch, scale, true);
    
    bunch->update();
    
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
        msg << "Max. potential level " << i << ": "<< phi[i]->max(0) << endl
            << "Min. potential level " << i << ": " << phi[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(0) << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(0) << endl;
    }
    
    double fieldenergy = totalFieldEnergy(efield, rrr);
    
    msg << "Total field energy: " << fieldenergy << endl;
    
    if (params.isWriteCSV && Ippl::getNodes() == 1 && myAmrOpal.maxGridSize(0) == (int)params.nr[0] )
        writeCSV(phi, efield, domain.lo(0) / scale, geom[0].CellSize(0) / scale);
    
    if ( params.isWriteYt )
        writeYt(rhs, phi, efield, geom, rrr, scale);
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg("Solver");
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);
    
    param_t params;
    
    amrex::Initialize(argc,argv, false);
    
    try {
        if ( !parseProgOptions(argc, argv, params, msg) && !params.isHelp )
            throw std::runtime_error("\033[1;31mError: Check the program options.\033[0m");
        else if ( params.isHelp )
            return 0;
    
    
        msg << "Particle test running with" << endl
            << "- grid                  = " << params.nr << endl
            << "- max. grid             = " << params.maxBoxSize << endl
            << "- #level                = " << params.nLevels - 1 << endl
            << "- sphere radius [m]     = " << params.radius << endl
            << "- cube side length [m]  = " << params.length << endl
            << "- #particles            = " << params.nParticles << endl;

        if ( params.isFixedCharge )
            msg << "- charge per particle   = " << params.pCharge << endl;
        
        if ( params.useMgtSolver )
            msg << "- MGT solver is used" << endl;
            
        doAMReX(params, msg);
        
    } catch(const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
    }
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();

    return 0;
}
