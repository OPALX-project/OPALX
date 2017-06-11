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

void initSphere(double r, std::unique_ptr<amrbunch_t>& bunch, int nParticles) {
    
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
    
    if ( Ippl::getNodes() == 1 )
        out.open(outfile, std::ios::out);
    
    long double qi = 4.0 * Physics::pi * Physics::epsilon_0 * r * r / double(nParticles);
    
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
        
        if ( Ippl::getNodes() == 1 )
            out << x << " " << y << " " << z << std::endl;
        
        bunch->R[i] = Vector_t( x, y, z );    // m
        bunch->qm[i] = qi;   // C
    }
    
    if (Ippl::getNodes() == 1 )
        out.close();
}

void doSolve(AmrOpal& myAmrOpal, amrbunch_t* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& efield,
             const Array<int>& rr, int nLevels,
             Inform& msg,
             const double& scale)
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
    
    const Array<Geometry>& geom = myAmrOpal.Geom();
    
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
    
//     for (int lev = base_level; lev <= finest_level; lev++) {
//         phi[lev]->setBndry(0.0);
// //         rhs[lev]->FillBoundary();
//         phi[lev]->FillBoundary();
//     }

    // solve                                                                                                                                                                                                     
    Solver sol;
    IpplTimings::startTimer(solvTimer);
    sol.solve_for_accel(rhs,            // [V m]
                        phi,            // [V m^3]
                        efield,       // [V m^2]
                        geom,
                        base_level,
                        finest_level,
                        offset);
    
//     for (int i = 0; i <=finest_level; ++i) {
//         rhs[i]->mult(1.0 / constant, 0, 1);
//     }
    
    // undo scale
    for (int i = 0; i <= finest_level; ++i)
        efield[i]->mult(scale, 0, 3);
    
    IpplTimings::stopTimer(solvTimer);
}

void doAMReX(const Vektor<size_t, 3>& nr,
             int nLevels,
             size_t maxBoxSize,
             double radius,
             double length,
             Inform& msg)
{
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    double halflength = 0.5 * length;
    
    std::array<double, BL_SPACEDIM> lower = {{-halflength, -halflength, -halflength}}; // m
    std::array<double, BL_SPACEDIM> upper = {{ halflength,  halflength,  halflength}}; // m
    
    RealBox domain;
    
    // in helper_functions.h
    init(domain, nr, lower, upper);
    
    msg << "Domain: " << domain << endl;
    
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
    int nParticles = 1e6;
    initSphere(radius, bunch, nParticles);
    
    msg << "Bunch radius: " << radius << " m" << endl
        << "#Particles: " << bunch->getTotalNum() << endl
        << "Charge per particle: " << bunch->qm[0] << " C" << endl
        << "Total charge: " << nParticles * bunch->qm[0] << " C" << endl
        << "#Cells per dim for bunch: " << 2.0 * radius / *(geom[0].CellSize()) << endl;
    
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
    
    std::string plotsolve = amrex::Concatenate("plt", 0, 4);
    
    
    
    msg << endl << "Transformed positions" << endl << endl;
    
    msg << "Bunch radius: " << radius * scale << " m" << endl
        << "#Particles: " << nParticles << endl
        << "Charge per particle: " << bunch->qm[0] << " C" << endl
        << "Total charge: " << nParticles * bunch->qm[0] << " C" << endl
        << "#Cells per dim for bunch: " << 2.0 * radius * scale / *(geom[0].CellSize()) << endl;
    
    bunch->update();
    
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    
    doSolve(myAmrOpal, bunch.get(), rhs, phi, efield, rrr, nLevels, msg, scale);
    
    msg << endl << "Back to normal positions" << endl << endl;
    
    domainMapping(*bunch, scale, true);
    
    bunch->update();
    
    
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
    
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
        msg << "Max. potential level " << i << ": "<< phi[i]->max(0) << endl
            << "Min. potential level " << i << ": " << phi[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(0) << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(0) << endl;
    }
    
//     writeScalarField(phi, *(geom[0].CellSize()), lower[0], "amr-phi_scalar-level-");
//     writeVectorField(efield, *(geom[0].CellSize()), lower[0]);
    
//     writePlotFile(plotsolve, rhs, phi, efield, rr, geom, 0);

    if (Ippl::getNodes() == 1 && myAmrOpal.maxGridSize(0) == (int)nr[0] )
        writeCSV(phi, efield, domain.lo(0) / scale, geom[0].CellSize(0) / scale);
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg("Solver");
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);

    std::stringstream call;
    call << "Call: mpirun -np [#procs] " << argv[0]
         << " [#gridpoints x] [#gridpoints y] [#gridpoints z] [#levels] [max. box size] [radius] [side length]";
    
    if ( argc < 8 ) {
        msg << call.str() << endl;
        return -1;
    }
    
    // number of grid points in each direction
    Vektor<size_t, 3> nr(std::atoi(argv[1]),
                         std::atoi(argv[2]),
                         std::atoi(argv[3]));
    
    
        
    amrex::Initialize(argc,argv, false);
    size_t nLevels = std::atoi(argv[4]) + 1; // i.e. nLevels = 0 --> only single level
    size_t maxBoxSize = std::atoi(argv[5]);
    double radius = std::atof(argv[6]);
    double length = std::atof(argv[7]);
    
    msg << "Particle test running with" << endl
        << "- grid                  = " << nr << endl
        << "- max. grid             = " << maxBoxSize << endl
        << "- #level                = " << nLevels - 1 << endl
        << "- sphere radius [m]     = " << radius << endl
        << "- cube side length [m]  = " << length << endl;
    
        
    doAMReX(nr, nLevels, maxBoxSize, radius, length, msg);
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();

    return 0;
}
