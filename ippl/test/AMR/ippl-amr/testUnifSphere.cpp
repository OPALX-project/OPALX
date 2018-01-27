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

#ifdef HAVE_AMR_MG_SOLVER
    #include "../trilinos/AmrMultiGrid.h"
#endif


#include <getopt.h>

typedef AmrOpal::amrplayout_t amrplayout_t;
typedef AmrOpal::amrbase_t amrbase_t;
typedef AmrOpal::amrbunch_t amrbunch_t;

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
#ifdef HAVE_AMR_MG_SOLVER
    bool useTrilinos;
    size_t nsweeps;
    AmrMultiGrid::Boundary bcx;
    AmrMultiGrid::Boundary bcy;
    AmrMultiGrid::Boundary bcz;
    AmrMultiGrid::BaseSolver bs;
    AmrMultiGrid::Smoother smoother;
    AmrMultiGrid::Preconditioner prec;
#endif
    AmrOpal::TaggingCriteria criteria;
    double tagfactor;
};


void getBC(AmrMultiGrid::Boundary& boundary, const char* optarg) {
    std::string bc = optarg;
    
    if ( bc == "dirichlet" )
        boundary = AmrMultiGrid::Boundary::DIRICHLET;
    else if ( bc == "open" )
        boundary = AmrMultiGrid::Boundary::OPEN;
    else if ( bc == "periodic" )
        boundary = AmrMultiGrid::Boundary::PERIODIC;
    else
        throw std::runtime_error("Error: Check boundary condition argument");
}


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
    params.criteria = AmrOpal::kChargeDensity;
    params.tagfactor = 1.0e-14; 
    
#ifdef HAVE_AMR_MG_SOLVER
    params.useTrilinos = false;
    params.nsweeps = 12;
    params.bcx = AmrMultiGrid::Boundary::DIRICHLET;
    params.bcy = AmrMultiGrid::Boundary::DIRICHLET;
    params.bcz = AmrMultiGrid::Boundary::DIRICHLET;
    params.bs = AmrMultiGrid::BaseSolver::CG;
    params.smoother = AmrMultiGrid::Smoother::GAUSS_SEIDEL;
    params.prec = AmrMultiGrid::Preconditioner::NONE;
#endif
    
    
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
#ifdef HAVE_AMR_MG_SOLVER
            { "use-trilinos",   no_argument,       0, 'a' },
            { "nsweeps",        required_argument, 0, 'g' },
            { "smoother",       required_argument, 0, 'q' },
            { "prec",           required_argument, 0, 'o' },
            { "bcx",            required_argument, 0, 'i' },
            { "bcy",            required_argument, 0, 'j' },
            { "bcz",            required_argument, 0, 'k' },
            { "basesolver",     required_argument, 0, 'u' },
#endif
            { "tagging",        required_argument, 0, 't' },
            { "tagging-factor", required_argument, 0, 'f' },
            { 0,                0,                 0,  0  }
        };
        
        int option_index = 0;
        
#ifdef HAVE_AMR_MG_SOLVER
        c = getopt_long(argc, argv, "x:y:z:l:m:r:b:n:whcvpst:f:a:g:q:o:u:i:j:k:", long_options, &option_index);
#else
        c = getopt_long(argc, argv, "x:y:z:l:m:r:b:n:whcvpst:f:", long_options, &option_index);
#endif
        
        if ( c == -1 )
            break;
        
        switch ( c ) {
#ifdef HAVE_AMR_MG_SOLVER
            case 'a':
                params.useTrilinos = true; break;
            case 'g':
                params.nsweeps = std::atoi(optarg); break;
            case 'q':
            {
                std::string smoother = optarg;
                if ( smoother == "SGS" )
                    params.smoother = AmrMultiGrid::Smoother::SGS;
                else if ( smoother == "JACOBI" )
                    params.smoother = AmrMultiGrid::Smoother::JACOBI;
                else
                    throw std::runtime_error("Error: Check smoother argument");
                break;
            }
            case 'o':
            {
                std::string prec = optarg;
                if ( prec == "ILUT" )
                    params.prec = AmrMultiGrid::Preconditioner::ILUT;
                else if ( prec == "CHEBYSHEV" )
                    params.prec = AmrMultiGrid::Preconditioner::CHEBYSHEV;
                else
                    throw std::runtime_error("Error: Check preconditioner argument");
                break;
            }
            case 'i':
            {
                getBC(params.bcx, optarg);
                break;
            }
            case 'j':
            {
                getBC(params.bcy, optarg);
                break;
            }
            case 'k':
            {
                getBC(params.bcz, optarg);
                break;
            }
            case 'u':
            {
                std::string bs = optarg;
                
                if ( bs == "bicgstab" )
                    params.bs = AmrMultiGrid::BaseSolver::BICGSTAB;
                else if ( bs == "minres" )
                    params.bs = AmrMultiGrid::BaseSolver::MINRES;
                else if ( bs == "pcpg" )
                    params.bs = AmrMultiGrid::BaseSolver::PCPG;
                else if ( bs == "gmres" )
                    params.bs = AmrMultiGrid::BaseSolver::GMRES;
                else if ( bs == "stochastic_cg" )
                    params.bs = AmrMultiGrid::BaseSolver::STOCHASTIC_CG;
                else if ( bs == "recycling_cg" )
                    params.bs = AmrMultiGrid::BaseSolver::RECYCLING_CG;
                else if ( bs == "recycling_gmres" )
                    params.bs = AmrMultiGrid::BaseSolver::RECYCLING_GMRES;
#ifdef HAVE_AMESOS2_KLU2
                else if ( bs == "klu2" )
                    params.bs = AmrMultiGrid::BaseSolver::KLU2;
#endif
#ifdef HAVE_AMESOS2_SUPERLU
                else if ( bs == "superlu" )
                    params.bs = AmrMultiGrid::BaseSolver::SUPERLU;
#endif
#ifdef HAVE_AMESOS2_UMFPACK
                else if ( bs == "umfpack" )
                    params.bs = AmrMultiGrid::BaseSolver::UMFPACK;
#endif
#ifdef HAVE_AMESOS2_PARDISO_MKL
                else if ( bs == "pardiso_mkl" )
                    params.bs = AmrMultiGrid::BaseSolver::PARDISO_MKL;
#endif
#ifdef HAVE_AMESOS2_MUMPS
                else if ( bs == "mumps" )
                    params.bs = AmrMultiGrid::BaseSolver::MUMPS;
#endif
#ifdef HAVE_AMESOS2_LAPACK
                else if ( bs == "lapack" )
                    params.bs = AmrMultiGrid::BaseSolver::LAPACK;
#endif
                else
                    throw std::runtime_error("Error: Check base solver argument");
                break;
            }
#endif
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
            case 't':
                if ( std::strcmp("efield", optarg) == 0 )
                    params.criteria = AmrOpal::kEfieldStrength;
                else if ( std::strcmp("potential", optarg) == 0 )
                    params.criteria = AmrOpal::kPotentialStrength;
                break;
            case 'f':
                params.tagfactor = std::atof(optarg);
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
                    << "--use-mgt-solver (optional)" << endl
#ifdef HAVE_AMR_MG_SOLVER
                    << "--use-trilinos (optional)" << endl
                    << "--nsweeps (optional, trilinos only, default: 12)" << endl
                    << "--smoother (optional, trilinos only, default: GAUSS_SEIDEL)" << endl
                    << "--prec (optional, trilinos only, default: NONE)" << endl
                    << "--bcx (optional, dirichlet or open, default: dirichlet)" << endl
                    << "--bcy (optional, dirichlet or open, default: dirichlet)" << endl
                    << "--bcz (optional, dirichlet or open, default: dirichlet)" << endl
                    << "--basesolver (optional, trilinos only, default: CG)" << endl
#endif
                    << "--tagging charge (default) / efield / potential (optional)" << endl
                    << "--tagfactor [charge value / 0 ... 1] (optiona)" << endl;
                params.isHelp = true;
                break;
            case '?':
                break;
            
            default:
                break;
            
        }
    }
    
#ifdef HAVE_AMR_MG_SOLVER
    if ( params.useMgtSolver && params.useTrilinos ) {
        params.useMgtSolver = false;
        msg << "Favouring Trilinos over MGT." << endl;
    }
#endif
    
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
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& lhs = (*phi[0])[mfi];
        
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
                    amrex::IntVect ivec(i, j, k);
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
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& lhs = (*efield[0])[mfi];
        
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
                    amrex::IntVect ivec(i, j, k);
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
             const amrex::Array<amrex::Geometry>& geom,
             const amrex::Array<int>& rr,
             const double& scalefactor)
{
    std::string dir = "yt-testUnifSphere";
    
    double time = 0.0;
    
    for (unsigned int i = 0; i < rho.size(); ++i)
        rho[i]->mult(- Physics::epsilon_0 / scalefactor, 0, 1);
    
    writePlotFile(dir, rho, phi, efield, rr, geom, time, scalefactor);
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
             const amrex::Array<int>& rr,
             Inform& msg,
             const double& scale, const param_t& params)
{
    static IpplTimings::TimerRef solvTimer = IpplTimings::getTimer("solve");
    
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
    
    const amrex::Array<amrex::Geometry>& geom = myAmrOpal.Geom();
    
    for (uint i = 0; i < rhs.size(); ++i)
        if ( rhs[i]->contains_nan() || rhs[i]->contains_inf() )
            throw std::runtime_error("\033[1;31mError: NANs or INFs on charge grid.\033[0m");
    
    
    // Check charge conservation
    double totCharge = totalCharge(rhs, finest_level, geom);
    
    msg << "Total Charge (computed): " << totCharge << " C" << endl
        << "Vacuum permittivity: " << Physics::epsilon_0 << " F/m (= C/(m V)" << endl;
    
    amrex::Real vol = (*(geom[0].CellSize()) * *(geom[0].CellSize()) * *(geom[0].CellSize()) );
    msg << "Cell volume: " << *(geom[0].CellSize()) << "^3 = " << vol << " m^3" << endl;
    
    // eps in C / (V * m)
    double constant = -1.0 / Physics::epsilon_0 ; //* scale;  // in [V m / C]
    for (int i = 0; i <= finest_level; ++i) {
        rhs[i]->mult(constant, 0, 1);       // in [V m]
    }
    
    
    // normalize each level
//     double l0norm[finest_level + 1];
    double l0norm = rhs[finest_level]->norm0(0);
    msg << "l0norm = " << l0norm << endl;
    for (int i = 0; i <= finest_level; ++i) {
//         l0norm[i] = rhs[i]->norm0(0);
        rhs[i]->mult(1.0 / l0norm/*[i]*/, 0, 1);
    }
    
    // **************************************************************************                                                                                                                                
    // Compute the total charge of all particles in order to compute the offset                                                                                                                                  
    //     to make the Poisson equations solvable                                                                                                                                                                
    // **************************************************************************                                                                                                                                

    amrex::Real offset = 0.;

    // solve
#ifdef HAVE_AMR_MG_SOLVER
    if ( params.useTrilinos ) {
        AmrMultiGrid sol(params.bcx, params.bcy, params.bcz,
                         AmrMultiGrid::Interpolater::PIECEWISE_CONST,
                         AmrMultiGrid::Interpolater::LAGRANGE, params.bs,
                         params.prec, params.smoother);
        
        sol.setNumberOfSweeps(params.nsweeps);
    
        IpplTimings::startTimer(solvTimer);
        
        sol.solve(rhs,            // [V m]
                  phi,            // [V m^3]
                  efield,       // [V m^2]
                  geom,
                  base_level,
                  finest_level);
    
        IpplTimings::stopTimer(solvTimer);
        
        msg << "#iterations: " << sol.getNumIters() << endl;
        
    } else if ( params.useMgtSolver ) {
#else
    if ( params.useMgtSolver ) {
#endif
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
    
    // undo normalization
    for (int i = 0; i <= finest_level; ++i) {
        phi[i]->mult(scale * l0norm/*[i]*/, 0, 1);
    }
    
    // undo scale
    for (int i = 0; i <= finest_level; ++i)
        efield[i]->mult(scale * scale * l0norm/*[i]*/, 0, 3);
    
    IpplTimings::stopTimer(solvTimer);
}

void doAMReX(const param_t& params, Inform& msg)
{
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
//     double halflength = 0.5 * params.length;
//     
//     std::array<double, AMREX_SPACEDIM> lower = {{-halflength, -halflength, -halflength}}; // m
//     std::array<double, AMREX_SPACEDIM> upper = {{ halflength,  halflength,  halflength}}; // m
//     
//     amrex::RealBox domain;
//     
//     // in helper_functions.h
//     init(domain, params.nr, lower, upper);
//     
//     msg << "Domain: " << domain << endl;
    
    /*
     * create an Amr object
     */
    amrex::ParmParse pp("amr");
    pp.add("max_grid_size", int(params.maxBoxSize));
    
    amrex::Array<int> error_buf(params.nLevels, 0);
    
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
    amrex::ParmParse pgeom("geometry");
    amrex::Array<int> is_per = { 0, 0, 0};
    pgeom.addarr("is_periodic", is_per);
    
    amrex::Array<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = params.nr[i];
    
    
    std::vector<int> rr(params.nLevels);
    amrex::Array<int> rrr(params.nLevels);
    for (unsigned int i = 0; i < params.nLevels; ++i) {
        rr[i] = 2;
        rrr[i] = 2;
    }
    
    amrex::RealBox amr_domain;
    
    std::array<double, AMREX_SPACEDIM> amr_lower = {{-1.02, -1.02, -1.02}}; // m
    std::array<double, AMREX_SPACEDIM> amr_upper = {{ 1.02,  1.02,  1.02}}; // m
    
    init(amr_domain, params.nr, amr_lower, amr_upper);
    
    AmrOpal myAmrOpal(&amr_domain, params.nLevels - 1, nCells, 0 /* cartesian */, rr);
    
    myAmrOpal.setTagging(params.criteria);
    
    if ( params.criteria == AmrOpal::kChargeDensity )
        myAmrOpal.setCharge(params.tagfactor);
    else
        myAmrOpal.setScalingFactor(params.tagfactor);
    
    // ========================================================================
    // 2. initialize all particles (just single-level)
    // ========================================================================
    
    const amrex::Array<amrex::BoxArray>& ba = myAmrOpal.boxArray();
    const amrex::Array<amrex::DistributionMapping>& dmap = myAmrOpal.DistributionMap();
    amrex::Array<amrex::Geometry>& geom = myAmrOpal.Geom();
    
    
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
    
    bunch->update();
    
    msg << "Bunch radius: " << params.radius << " m" << endl
        << "#Particles: " << bunch->getTotalNum() << endl
        << "Charge per particle: " << bunch->qm[0] << " C" << endl
        << "Total charge: " << params.nParticles * bunch->qm[0] << " C" << endl
        << "#Cells per dim for bunch: " << 2.0 * params.radius / *(geom[0].CellSize()) << endl;
    
    // map particles
    double scale = 1.0;
    
    scale = domainMapping(*bunch, scale);
    
    msg << "Scale: " << scale << endl;
    
    
    // redistribute on single-level
    bunch->update();
    
    myAmrOpal.setBunch(bunch.get());
    
    msg << "//\n//  Process Statistic BEFORE regrid\n//" << endl;
    
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
    
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, scale/*0.0*/ /*time*/);
    
    bunch->gatherLevelStatistics();
    
    msg << "//\n//  Process Statistic AFTER regrid\n//" << endl;
    
    bunch->gatherStatistics();
    
    static IpplTimings::TimerRef statisticsTimer = IpplTimings::getTimer("dump-statistics");
    std::string statistics = "particle-statistics-ncores-" + std::to_string(Ippl::getNodes()) + ".dat";
    IpplTimings::startTimer(statisticsTimer);
    bunch->dumpStatistics(statistics);
    IpplTimings::stopTimer(statisticsTimer);
    
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i)
        msg << "#boxes: " <<  myAmrOpal.boxArray(i).size() << endl;
    
    doSolve(myAmrOpal, bunch.get(), rhs, phi, efield, rrr, msg, scale, params);
    
    msg << endl << "Back to normal positions" << endl << endl;
    
    domainMapping(*bunch, scale, true);
    
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
        msg << "Max. potential level " << i << ": "<< phi[i]->max(0) << endl
            << "Min. potential level " << i << ": " << phi[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(0) << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(0) << endl;
    }
    
    double fieldenergy = totalFieldEnergy(efield, rrr);
    
    msg << "Total field energy: " << fieldenergy << endl;
    
    if (params.isWriteCSV && Ippl::getNodes() == 1 && myAmrOpal.maxGridSize(0) == (int)params.nr[0] )
        writeCSV(phi, efield, amr_domain.lo(0) / scale, geom[0].CellSize(0) / scale);
    
    if ( params.isWriteYt ) {
//         double halflength = 0.5 * params.length;
//     
//         std::array<double, AMREX_SPACEDIM> lower = {{-halflength, -halflength, -halflength}}; // m
//         std::array<double, AMREX_SPACEDIM> upper = {{ halflength,  halflength,  halflength}}; // m
//     
//         amrex::RealBox domain;
//         
//         init(domain, params.nr, lower, upper);
//         
//         amrex::RealBox orig = geom[0].ProbDomain();
//         geom[0].ProbDomain(domain);
//         
//         amrex::IntVect low(0, 0, 0);
//         amrex::IntVect high(params.nr[0] - 1, params.nr[1] - 1, params.nr[2] - 1);    
//         Box bx(low, high);
//         geom[0].Domain(bx);
        
        writeYt(rhs, phi, efield, geom, rrr, scale);
    }
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg(argv[0]);
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);
    
    param_t params;
    
    amrex::Initialize(argc,argv, false);
    
    try {
        if ( !parseProgOptions(argc, argv, params, msg) && !params.isHelp )
            throw std::runtime_error("\033[1;31mError: Check the program options.\033[0m");
        else if ( params.isHelp )
            return 0;
        
        
        std::string tagging = "charge";
        if ( params.criteria == AmrOpal::kEfieldStrength )
            tagging = "efield";
        else if ( params.criteria == AmrOpal::kPotentialStrength )
            tagging = "potential";
    
        msg << "Particle test running with" << endl
            << "- grid                  = " << params.nr << endl
            << "- max. grid             = " << params.maxBoxSize << endl
            << "- #level                = " << params.nLevels - 1 << endl
            << "- sphere radius [m]     = " << params.radius << endl
            << "- cube side length [m]  = " << params.length << endl
            << "- #particles            = " << params.nParticles << endl
            << "- tagging               = " << tagging << endl
            << "- tagging factor        = " << params.tagfactor << endl;

        if ( params.isFixedCharge )
            msg << "- charge per particle   = " << params.pCharge << endl;
        
        if ( params.useMgtSolver )
            msg << "- MGT solver is used" << endl;
        
#ifdef HAVE_AMR_MG_SOLVER
        if ( params.useTrilinos ) {
            std::string smoother = "Gauss-Seidel";
            if ( params.smoother == AmrMultiGrid::Smoother::SGS )
                smoother = "symmetric" + smoother;
            else if ( params.smoother == AmrMultiGrid::Smoother::JACOBI )
                smoother = "Jacobi";
            msg << "- Trilinos solver is used with: " << endl
                << "    - nsweeps:     " << params.nsweeps << endl
                << "    - smoother:    " << smoother
                << endl;
        }
#endif
        doAMReX(params, msg);
        
    } catch(const std::exception& ex) {
        msg << ex.what() << endl;
    }
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();
    
    std::string fn = std::string(argv[0]) + "-timing.dat";
    
    std::map<std::string, unsigned int> problemSize;
    
    problemSize["level"] = params.nLevels;
    problemSize["maxgrid"] = params.maxBoxSize;
    problemSize["gridx"] = params.nr[0];
    problemSize["gridy"] = params.nr[1];
    problemSize["gridz"] = params.nr[2];
    
    IpplTimings::print(fn, problemSize);

    return 0;
}
