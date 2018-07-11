#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <sstream>
#include <iomanip>

#include <tuple>

#include <boost/filesystem.hpp>


#include <AMReX_ParmParse.H>

#include "../Distribution.h"
#include "../Solver.h"
#include "../MGTSolver.h"
#include "../AmrOpal.h"

#include "../helper_functions.h"


#ifdef HAVE_AMR_MG_SOLVER
    #include "../trilinos/AmrMultiGrid.h"
#endif


#include <getopt.h>

#include <cmath>

#include "Physics/Physics.h"

typedef AmrOpal::amrplayout_t amrplayout_t;
typedef AmrOpal::amrbase_t amrbase_t;
typedef AmrOpal::amrbunch_t amrbunch_t;

typedef Vektor<double, AMREX_SPACEDIM> Vector_t;

struct param_t {
    Vektor<size_t, 3> nr;
    size_t nLevels;
    size_t maxBoxSize;
    size_t blocking_factor;
    double timestep;
    size_t nIterations;
    Distribution::Type type;
    bool isHelp;
    bool useMgtSolver;
#ifdef HAVE_AMR_MG_SOLVER
    bool useTrilinos;
    size_t nsweeps;
    std::string bs;
    std::string smoother;
    std::string prec;
    bool rebalance;
    double bboxincr;
#endif
    AmrOpal::TaggingCriteria criteria;
    double tagfactor;
};


bool parseProgOptions(int argc, char* argv[], param_t& params, Inform& msg) {
    /* Parsing Command Line Arguments
     * 
     * 26. June 2017
     * https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example
     */
    params.isHelp = false;
    params.useMgtSolver = false;
    params.criteria = AmrOpal::kChargeDensity;
    params.tagfactor = 1.0e-14; 
    
#ifdef HAVE_AMR_MG_SOLVER
    params.useTrilinos = false;
    params.nsweeps = 12;
    params.bs = "CG";
    params.smoother = "GS";
    params.prec = "NONE";
    params.rebalance = false;
#endif
    
    
    int c = 0;
    
    int cnt = 0;
    
    int required = 10;
    
    while ( true ) {
        static struct option long_options[] = {
            { "gridx",           required_argument, 0, 'x' },
            { "gridy",           required_argument, 0, 'y' },
            { "gridz",           required_argument, 0, 'z' },
            { "level",           required_argument, 0, 'l' },
            { "maxgrid",         required_argument, 0, 'm' },
            { "blocking_factor", required_argument, 0, 'd' },
            { "help",            no_argument,       0, 'h' },
            { "use-mgt-solver",  no_argument,       0, 's' },
            { "timestep",        required_argument, 0, 'S' },
            { "iterations",      required_argument, 0, 'I' },
            { "test",            required_argument, 0, 'T' },
#ifdef HAVE_AMR_MG_SOLVER
            { "use-trilinos",    no_argument,       0, 'a' },
            { "nsweeps",         required_argument, 0, 'g' },
            { "smoother",        required_argument, 0, 'q' },
            { "prec",            required_argument, 0, 'o' },
            { "basesolver",      required_argument, 0, 'u' },
            { "rebalance",       no_argument,       0, 'r' },
#endif
            { "tagging",         required_argument, 0, 't' },
            { "tagging-factor",  required_argument, 0, 'f' },
            { "bboxincr",        required_argument, 0, 'b' },
            { 0,                 0,                 0,  0  }
        };
        
        int option_index = 0;
        
#ifdef HAVE_AMR_MG_SOLVER
        c = getopt_long(argc, argv, "ab:d:f:g:hl:m:o:q:rst:u:x:y:z:", long_options, &option_index);
#else
        c = getopt_long(argc, argv, "b:x:y:z:l:m:r:hst:f:", long_options, &option_index);
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
                params.smoother = optarg;
                break;
            }
            case 'o':
            {
                params.prec = optarg;
                break;
            }
            case 'u':
            {
                params.bs = optarg;
                break;
            }
            case 'r':
            {
                params.rebalance = true;
                break;
            }
#endif
            case 'b':
                params.bboxincr = std::atof(optarg); ++cnt; break;
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
            case 'd':
                params.blocking_factor = std::atoi(optarg); ++cnt; break;
            case 'S':
                params.timestep = std::atof(optarg); ++cnt; break;
            case 'I':
                params.nIterations = std::atoi(optarg); ++cnt; break;
            case 'T':
            {
                if ( std::strcmp("twostream", optarg) == 0 )
                    params.type = Distribution::Type::kTwoStream;
                else if ( std::strcmp("recurrence", optarg) == 0 )
                    params.type = Distribution::Type::kRecurrence;
                else if ( std::strcmp("landau", optarg) == 0 )
                    params.type = Distribution::Type::kLandauDamping;
                else {
                    msg << "No such plasma test. Check '--test' option" << endl;
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
                ++cnt;
                break;
            }
            case 's':
                params.useMgtSolver = true;
                break;
            case 't':
            {
                if ( std::strcmp("efield", optarg) == 0 )
                    params.criteria = AmrOpal::kEfieldStrength;
                else if ( std::strcmp("potential", optarg) == 0 )
                    params.criteria = AmrOpal::kPotentialStrength;
                break;
            }
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
                    << "--blocking_factor [val] (only grids modulo bf == 0 allowed)" << endl
                    << "--bboxincr [value] (increase box size by [value] percent)" << endl
                    << "--use-mgt-solver (optional)" << endl
                    << "--timstep [val > 0]" << endl
                    << "--iterations [val > 0]" << endl
                    << "--test [twostream, recurrence, landau]" << endl
#ifdef HAVE_AMR_MG_SOLVER
                    << "--use-trilinos (optional)" << endl
                    << "--nsweeps (optional, trilinos only, default: 12)" << endl
                    << "--smoother (optional, trilinos only, default: GAUSS_SEIDEL)" << endl
                    << "--prec (optional, trilinos only, default: NONE)" << endl
                    << "--basesolver (optional, trilinos only, default: CG)" << endl
                    << "--rebalance (optional, trilinos only)" << endl
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



using namespace amrex;

// ----------------------------------------------------------------------------

#include "Particle/BoxParticleCachingPolicy.h"
#define Dim 3
typedef Cell                                                          Center_t;
typedef IntCIC                                                        IntrplCIC_t;
typedef UniformCartesian<2, double>                                   Mesh2d_t;
typedef CenteredFieldLayout<2, Mesh2d_t, Center_t>                    FieldLayout2d_t;
typedef Field<double, 2, Mesh2d_t, Center_t>                          Field2d_t;

void writeParticles(amrbunch_t* bunch, int step) {
    std::ofstream csvout;
    csvout.precision(10);
    csvout.setf(std::ios::scientific, std::ios::floatfield);

    std::stringstream fname;
    fname << "data/charges_nod_";
    fname << std::setw(1) << Ippl::myNode();
    fname << std::setw(5) << "_it_";
    fname << std::setw(4) << std::setfill('0') << step;
    fname << ".csv";
    
    // open a new data file for this iteration
    // and start with header
    csvout.open(fname.str().c_str(), std::ios::out);
    csvout << "x coord, y coord, z coord, charge, EfieldMagnitude, ID, vx, vy, vz, f" << std::endl;
    
    for (std::size_t i = 0; i< bunch->getLocalNum(); ++i) {
        double distributionf = bunch->P[i][2] * bunch->P[i][2] * 
                               std::exp( - ( bunch->P[i][0] * bunch->P[i][0] +
                                             bunch->P[i][1] * bunch->P[i][1] +
                                             bunch->P[i][2] * bunch->P[i][2] ) /2.
                                       );
        
        csvout << bunch->R[i](0) << ","
               << bunch->R[i](1) << ","
               << bunch->R[i](2) << ","
               << bunch->qm[i] << ","
               << std::sqrt( bunch->E[i][0] * bunch->E[i][0] +
                             bunch->E[i][1] * bunch->E[i][1] +
                             bunch->E[i][2] * bunch->E[i][2]
                           )
               << ","
               << bunch->ID[i] << ","
               << bunch->P[i][0] << ","
               << bunch->P[i][1] << ","
               << bunch->P[i][2] << ","
               << distributionf << std::endl;
    }
    csvout << std::endl;
    
    csvout.close();
}

void writeEnergy(amrbunch_t* bunch,
                 container_t& rho,
                 container_t& phi,
                 container_t& efield,
                 const std::vector<int>& rr,
                 double cell_volume,
                 int step,
                 std::string dir = "./")
{
    for (int lev = efield.size() - 2; lev >= 0; lev--)
        amrex::average_down(*(efield[lev+1].get()), *(efield[lev].get()), 0, 3, rr[lev]);
    
    // field energy (Ulmer version, i.e. cell_volume instead #points)
    double field_energy = 0.5 * cell_volume * MultiFab::Dot(*(efield[0].get()), 0, *(efield[0].get()), 0, 3, 0);
    
    // kinetic energy
    double ekin = 0.5 * sum( dot(bunch->P, bunch->P) );
    
    // potential energy
    rho[0]->mult(cell_volume, 0, 1);
    MultiFab::Multiply(*(phi[0].get()), *(rho[0].get()), 0, 0, 1, 0);
    double integral_phi_m = 0.5 * phi[0]->sum(0);
    
    if(Ippl::myNode()==0) {
        std::ofstream csvout;
        csvout.precision(10);
        csvout.setf(std::ios::scientific, std::ios::floatfield);

        std::stringstream fname;
        fname << dir << "/energy";
        fname << ".csv";

        // open a new data file for this iteration
        // and start with header
        csvout.open(fname.str().c_str(), std::ios::out | std::ofstream::app);
        
        if (step == 0) {
            csvout << "it,Efield,Ekin,Etot,Epot" << std::endl;
        }
        
        csvout << step << ", "
               << field_energy << ","
               << ekin << ","
               << field_energy + ekin << "," 
               << integral_phi_m << std::endl;
        
        csvout.close();
    }
}

void writeGridSum(container_t& container, int step, std::string label, std::string dir = "./") {
    std::ofstream csvout;
    csvout.precision(10);
    csvout.setf(std::ios::scientific, std::ios::floatfield);
    
    std::stringstream fname;
    fname << dir << "/" << label << "Sum";
    fname << ".csv";

    // open a new data file for this iteration
    // and start with header
    csvout.open(fname.str().c_str(), std::ios::out | std::ofstream::app);
    
    if ( step == 0 )
        csvout << "it,FieldSum" << std::endl;
    
    csvout << step << ", "<< container[0]->sum(0) << std::endl;
    csvout.close();
}

void ipplProjection(Field2d_t& field,
                    const Vector_t& dx,
                    const Vector_t& dv,
                    const Vector_t& Vmax,
                    const NDIndex<2>& lDom,
                    amrbunch_t* bunch, int step,
                    std::string dir = "./")
{
    field = 0; // otherwise values are accumulated over steps
    
    for (unsigned i=0; i< bunch->getLocalNum(); ++i)
        bunch->Rphase[i]= Vektor<double,2>(bunch->R[i][2],bunch->P[i][2]);
    
    bunch->qm.scatter(field, bunch->Rphase, IntrplCIC_t());
    
    std::ofstream csvout;
    csvout.precision(10);
    csvout.setf(std::ios::scientific, std::ios::floatfield);

    std::stringstream fname;
    fname << dir << "/f_mesh_";
    fname << std::setw(4) << std::setfill('0') << step;
    fname << ".csv";
    
    // open a new data file for this iteration
    // and start with header
    csvout.open(fname.str().c_str(), std::ios::out);
    csvout << "z, vz, f" << std::endl;
    
    for (int i=lDom[0].first(); i<=lDom[0].last(); i++) {
    
        for (int j=lDom[1].first(); j<=lDom[1].last(); j++) {
        
            csvout << (i+0.5) * dx[2] << ","
                   << (j+0.5) * dv[2] - Vmax[2]
                   << "," << field[i][j].get() << std::endl;
        }
    }
    
    csvout.close();
}

// ------------------------------------------------------------------------------------------


void doSolve(AmrOpal& myAmrOpal, amrbunch_t* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& efield,
             const amrex::Vector<int>& rr,
             Inform& msg,
             const double& scale, const param_t& params,
             std::string dir = "./")
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
    
    const amrex::Vector<amrex::Geometry>& geom = myAmrOpal.Geom();
    
    for (uint i = 0; i < rhs.size(); ++i)
        if ( rhs[i]->contains_nan() || rhs[i]->contains_inf() )
            throw std::runtime_error("\033[1;31mError: NANs or INFs on charge grid.\033[0m");
    
    
    // Check charge conservation
    double totCharge = totalCharge(rhs, finest_level, geom);
    double totCharge_composite = totalCharge_composite(rhs, finest_level, geom);
    
    msg << "Total Charge (computed):  " << totCharge << " C" << endl
        << "Total Charge (composite): " << totCharge_composite << " C" << endl
        << "Vacuum permittivity:      " << Physics::epsilon_0 << " F/m (= C/(m V)" << endl;
    
    amrex::Real vol = (*(geom[0].CellSize()) * *(geom[0].CellSize()) * *(geom[0].CellSize()) );
    msg << "Cell volume: " << *(geom[0].CellSize()) << "^3 = " << vol << " m^3" << endl;
    
    // eps in C / (V * m)
    double constant = -1.0; // / Physics::epsilon_0 ; //* scale;  // in [V m / C]
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
    
    if ( geom[0].isAllPeriodic() ) {
        double sum = rhs[0]->sum(0);
        double max = rhs[0]->max(0);
        msg << "total charge in density field before ion subtraction is " << sum << endl;
        msg << "max total charge in density field before ion subtraction is " << max << endl;
//         for (std::size_t i = 0; i < bunch->getLocalNum(); ++i)
//             offset += bunch->qm[i];
        
//         offset /= geom[0].ProbSize();
        offset = -1.0;
    }

    // solve
#ifdef HAVE_AMR_MG_SOLVER
    if ( params.useTrilinos ) {
	std::string interp = "PC";
        std::string norm = "LINF";
        AmrMultiGrid sol(&myAmrOpal, params.bs, params.prec,
                         params.rebalance, "periodic",
                         "periodic", "periodic",
                         params.smoother, params.nsweeps,
                         interp, norm, 0, 0);
        
        IpplTimings::startTimer(solvTimer);
        
        sol.solve(rhs,            // [V m]
                  phi,            // [V m^3]
                  efield,       // [V m^2]
                  base_level,
                  finest_level, false);
    
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

// void doSolve(AmrOpal& myAmrOpal, amrbunch_t* bunch,
//              container_t& rhs,
//              container_t& phi,
//              container_t& efield,
//              const Vector<Geometry>& geom,
//              const Vector<int>& rr,
//              int nLevels,
//              int step,
//              Inform& msg,
//              std::string dir = "./")
// {
//     // we need to write in the format of Ulmer
//     for (int i = 0; i <=finest_level; ++i) {
//         double cell_volume = geom[i].CellSize(0) *
//                              geom[i].CellSize(1) *
//                              geom[i].CellSize(2);
// #ifdef UNIQUE_PTR
//         rhs[i]->mult(cell_volume, 0, 1);       // in [V m]
// #else
//         rhs[i].mult(cell_volume, 0, 1);
// #endif
//     }
//     
//     writeScalarField(rhs, dir + "/rho", step);
//     
//     // undo Ulmer scaling
//     for (int i = 0; i <=finest_level; ++i) {
//         double cell_volume = geom[i].CellSize(0) *
//                              geom[i].CellSize(1) *
//                              geom[i].CellSize(2);
// #ifdef UNIQUE_PTR
//         rhs[i]->mult(1.0 / cell_volume, 0, 1);       // in [V m]
// #else
//         rhs[i].mult(1.0 / cell_volume, 0, 1);
// #endif
//     }
//     
//     // **************************************************************************                                                                                                                                
//     // Compute the total charge of all particles in order to compute the offset                                                                                                                                  
//     //     to make the Poisson equations solvable                                                                                                                                                                
//     // **************************************************************************                                                                                                                                
// 
//     Real offset = 0.0;
//     
//     if ( geom[0].isAllPeriodic() ) {
//         double sum = rhs[0]->sum(0);
//         double max = rhs[0]->max(0);
//         msg << "total charge in density field before ion subtraction is " << sum << endl;
//         msg << "max total charge in densitty field before ion subtraction is " << max << endl;
// //         for (std::size_t i = 0; i < bunch->getLocalNum(); ++i)
// //             offset += bunch->qm[i];
//         
// //         offset /= geom[0].ProbSize();
//         offset = -1.0;
//     }
// 
//     // solve                                                                                                                                                                                                     
//     Solver sol;
//     sol.solve_for_accel(rhs,
//                         phi,
//                         efield,
//                         geom,
//                         base_level,
//                         finest_level,
//                         offset,
//                         false);
//     
//     if ( geom[0].isAllPeriodic() ) {
//         double sum = rhs[0]->sum(0);
//         msg << "total charge in density field after ion subtraction is " << sum << endl;
//     }
//     
//     // for plotting unnormalize
//     for (int i = 0; i <=finest_level; ++i) {
// #ifdef UNIQUE_PTR
//         rhs[i]->mult(1.0 / constant, 0, 1);       // in [V m]
// #else
//         rhs[i].mult(1.0 / constant, 0, 1);
// #endif
//     }
// }


std::tuple<Vektor<std::size_t, 3>,
           Vektor<std::size_t, 3>,
           Vector_t,
           std::string
           >
initDistribution(const param_t& params,
                 std::unique_ptr<amrbunch_t>& bunch,
                 const Vector_t& extend_l,
                      const Vector_t& extend_r)
{
    Distribution dist;
    
    Vektor<std::size_t, 3> Nx, Nv;
    Vector_t Vmax;
    
    std::string dirname = "";
    
    if ( params.type == Distribution::Type::kTwoStream ) {
        dirname = "twostream";
        Nx = Vektor<std::size_t, 3>(4, 4, 32);
        Nv = Vektor<std::size_t, 3>(8, 8, 128);
        Vmax = Vector_t(6.0, 6.0, 6.0);
        
        dist.special(extend_l,
                     extend_r,
                     Nx,
                     Nv,
                     Vmax,
                     params.type,
                     0.05);
    } else if ( params.type == Distribution::Type::kRecurrence ) {
        dirname = "recurrence";
        Nx = Vektor<std::size_t, 3>(8, 8, 8);
        Nv = Vektor<std::size_t, 3>(32, 32, 32);
        Vmax = Vector_t(6.0, 6.0, 6.0);
        
        dist.special(extend_l,
                     extend_r,
                     Nx,
                     Nv,
                     Vmax,
                     params.type,
                     0.01);
    } else if ( params.type == Distribution::Type::kLandauDamping ) {
        dirname = "landau";
        Nx = Vektor<std::size_t, 3>(8, 8, 8);
        Nv = Vektor<std::size_t, 3>(32, 32, 32);
        Vmax = Vector_t(6.0, 6.0, 6.0);
        
        dist.special(extend_l,
                     extend_r,
                     Nx,
                     Nv,
                     Vmax,
                     params.type,
                     0.05);
    }
    
    dirname += (
        "-data-grid-" +
        std::to_string(params.nr[0]) + "-" + 
        std::to_string(params.nr[1]) + "-" + 
        std::to_string(params.nr[2])
    );
    
    boost::filesystem::path dir(dirname);
    if ( Ippl::myNode() == 0 )
        boost::filesystem::create_directory(dir);
    
    
    // copy particles to the PartBunchAmr object.
    dist.injectBeam( *(bunch.get()) );
    
    return std::make_tuple(Nx, Nv, Vmax, dir.string());
}


void doPlasma(const param_t& params, Inform& msg)
{
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    /*
     * create an Amr object
     */
    amrex::ParmParse pp("amr");
    pp.add("max_grid_size", int(params.maxBoxSize));
    
    amrex::Vector<int> error_buf(params.nLevels, 0);
    
    amrex::Vector<int> bf(params.nLevels, int(params.blocking_factor));
    pp.addarr("blocking_factor", bf);
    
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
    amrex::ParmParse pgeom("geometry");
    amrex::Vector<int> is_per = { 1, 1, 1};
    pgeom.addarr("is_periodic", is_per);
    
    amrex::Vector<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = params.nr[i];
    
    
    std::vector<int> rr(params.nLevels);
    amrex::Vector<int> rrr(params.nLevels);
    for (unsigned int i = 0; i < params.nLevels; ++i) {
        rr[i] = 2;
        rrr[i] = 2;
    }
    
    amrex::RealBox amr_domain;
    
    double incr = params.bboxincr * 0.01;
    double blen = 1.0 + incr;
    
    std::array<double, AMREX_SPACEDIM> amr_lower = {{-blen, -blen, -blen}}; // m
    std::array<double, AMREX_SPACEDIM> amr_upper = {{ blen,  blen,  blen}}; // m
    
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
    
    const amrex::Vector<amrex::BoxArray>& ba = myAmrOpal.boxArray();
    const amrex::Vector<amrex::DistributionMapping>& dmap = myAmrOpal.DistributionMap();
    amrex::Vector<amrex::Geometry>& geom = myAmrOpal.Geom();
    
    
    amrplayout_t* playout = new amrplayout_t(geom, dmap, ba, rrr);
    
    std::unique_ptr<amrbunch_t> bunch( new amrbunch_t() );
    bunch->initialize(playout);
    bunch->initializeAmr(); // add attributes: level, grid
    
    bunch->setAllowParticlesNearBoundary(true);
        
    // initialize a particle distribution
    Vector_t extend_l = Vector_t(0.0, 0.0, 0.0);
    
    Vector_t extend_r = Vector_t(4.0 * Physics::pi,
                                 4.0 * Physics::pi,
                                 4.0 * Physics::pi);
    
    auto tuple = initDistribution(params, bunch, extend_l, extend_r);
    
    // --------------------------------------------------------------------
    
    Vektor<std::size_t, 3> Nx = std::get<0>(tuple);
    Vektor<std::size_t, 3> Nv = std::get<1>(tuple);
    Vector_t Vmax = std::get<2>(tuple);
    std::string dir = std::get<3>(tuple);
    
    double spacings[2] = {
        ( extend_r[2] - extend_l[2] ) / Nx[2],
        2. * Vmax[2] / Nv[2]
    };
    
    Vektor<double,2> origin = {
        extend_l[2],
        -Vmax[2]
    };
    
    Index I(Nx[2]+1);
    Index J(Nv[2]+1);
    NDIndex<2> domain2d;
    domain2d[0]=I;
    domain2d[1]=J;
    
    Mesh2d_t mesh2d = Mesh2d_t(domain2d, spacings, origin);
//     std::unique_ptr<FieldLayout2d_t> layout2d = std::unique_ptr<FieldLayout2d_t>( new FieldLayout2d_t(mesh2d) );
    FieldLayout2d_t* layout2d = new FieldLayout2d_t(mesh2d);
    
    mesh2d.set_meshSpacing(&(spacings[0]));
    mesh2d.set_origin(origin);

    domain2d = layout2d->getDomain();
    
    BConds<double, 2, Mesh2d_t, Center_t> BC;
    if (Ippl::getNodes()>1) {
        BC[0] = new ParallelInterpolationFace<double,2,Mesh2d_t, Center_t>(0);
        BC[1] = new ParallelInterpolationFace<double,2,Mesh2d_t, Center_t>(1);
        BC[2] = new ParallelInterpolationFace<double,2,Mesh2d_t, Center_t>(2);
        BC[3] = new ParallelInterpolationFace<double,2,Mesh2d_t, Center_t>(3);
    }
    else {
        BC[0] = new InterpolationFace<double,2,Mesh2d_t, Center_t>(0);
        BC[1] = new InterpolationFace<double,2,Mesh2d_t, Center_t>(1);
        BC[2] = new InterpolationFace<double,2,Mesh2d_t, Center_t>(2);
        BC[3] = new InterpolationFace<double,2,Mesh2d_t, Center_t>(3);
    }
    
    NDIndex<2> lDom = domain2d;
    Vektor<double,3> dx = (extend_r - extend_l) / Vector_t(Nx);
    Vektor<double,3> dv = 2. * Vmax / Vector_t(Nv);
    
    
    // field is used for twostream instability as 2D phase space mesh
    Field2d_t field;
//     field.initialize(mesh2d, *(layout2d.get()), GuardCellSizes<2>(1),BC);
    field.initialize(mesh2d, *layout2d, GuardCellSizes<2>(2),BC);
    
    
    
    
    // --------------------------------------------------------------------
    
    container_t rhs;
    container_t phi;
    container_t efield;
    
    
    // map particles
    double scale = 1.0;
    
    scale = domainMapping(*bunch, scale);
    
    msg << "Scale: " << scale << endl;
    msg << endl << "Transformed positions" << endl << endl;
    
    // redistribute on single-level
    bunch->update();
    
    myAmrOpal.setBunch( bunch.get() );
    
    const Vector<Geometry>& geoms = myAmrOpal.Geom();
    
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    
    msg << "Multi-level statistics" << endl;
    bunch->gatherStatistics();
    
    msg << "Total number of particles: " << bunch->getTotalNum() << endl;
    
    doSolve(myAmrOpal, bunch.get(), rhs, phi, efield, rrr, msg, scale, params, dir);
    
    msg << endl << "Back to normal positions" << endl << endl;
    
    domainMapping(*bunch, scale, true);
    
    
//     writeScalarField(rhs, dir + "/rho_0.dat");
    
//     std::string plotsolve = amrex::Concatenate(dir + "/plt", 0, 4);
    
//     writePlotFile(plotsolve, rhs, phi, efield, rr, geoms, 0);
    
    Vector_t hr = ( extend_r - extend_l ) / Vector_t(params.nr);
    double cell_volume = hr[0] * hr[1] * hr[2];
    writeEnergy(bunch.get(), rhs, phi, efield, rr, cell_volume, 0, dir);
    
//     writeGridSum(rhs, 0, "RhoInterpol", dir);
//     writeGridSum(phi, 0, "Phi_m", dir);
    
    for (std::size_t i = 0; i < params.nIterations; ++i) {
        msg << "Processing step " << i << endl;

        if ( params.type == Distribution::Type::kTwoStream )
            ipplProjection(field, dx, dv, Vmax, lDom, bunch.get(), i, dir);
        
        assign(bunch->R, bunch->R + params.timestep * bunch->P);
        
        // periodic shift
        double up = 4.0 * Physics::pi;
        for (std::size_t j = 0; j < bunch->getLocalNum(); ++j) {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                if ( bunch->R[j](d) > up )
                    bunch->R[j](d) = bunch->R[j](d) - up;
                else if ( bunch->R[j](d) < 0.0 )
                    bunch->R[j](d) = up + bunch->R[j](d);
            }
        }
        
        msg << endl << "Transformed positions" << endl << endl;
        scale = domainMapping(*bunch, scale);
        
        if ( myAmrOpal.maxLevel() > 0 )
            for (int k = 0; k <= myAmrOpal.finestLevel() && k < myAmrOpal.maxLevel(); ++k)
                myAmrOpal.regrid(k /*lbase*/, i /*time*/);
        else
            bunch->update();
        
        
        doSolve(myAmrOpal, bunch.get(), rhs, phi, efield, rrr, msg, scale, params, dir);
        
        bunch->GetGravity(bunch->E, efield);
        
        domainMapping(*bunch, scale, true);
        msg << endl << "Back to normal positions" << endl << endl;
    
        
        
        /* epsilon_0 not used by Ulmer --> multiply it away */
        assign(bunch->P, bunch->P + params.timestep * bunch->qm / bunch->mass * bunch->E ); //* Physics::epsilon_0);
        
        writeEnergy(bunch.get(), rhs, phi, efield, rr, cell_volume, i + 1, dir);
//         writeGridSum(rhs, i + 1, "RhoInterpol", dir);
//         writeGridSum(phi, i + 1, "Phi_m", dir);
        
        msg << "Done with step " << i << endl;
    }
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    amrex::Initialize(argc,argv, false);
    
    Inform msg(argv[0]);
    
    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    
    IpplTimings::startTimer(mainTimer);
    
    std::string test = "";
    
    param_t params;
    
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
        
        if ( params.type == Distribution::Type::kTwoStream )
            test = "twostream";
        else if ( params.type == Distribution::Type::kRecurrence )
            test = "recurrence";
        else if ( params.type == Distribution::Type::kLandauDamping )
            test = "landau";
    
        msg << "Plasma test running with" << endl
            << "- grid                  = " << params.nr << endl
            << "- max. grid             = " << params.maxBoxSize << endl
            << "- blocking factor       = " << params.blocking_factor << endl
            << "- #level                = " << params.nLevels - 1 << endl
            << "- time step             = " << params.timestep << endl
            << "- #steps                = " << params.nIterations << endl
            << "- test                  = " << test << endl
            << "- bboxincr              = " << params.bboxincr << endl
            << "- tagging               = " << tagging << endl
            << "- tagging factor        = " << params.tagfactor << endl;

        if ( params.useMgtSolver )
            msg << "- MGT solver is used" << endl;
        
#ifdef HAVE_AMR_MG_SOLVER
        if ( params.useTrilinos ) {
            msg << "- Trilinos solver is used with: " << endl
                << "    - nsweeps:     " << params.nsweeps << endl
                << "    - smoother:    " << params.smoother
                << endl;
        }
#endif
        doPlasma(params, msg);
        
    } catch(const std::exception& ex) {
        msg << ex.what() << endl;
    }
    
    IpplTimings::stopTimer(mainTimer);

    std::string fn = std::string(argv[0]) + "-timing.dat";
    
    std::map<std::string, unsigned int> problemSize;
    
    problemSize["level"]       = params.nLevels;
    problemSize["maxgrid"]     = params.maxBoxSize;
    problemSize["gridx"]       = params.nr[0];
    problemSize["gridy"]       = params.nr[1];
    problemSize["gridz"]       = params.nr[2];
    problemSize["#iterations"] = params.nIterations;
    
    IpplTimings::print(fn, problemSize);
    
    return 0;
}
