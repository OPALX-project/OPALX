#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <sstream>
#include <iomanip>
#include <functional>

#include <tuple>

#include <boost/filesystem.hpp>


#include <AMReX_ParmParse.H>

#include "../Distribution.h"
#include "../Solver.h"
#include "../MGTSolver.h"
#include "../AmrOpal.h"

#include "../helper_functions.h"

#include "../writePlotFile.H"

#include "../H5Reader.h"


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
    Vektor<size_t, AMREX_SPACEDIM> nr;
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
#endif
    double bboxincr;
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
    
    int required = 9;
    
#if AMREX_SPACEDIM == 3
    ++required;
#endif
    
    while ( true ) {
        static struct option long_options[] = {
            { "gridx",           required_argument, 0, 'x' },
            { "gridy",           required_argument, 0, 'y' },
#if AMREX_SPACEDIM == 3
            { "gridz",           required_argument, 0, 'z' },
#endif
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
#if AMREX_SPACEDIM == 3
            case 'z':
                params.nr[2] = std::atoi(optarg); ++cnt; break;
#endif
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
#if AMREX_SPACEDIM == 3
                    << "--gridz [#gridpoints in z]" << endl
#endif
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

void writeParticles(std::string file, amrbunch_t* bunch, unsigned int step) {
    
    if ( step == 0 ) {
        H5Reader h5(file);
        h5.open(step, H5_O_WRONLY);
        h5.writeHeader();
        h5.write(bunch, step);
        h5.close();
    } else {
        H5Reader h5(file);
        h5.open(step, H5_O_APPENDONLY);
        h5.write(bunch, step);
        h5.close();
    }
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
        amrex::average_down(*(efield[lev+1].get()), *(efield[lev].get()), 0,
                            AMREX_SPACEDIM, rr[lev]);
    
    // field energy (Ulmer version, i.e. cell_volume instead #points)
    double field_energy = 0.5 * cell_volume * MultiFab::Dot(*(efield[0].get()), 0, *(efield[0].get()),
                                                            0, AMREX_SPACEDIM, 0);
//     double field_energy = 0.5 * cell_volume * sum(dot(bunch->E, bunch->E));
    
    // kinetic energy
    double ekin = 0.5 * sum( dot(bunch->P, bunch->P) );
    
    // potential energy
    rho[0]->mult(cell_volume, 0, 1);
    MultiFab::Multiply(*(phi[0].get()), *(rho[0].get()), 0, 0, 1, 0);
    double integral_phi_m = 0.5 * phi[0]->sum(0);
    
    double potential_energy=0;
    for (unsigned i=0; i<bunch->getLocalNum(); ++i) {
        potential_energy += 0.5 * (bunch->qm[i] *  bunch->phi[i]);
    }
    
        // open a new data file for this iteration
        // and start with header
    double AmplitudeEfield = max(sqrt(dot(bunch->E,bunch->E)));
    bunch->E1 = 0;
    assign(bunch->E1, bunch->E * AMREX_D_PICK(1.0, Vector_t(0.0, 1.0), Vector_t(0,0,1)));
    double AmplitudeEFz=max(sqrt(dot(bunch->E1,bunch->E1)));
    
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
               << /*field_energy*/potential_energy + ekin << "," 
               << /*integral_phi_m*/potential_energy << std::endl;
        
        csvout.close();
        
        
        fname.str("");
        fname << dir << "/amplitude";
        fname << ".csv";

        
        csvout.open(fname.str().c_str(), std::ios::out | std::ofstream::app);
        if (step == 0){
                csvout << "it,max(|E|),max(|Ez|)" << std::endl;
        }
        csvout << step << ", "
                << AmplitudeEfield << ", "
                << AmplitudeEFz << std::endl;
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
        bunch->Rphase[i]= Vektor<double,2>(bunch->R[i][AMREX_SPACEDIM-1], bunch->P[i][AMREX_SPACEDIM-1]);
    
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
        
            csvout << (i+0.5) * dx[AMREX_SPACEDIM-1] << ","
                   << (j+0.5) * dv[AMREX_SPACEDIM-1] - Vmax[AMREX_SPACEDIM-1]
                   << "," << field[i][j].get() << std::endl;
        }
    }
    
    csvout.close();
}

// ------------------------------------------------------------------------------------------


void doSolve(AmrOpal& myAmrOpal, amrbunch_t* bunch,
//              amrplayout_t::ParticlePos_t& E,
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
    
    amrex::Real vol = AMREX_D_TERM(*(geom[0].CellSize()), * *(geom[0].CellSize()), * *(geom[0].CellSize()));
    msg << "Cell volume: " << *(geom[0].CellSize()) << "^"
        << AMREX_SPACEDIM << " = " << vol << " m^" << AMREX_SPACEDIM << endl;
    
    // eps in C / (V * m)
    double constant = -1.0; // / Physics::epsilon_0 ; //* scale;  // in [V m / C]
    
    msg << "total charge in density field before ion subtraction is "
        << totCharge_composite / vol
        << endl;
    
    //subtract the background charge of the ions
    for (int i = 0; i <= finest_level; ++i) {
//         rhs[i]->mult(-1.0, 0, 1);
        rhs[i]->plus(1.0, 0, 1);
    }
    
    totCharge_composite = totalCharge_composite(rhs, finest_level, geom);
    msg << "total charge in density field after ion subtraction is "
        << totCharge_composite / vol
        << endl;
    
    for (int i = 0; i <= finest_level; ++i) {
        rhs[i]->mult(constant, 0, 1);
    }
    
    
    
    // normalize each level
    double l0norm = rhs[finest_level]->norm0(0);
    msg << "l0norm = " << l0norm << endl;
    for (int i = 0; i <= finest_level; ++i) {
        rhs[i]->mult(1.0 / l0norm, 0, 1);
    }
    
    amrex::Real offset = 0.;

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
    
    bunch->phi = 0;
    
    bunch->InterpolateFort(bunch->phi, phi, base_level, finest_level);
    
    // undo scale
    for (int i = 0; i <= finest_level; ++i)
        efield[i]->mult(scale * scale * l0norm/*[i]*/, 0, AMREX_SPACEDIM);
    
    
    bunch->E = 0;
    
    bunch->InterpolateFort(bunch->E, efield, base_level, finest_level);
    
    for (int i = 0; i <= finest_level; ++i) {
        rhs[i]->mult(l0norm / constant, 0, 1);
    }
    
    IpplTimings::stopTimer(solvTimer);
    
    
}


std::tuple<Vektor<std::size_t, AMREX_SPACEDIM>,
           Vektor<std::size_t, AMREX_SPACEDIM>,
           Vector_t,
           std::string
           >
initDistribution(const param_t& params,
                 std::unique_ptr<amrbunch_t>& bunch,
                 const Vector_t& extend_l,
                 const Vector_t& extend_r)
{
    Distribution dist;
    
    Vektor<std::size_t, AMREX_SPACEDIM> Nx, Nv;
    Vector_t Vmax;
    
    std::string dirname = "";
    
    if ( params.type == Distribution::Type::kTwoStream ) {
        dirname = "twostream";
//         Nx = Vektor<std::size_t, AMREX_SPACEDIM>(D_DECL(4, 4, 32)); // 4, 4, 32
//         Nv = Vektor<std::size_t, AMREX_SPACEDIM>(D_DECL(8, 8, 128)); // 8, 8, 128
//         Vmax = Vector_t(D_DECL(6.0, 6.0, 6.0));
        
        Nx = Vektor<std::size_t, AMREX_SPACEDIM>(D_DECL(4, 384, 32)); // 4, 4, 32
        Nv = Vektor<std::size_t, AMREX_SPACEDIM>(D_DECL(4, 384, 64)); // 8, 8, 128
        Vmax = Vector_t(D_DECL(9.0, 9.0, 9.0));
        
        dist.special(extend_l,
                     extend_r,
                     Nx,
                     Nv,
                     Vmax,
                     params.type,
                     0.05);
    } else if ( params.type == Distribution::Type::kRecurrence ) {
        dirname = "recurrence";
        Nx = Vektor<std::size_t, AMREX_SPACEDIM>(D_DECL(8, 8, 8));
        Nv = Vektor<std::size_t, AMREX_SPACEDIM>(D_DECL(32, 32, 32));
        Vmax = Vector_t(D_DECL(6.0, 6.0, 6.0));
        
        dist.special(extend_l,
                     extend_r,
                     Nx,
                     Nv,
                     Vmax,
                     params.type,
                     0.01);
    } else if ( params.type == Distribution::Type::kLandauDamping ) {
        dirname = "landau";
//         Nx = Vektor<std::size_t, AMREX_SPACEDIM>(D_DECL(8, 8, 8));
//         Nv = Vektor<std::size_t, AMREX_SPACEDIM>(D_DECL(32, 32, 32));
        Nx = Vektor<std::size_t, AMREX_SPACEDIM>(D_DECL(8, 256, 32));
        Nv = Vektor<std::size_t, AMREX_SPACEDIM>(D_DECL(8, 256, 64));
        Vmax = Vector_t(D_DECL(6.0, 6.0, 6.0));
        
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
        AMREX_D_TERM(std::to_string(params.nr[0]),
                     + "-" + std::to_string(params.nr[1]),
                     + "-" + std::to_string(params.nr[2]))
    );
    
    boost::filesystem::path dir(dirname);
    if ( Ippl::myNode() == 0 )
        boost::filesystem::create_directory(dir);
    
    
    // copy particles to the PlasmaBunch object.
    dist.injectBeam( *(bunch.get()) );
    
    return std::make_tuple(Nx, Nv, Vmax, dir.string());
}



// void remapping(std::unique_ptr<amrbunch_t>& bunch,
//                const Vector_t& lower,
//                const Vector_t& upper)
// {
//     double vmax = max(bunch->P)[1];
//     
//     Vector_t nx = Vektor<std::size_t, AMREX_SPACEDIM>(D_DECL(4, 32, 32));
//     Vector_t nv = Vektor<std::size_t, AMREX_SPACEDIM>(D_DECL(4, 32, 64));
//     
//     Vektor<double, AMREX_SPACEDIM> hx = (upper - lower) / Vector_t(nx);
//     Vektor<double, AMREX_SPACEDIM> hv = 2.0 * vmax / Vector_t(nv);
//     
//     std::function<double(const Vector_t&)> W3 = [&](const Vector_t& y) {
//         double res = 1.0;
//         for (int d = 0; d < AMREX_SPACEDIM; ++d) {
//             double x = std::abs(y[d]);
//             if ( x <= 1 && x >= 0 ) {
//                 res *= (1.0 - 2.5 * x * x + 1.5 * x * x * x);
//             } else if (  x <= 2 && x >= 1 ) {
//                 res *= 0.5 * (2.0 - x) * (2.0 - x) * (1.0 - x);
//             } else {
//                 res *= 0.0;
//             }
//         }
//         return res;
//     };
//     
//     Distribution::container_t x_m;    ///< Horizontal particle positions [m]
//     Distribution::container_t y_m;    ///< Vertical particle positions [m]
//     Distribution::container_t z_m;    ///< Longitudinal particle positions [m]
// 
//     Distribution::container_t px_m;   ///< Horizontal particle momentum
//     Distribution::container_t py_m;   ///< Vertical particle momentum
//     Distribution::container_t pz_m;   ///< Longitudinal particle momentum
//     Distribution::container_t q_m;    ///< Particle charge (always set to 1.0, except for Distribution::readH5)
//     Distribution::container_t mass_m;
//     
//     std::cout << "    Start ..." << std::endl;
//     
//     for (std::size_t i = 0; i < nx[0]; ++i) {
//             for (std::size_t j = 0; j < nx[1]; ++j) {
// #if AMREX_SPACEDIM == 3
//                 for ( std::size_t k = 0; k < nx[2]; ++k) {
// #endif
//                     Vektor<double, AMREX_SPACEDIM> pos = Vektor<double,AMREX_SPACEDIM>(
//                                         D_DECL(
//                                                 (0.5 + i) * hx[0] + lower[0],
//                                                 (0.5 + j) * hx[1] + lower[1],
//                                                 (0.5 + k) * hx[2] + lower[2]
//                                         )
//                     );
//                     
//                     for (std::size_t iv = 0; iv < nv[0]; ++iv) {
//                         for (std::size_t jv = 0; jv < nv[1]; ++jv) {
// #if AMREX_SPACEDIM == 3
//                             for (std::size_t kv = 0; kv < nv[2]; ++kv) {
// #endif
//                                 Vektor<double, AMREX_SPACEDIM> vel = -vmax + hv *
//                                                         Vektor<double, AMREX_SPACEDIM>(
//                                                             D_DECL(iv + 0.5,
//                                                                    jv + 0.5,
//                                                                    kv + 0.5
//                                                                    )
//                                 );
//                                 
//                                     x_m.push_back( pos[0] );
//                                     y_m.push_back( pos[1] );
// #if AMREX_SPACEDIM == 3
//                                     z_m.push_back( pos[2] );
// #endif
//                                     px_m.push_back( vel[0] );
//                                     py_m.push_back( vel[1] );
// #if AMREX_SPACEDIM == 3
//                                     pz_m.push_back( vel[2] );
// #endif             
//                                     double tmp = 0;
//                                     for (std::size_t ip = 0; ip < bunch->getLocalNum(); ++ip) {
//                                         tmp += bunch->qm[ip] *
//                                                W3((bunch->R[ip] - pos) / hx) *
//                                                W3((bunch->P[ip] - vel) / hv);
//                                         
//                                     }
//                                     
//                                     q_m.push_back(tmp);
//                                     mass_m.push_back( -tmp );
// #if AMREX_SPACEDIM == 3
//                             }
// #endif
//                         }
//                     }
// #if AMREX_SPACEDIM == 3
//                 }
// #endif
//             }
//         }
//         
//         std::cout << "    Done." << std::endl;
//         
//         for (int i = bunch->getLocalNum()-1; i > -1; --i) {
//             bunch->R[i] = Vector_t(D_DECL(x_m[i], y_m[i], z_m[i]));
//             bunch->P[i] = Vector_t(D_DECL(px_m[i], py_m[i], pz_m[i]));
//             bunch->qm[i] = q_m[i];
//             bunch->mass[i] = mass_m[i];
//             
//             x_m.pop_back();
//             y_m.pop_back();
// #if AMREX_SPACEDIM == 3
//         z_m.pop_back();
// #endif
//         px_m.pop_back();
//         py_m.pop_back();
// #if AMREX_SPACEDIM == 3
//         pz_m.pop_back();
// #endif
//         q_m.pop_back();
//         mass_m.pop_back();
//     }
//     
//     bunch->update();
// }


void updateIpplMesh(Field2d_t* field,
                    FieldLayout2d_t* layout2d,
                    BConds<double, 2, Mesh2d_t, Center_t>& BC,
                    std::unique_ptr<amrbunch_t> &bunch,
                    int nx, int nv)
{
    Mesh2d_t& mesh = field->get_mesh();
    
    
//     Mesh2d_t(domain2d, spacings, origin);
    
//     Vector_t rmin = min(bunch->R);
    Vector_t pmin = min(bunch->P);
    Vector_t pmax = max(bunch->P);
    
    Vektor<double, 2> origin(0, pmin[AMREX_SPACEDIM-1]);
    
    double spacings[2] = {
        ( 4.0 * Physics::pi ) / nx,
        (pmax[AMREX_SPACEDIM-1] - pmin[AMREX_SPACEDIM-1]) / nv
    };
    
    mesh.set_meshSpacing(&(spacings[0]));
    mesh.set_origin(origin);
    
    field->initialize(mesh,
                      *layout2d,
                      GuardCellSizes<2>(1),
                      BC);
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
    amrex::Vector<int> is_per = { D_DECL(1, 1, 1) };
    pgeom.addarr("is_periodic", is_per);
    
    amrex::Vector<int> nCells(AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; ++i)
        nCells[i] = params.nr[i];
    
    
    std::vector<int> rr(params.nLevels);
    amrex::Vector<int> rrr(params.nLevels);
    for (unsigned int i = 0; i < params.nLevels; ++i) {
        rr[i] = 2;
        rrr[i] = 2;
    }
    
    amrex::RealBox amr_domain;
    
//     double incr = params.bboxincr * 0.01;
//     double blen = 1.0 + incr;
    
//     std::array<double, AMREX_SPACEDIM> amr_lower = {{D_DECL(-blen, -blen, -blen)}}; // m
//     std::array<double, AMREX_SPACEDIM> amr_upper = {{D_DECL( blen,  blen,  blen)}}; // m
    
    std::array<double, AMREX_SPACEDIM> amr_lower = {{D_DECL(0.0, 0.0, 0.0)}}; // m
    std::array<double, AMREX_SPACEDIM> amr_upper = {{
        D_DECL(4.0 * Physics::pi,
               4.0 * Physics::pi,
               4.0 * Physics::pi)
    }}; // m
    
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
    Vector_t extend_l = Vector_t(D_DECL(0.0, 0.0, 0.0));
    
    Vector_t extend_r = Vector_t(D_DECL(4.0 * Physics::pi,
                                        4.0 * Physics::pi,
                                        4.0 * Physics::pi));
    
    auto tuple = initDistribution(params, bunch, extend_l, extend_r);
    
    // --------------------------------------------------------------------
    
    Vektor<std::size_t, AMREX_SPACEDIM> Nx = std::get<0>(tuple);
    Vektor<std::size_t, AMREX_SPACEDIM> Nv = std::get<1>(tuple);
    Vector_t Vmax = std::get<2>(tuple);
    std::string dir = std::get<3>(tuple);
    
    double spacings[2] = {
        ( extend_r[AMREX_SPACEDIM-1] - extend_l[AMREX_SPACEDIM-1] ) / Nx[AMREX_SPACEDIM-1],
        2. * Vmax[AMREX_SPACEDIM-1] / Nv[AMREX_SPACEDIM-1]
    };
    
    Vektor<double,2> origin = {
        extend_l[AMREX_SPACEDIM-1],
        -Vmax[AMREX_SPACEDIM-1]
    };
    
    Index I(Nx[AMREX_SPACEDIM-1]+1);
    Index J(Nv[AMREX_SPACEDIM-1]+1);
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
    Vektor<double,AMREX_SPACEDIM> dx = (extend_r - extend_l) / Vector_t(Nx);
//     Vektor<double,3> dv = 2. * Vmax / Vector_t(Nv);
    
    
    // field is used for twostream instability as 2D phase space mesh
    Field2d_t field;
//     field.initialize(mesh2d, *(layout2d.get()), GuardCellSizes<2>(1),BC);
    field.initialize(mesh2d, *layout2d, GuardCellSizes<2>(1),BC);
    
    
    
    
    // --------------------------------------------------------------------
    
    container_t rhs;
    container_t phi;
    container_t efield;
    
    
    // map particles
    double scale = 1.0;
    
//     scale = domainMapping(*bunch, scale);
    
//     msg << "Scale: " << scale << endl;
//     msg << endl << "Transformed positions" << endl << endl;
    
    // redistribute on single-level
    bunch->update();
    
    myAmrOpal.setBunch( bunch.get() );
    
    const Vector<Geometry>& geoms = myAmrOpal.Geom();
    
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    
//     msg << "Multi-level statistics" << endl;
//     bunch->gatherStatistics();
    
//     msg << "Total number of particles: " << bunch->getTotalNum() << endl;
    
//     msg << endl << "Back to normal positions" << endl << endl;
    
//     domainMapping(*bunch, scale, true);
    
    Vector_t hr = ( extend_r - extend_l ) / Vector_t(params.nr);
    double cell_volume = AMREX_D_TERM(*(geom[0].CellSize()),
                                      * *(geom[0].CellSize()),
                                      * *(geom[0].CellSize()) );//hr[0] * hr[1] * hr[2];
    
    for (std::size_t i = 0; i < params.nIterations; ++i) {
        msg << "Processing step " << i << endl;
        
//         if ( i != 0 && i % 5 == 0) {
//             msg << "    Particle remapping" << endl;
//             remapping(bunch, extend_l, extend_r);
//         }

        if ( params.type == Distribution::Type::kTwoStream ) {
            writeParticles(dir + "/particles.h5", bunch.get(), i);
            
            updateIpplMesh(&field, layout2d, BC, bunch, Nx[AMREX_SPACEDIM-1], Nv[AMREX_SPACEDIM-1]);
            
            Vmax = max(bunch->P);
            Vector_t dv = (Vmax - min(bunch->P)) / Vector_t(Nv);
            domain2d = layout2d->getDomain();
            NDIndex<2> lDom = domain2d;
            ipplProjection(field, dx, dv, Vmax, lDom, bunch.get(), i, dir);
        }
        
//         msg << endl << "Transformed positions" << endl << endl;
//         scale = domainMapping(*bunch, scale);
        
        doSolve(myAmrOpal, bunch.get(), /*bunch->E,*/ rhs, phi, efield, rrr, msg, scale, params, dir);
        
//         domainMapping(*bunch, scale, true);
//         msg << endl << "Back to normal positions" << endl << endl;
        
        writeEnergy(bunch.get(), rhs, phi, efield, rr, cell_volume, i, dir);
        /* Leap-Frog
         * 
         * v_{i+1/2} = v_i + a_i * dt / 2    (a_i = F_i = q * E_i )
         * x_{i+1} = x_i + v_{i+1/2} * dt
         * v_{i+1} = v_{i+1/2} + a_{i+1} * dt / 2
         */
//         assign(bunch->P, bunch->P + 0.5 * params.timestep * bunch->qm / bunch->mass * bunch->E);
        
        
        double up = 4.0 * Physics::pi;
        for (std::size_t j = 0; j < bunch->getLocalNum(); ++j) {
            bunch->P[j] = bunch->P[j] + 0.5 * params.timestep * bunch->qm[j] / bunch->mass[j] * bunch->E[j];
            bunch->R[j] = bunch->R[j] + params.timestep * bunch->P[j];
            
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                if ( bunch->R[j](d) > up )
                    bunch->R[j](d) = bunch->R[j](d) - up;
                else if ( bunch->R[j](d) < 0.0 )
                    bunch->R[j](d) = up + bunch->R[j](d);
            }
        }
        
//         msg << endl << "Transformed positions" << endl << endl;
//         scale = domainMapping(*bunch, scale);
        
        if ( myAmrOpal.maxLevel() > 0 )
            for (int k = 0; k <= myAmrOpal.finestLevel() && k < myAmrOpal.maxLevel(); ++k)
                myAmrOpal.regrid(k /*lbase*/, i /*time*/);
        else
            bunch->update();
        
        doSolve(myAmrOpal, bunch.get(), /*bunch->E,*/ rhs, phi, efield, rrr, msg, scale, params, dir);
        
//         domainMapping(*bunch, scale, true);
//         msg << endl << "Back to normal positions" << endl << endl;
    
        
        /* epsilon_0 not used by Ulmer --> multiply it away */
        for (std::size_t j = 0; j < bunch->getLocalNum(); ++j) {
            bunch->P[j] = bunch->P[j] + 0.5 * params.timestep * bunch->qm[j] / bunch->mass[j] * bunch->E[j];
            //* Physics::epsilon_0);
        }

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
#if AMREX_SPACEDIM == 3
    problemSize["gridz"]       = params.nr[2];
#endif
    problemSize["#iterations"] = params.nIterations;
    
    IpplTimings::print(fn, problemSize);
    
    amrex::Finalize(true);
    
    return 0;
}
