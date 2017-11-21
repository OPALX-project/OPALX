#include "Ippl.h"

#include <AMReX_ParmParse.H>

#include "AmrMultiGrid.h"

#include "../Solver.h"

#include "../writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"
#include <random>

#include <getopt.h>

#include <AMReX_MultiFab.H>

using namespace amrex;

typedef Array<std::unique_ptr<MultiFab> > container_t;

struct param_t {
    Vektor<size_t, AMREX_SPACEDIM> nr;
    size_t maxBoxSize;
    bool isWriteYt;
    bool isHelp;
    size_t nLevels;
};

bool parseProgOptions(int argc, char* argv[], param_t& params, Inform& msg) {
    /* Parsing Command Line Arguments
     * 
     * 26. June 2017
     * https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example
     */
    
    params.isWriteYt = false;
    params.isHelp = false;
    
    int c = 0;
    
    int cnt = 0;
    
    int required = 2 +  AMREX_SPACEDIM;
    
    while ( true ) {
        static struct option long_options[] = {
            { "gridx",          required_argument, 0, 'x' },
            { "gridy",          required_argument, 0, 'y' },
#if AMREX_SPACEDIM == 3
            { "gridz",          required_argument, 0, 'z' },
#endif
            { "level",          required_argument, 0, 'l' },
            { "maxgrid",        required_argument, 0, 'm' },
            { "writeYt",        no_argument,       0, 'w' },
            { "help",           no_argument,       0, 'h' },
            { 0,                0,                 0,  0  }
        };
        
        int option_index = 0;
        
        c = getopt_long(argc, argv,
#if AMREX_SPACEDIM == 3
                        "x:y:z:l:m:wh",
#else
                        "x:y:l:m:wh",
#endif
                        long_options, &option_index);
        
        if ( c == -1 )
            break;
        
        switch ( c ) {
            case 'x':
                params.nr[0] = std::atoi(optarg); ++cnt; break;
            case 'y':
                params.nr[1] = std::atoi(optarg); ++cnt; break;
#if AMREX_SPACEDIM == 3
            case 'z':
                params.nr[2] = std::atoi(optarg); ++cnt; break;
#endif
            case 'l':
                params.nLevels = std::atoi(optarg); ++cnt; break;
            case 'm':
                params.maxBoxSize = std::atoi(optarg); ++cnt; break;
            case 'w':
                params.isWriteYt = true;
                break;
            case 'h':
                msg << "Usage: " << argv[0]
                    << endl
                    << "--gridx [#gridpoints in x]" << endl
                    << "--gridy [#gridpoints in y]" << endl
#if AMREX_SPACEDIM == 3
                    << "--gridz [#gridpoints in z]" << endl
#endif
                    << "--level [#level (<= 1)]" << endl
                    << "--maxgrid [max. grid]" << endl
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
             const std::string dir)
{
    double time = 0.0;
    
    for (unsigned int i = 0; i < rho.size(); ++i)
        rho[i]->mult(/*Physics::epsilon_0 / */scalefactor, 0, 1);
    
    writePlotFile(dir, rho, phi, efield, rr, geom, time, scalefactor);
}

void doSolve(const Array<BoxArray>& ba,
             const Array<DistributionMapping>& dmap,
             const Array<Geometry>& geom,
             container_t& rhs,
             container_t& phi,
             container_t& efield,
             const Array<int>& rr,
             Inform& msg,
             const param_t& params, bool trilinos)
{
    static IpplTimings::TimerRef solvTimer = IpplTimings::getTimer("solve");
    
    rhs.resize(params.nLevels + 1);
    phi.resize(params.nLevels + 1);
    efield.resize(params.nLevels + 1);
    

    // Define the density on level 0 from all particles at all levels                                                                                                                                            
    int base_level   = 0;
    int finest_level = params.nLevels;
    
    for (uint lev = 0; lev < params.nLevels + 1; ++lev) {
        //                                                                       # component # ghost cells                                                                                                                                          
        rhs[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev],    1          , 0));
        phi[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev],    1          , 1));
        efield[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev], AMREX_SPACEDIM, 1));
        
        rhs[lev]->setVal(-1.0 - lev);
        phi[lev]->setVal(0.0, 1);
        efield[lev]->setVal(0.0, 1);
    }
    
    // normalize each level
    double l0norm = rhs[finest_level]->norm0(0);
    for (int i = 0; i <= finest_level; ++i) {
        rhs[i]->mult(1.0 / l0norm, 0, 1);
    }
    
    // solve
    
    IpplTimings::startTimer(solvTimer);
    
    if ( trilinos ) {
        AmrMultiGrid sol1;
        sol1.solve(rhs,            // [V m]
                   phi,            // [V m^3]
                   efield,       // [V m^2]
                   geom,
                   base_level,
                   finest_level);
    } else {
        Real offset = 0.;
        Solver sol2;
        sol2.solve_for_accel(rhs,
                             phi,
                             efield,
                             geom,
                             base_level,
                             finest_level,
                             offset);
    }
    
    IpplTimings::stopTimer(solvTimer);
    
    
    // undo normalization
    for (int i = 0; i <= finest_level; ++i) {
        rhs[i]->mult(l0norm, 0, 1);
        phi[i]->mult(l0norm, 0, 1);
    }
}


void doAMReX(const param_t& params, Inform& msg)
{
    int nlevs = params.nLevels + 1;
    
    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    IntVect domain_lo(D_DECL(0 , 0, 0)); 
    
    IntVect domain_hi(D_DECL(params.nr[0] - 1, params.nr[1] - 1, params.nr[2]-1)); 
    const Box domain(domain_lo, domain_hi);

    // Define the refinement ratio
    Array<int> rr(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++)
        rr[lev-1] = 2;

    // This says we are using Cartesian coordinates
    int coord = 0;

    // This sets the boundary conditions to Dirichlet
    int is_per[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++) 
        is_per[i] = 0; 

    // This defines a Geometry object which is useful for writing the plotfiles  
    Array<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, coord, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, coord, is_per);
    }

    Array<BoxArray> ba(nlevs);
    ba[0].define(domain);
    
    Array<DistributionMapping> dmap(nlevs);
    
    // Now we make the refined level be the center eighth of the domain
    if ( nlevs > 1 &&
         IntVect(D_DECL(params.nr[0], params.nr[1], params.nr[2])) == IntVect(D_DECL(8, 8, 8)) )
    {
//         int n_fine = params.nr[0]*rr[0];
//         IntVect refined_lo(D_DECL(n_fine/4,n_fine/4,n_fine/4)); 
//         IntVect refined_hi(D_DECL(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1));

        // Build a box for the level 1 domain
//         Box refined_patch(refined_lo, refined_hi);
        
        
        BoxList bl;
//         Box b1(IntVect(D_DECL(0, 4, 4)), IntVect(D_DECL(3, 7, 7)));
        
//         bl.push_back(b1);
        
        Box b2(IntVect(D_DECL(4, 4, 4)), IntVect(D_DECL(11, 11, 11)));
        
        bl.push_back(b2);
        
//         Box b3(IntVect(D_DECL(14, 6, 6)), IntVect(D_DECL(15, 9, 9)));
        
//         bl.push_back(b3);
        
        
        ba[1].define(bl);//define(refined_patch);
    } else if ( nlevs > 1 &&
                IntVect(D_DECL(params.nr[0], params.nr[1], params.nr[2])) == IntVect(D_DECL(32, 32, 32)) )
    {
        BoxList bl;
//         Box b1(IntVect(D_DECL(0, 4, 4)), IntVect(D_DECL(3, 7, 7)));
        
//         bl.push_back(b1);
        
        Box b2(IntVect(D_DECL(16, 16, 16)), IntVect(D_DECL(47, 47, 47)));
        
        bl.push_back(b2);
        
//         Box b3(IntVect(D_DECL(14, 6, 6)), IntVect(D_DECL(15, 9, 9)));
        
//         bl.push_back(b3);
        ba[1].define(bl);
    }
    
    
    if ( nlevs > 2 &&
         IntVect(D_DECL(params.nr[0], params.nr[1], params.nr[2])) == IntVect(D_DECL(32, 32, 32)) )
    {
        BoxList bl;
        
        Box b(IntVect(D_DECL(48, 48, 48)), IntVect(D_DECL(79, 79, 79)));
        
        bl.push_back(b);
        
        ba[2].define(bl);
    }
    
    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(params.maxBoxSize);
        dmap[lev].define(ba[lev]);
    }
    
    // ========================================================================
    // 3. multi-level redistribute
    // ========================================================================
    
    container_t rhs;
    container_t phi1, phi2;
    container_t efield1, efield2;
    
    // trilinos
    doSolve(ba, dmap, geom, rhs, phi1, efield1, rr, msg, params, true);
    
    for (int i = 0; i < nlevs; ++i) {
        msg << "Max. potential level " << i << ": "<< phi1[i]->max(0) << endl
            << "Min. potential level " << i << ": " << phi1[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << efield1[i]->max(0) << endl
            << "Min. ex-field level " << i << ": " << efield1[i]->min(0) << endl;
    }
    
    if ( params.isWriteYt )
        writeYt(rhs, phi1, efield1, geom, rr, 1.0, "Trilinos");
    
    // no trilinos
    doSolve(ba, dmap, geom, rhs, phi2, efield2, rr, msg, params, false);
    
    if ( params.isWriteYt )
        writeYt(rhs, phi2, efield2, geom, rr, 1.0, "FMGMultigrid");
    
    
    for (int i = 0; i < nlevs; ++i) {
        msg << "Max. potential level " << i << ": "<< phi2[i]->max(0) << endl
            << "Min. potential level " << i << ": " << phi2[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << efield2[i]->max(0) << endl
            << "Min. ex-field level " << i << ": " << efield2[i]->min(0) << endl;
    }
    
    
    double diff = 0.0;
    for (int i = 0; i < nlevs; ++i) {
        phi1[i]->minus(*phi2[i], 0, 1, 0);
        efield1[i]->minus(*efield2[i], 0, 3, 0);
    }
    
    if ( params.isWriteYt )
        writeYt(rhs, phi1, efield1, geom, rr, 1.0, "yt-testSolverComparison");
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
        
        msg << "Particle test running with" << endl
            << "- grid      = " << params.nr << endl
            << "- #level    = " << params.nLevels << endl
            << "- max. grid = " << params.maxBoxSize << endl;
        
        doAMReX(params, msg);
        
    } catch(const std::exception& ex) {
        msg << ex.what() << endl;
    }
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();

    return 0;
}