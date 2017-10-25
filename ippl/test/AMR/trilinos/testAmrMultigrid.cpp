#include "Ippl.h"

#include <AMReX_ParmParse.H>

#include "AmrMultiGrid.h"

#include "../writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"
#include <random>

#include <getopt.h>

#include <AMReX_MultiFab.H>

using namespace amrex;

typedef Array<std::unique_ptr<MultiFab> > container_t;

struct param_t {
    Vektor<size_t, BL_SPACEDIM> nr;
    size_t maxBoxSize;
    bool isWriteYt;
    bool isWriteCSV;
    bool isHelp;
    size_t nLevels;
    size_t smoothing;
};

bool parseProgOptions(int argc, char* argv[], param_t& params, Inform& msg) {
    /* Parsing Command Line Arguments
     * 
     * 26. June 2017
     * https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example
     */
    
    params.isWriteYt = false;
    params.isWriteCSV = false;
    params.isHelp = false;
    
    int c = 0;
    
    int cnt = 0;
    
    int required = 3 +  BL_SPACEDIM;
    
    while ( true ) {
        static struct option long_options[] = {
            { "gridx",          required_argument, 0, 'x' },
            { "gridy",          required_argument, 0, 'y' },
#if BL_SPACEDIM == 3
            { "gridz",          required_argument, 0, 'z' },
#endif
            { "smoothing",      required_argument, 0, 's' },
            { "level",          required_argument, 0, 'l' },
            { "maxgrid",        required_argument, 0, 'm' },
            { "writeYt",        no_argument,       0, 'w' },
            { "help",           no_argument,       0, 'h' },
            { "writeCSV",       no_argument,       0, 'v' },
            { 0,                0,                 0,  0  }
        };
        
        int option_index = 0;
        
        c = getopt_long(argc, argv,
#if BL_SPACEDIM == 3
                        "x:y:z:s:l:m:whv",
#else
                        "x:y:s:l:m:whv",
#endif
                        long_options, &option_index);
        
        if ( c == -1 )
            break;
        
        switch ( c ) {
            case 'x':
                params.nr[0] = std::atoi(optarg); ++cnt; break;
            case 'y':
                params.nr[1] = std::atoi(optarg); ++cnt; break;
#if BL_SPACEDIM == 3
            case 'z':
                params.nr[2] = std::atoi(optarg); ++cnt; break;
#endif
            case 's':
                params.smoothing = std::atoi(optarg); ++cnt; break;
            case 'l':
                params.nLevels = std::atoi(optarg); ++cnt; break;
            case 'm':
                params.maxBoxSize = std::atoi(optarg); ++cnt; break;
            case 'w':
                params.isWriteYt = true;
                break;
            case 'v':
                params.isWriteCSV = true;
                break;
            case 'h':
                msg << "Usage: " << argv[0]
                    << endl
                    << "--gridx [#gridpoints in x]" << endl
                    << "--gridy [#gridpoints in y]" << endl
#if BL_SPACEDIM == 3
                    << "--gridz [#gridpoints in z]" << endl
#endif
                    << "--smoothing [#steps]" << endl
                    << "--level [#level (<= 1)]" << endl
                    << "--maxgrid [max. grid]" << endl
                    << "--writeYt (optional)" << endl
                    << "--writeCSV (optional)" << endl;
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
//     // Immediate debug output:
//     // Output potential and e-field along axis
//     std::string outfile = "potential.grid";
//     std::ofstream out;
//     for (MFIter mfi(*phi[0]); mfi.isValid(); ++mfi) {
//         const Box& bx = mfi.validbox();
//         const FArrayBox& lhs = (*phi[0])[mfi];
//         
//         for (int proc = 0; proc < amrex::ParallelDescriptor::NProcs(); ++proc) {
//             if ( proc == amrex::ParallelDescriptor::MyProc() ) {
//                 
//                 if ( proc == 0 ) {
//                     out.open(outfile, std::ios::out);
//                     out << "$x$ [m], $\\Phi$ [V]" << std::endl;
//                 } else
//                     out.open(outfile, std::ios::app);
//                 
//                 int j = 0.5 * (bx.hiVect()[1] - bx.loVect()[1]);
//                 int k = 0.5 * (bx.hiVect()[2] - bx.loVect()[2]);
//                 
//                 for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
//                     IntVect ivec(i, j, k);
//                     // add one in order to have same convention as PartBunch::computeSelfField()
//                     out << lower + i * dx << ", " << lhs(ivec, 0)  << std::endl;
//                 }
//                 out.close();
//             }
//             amrex::ParallelDescriptor::Barrier();
//         }
//     }
//     
//     outfile = "efield.grid";
//     for (MFIter mfi(*efield[0]); mfi.isValid(); ++mfi) {
//         const Box& bx = mfi.validbox();
//         const FArrayBox& lhs = (*efield[0])[mfi];
//         
//         for (int proc = 0; proc < amrex::ParallelDescriptor::NProcs(); ++proc) {
//             if ( proc == amrex::ParallelDescriptor::MyProc() ) {
//                 
//                 if ( proc == 0 ) {
//                     out.open(outfile, std::ios::out);
//                     out << "$x$ [m], $E_x$ [V/m], $E_x$ [V/m], $E_x$ [V/m]" << std::endl;
//                 } else
//                     out.open(outfile, std::ios::app);
//                 
//                 int j = 0.5 * (bx.hiVect()[1] - bx.loVect()[1]);
//                 int k = 0.5 * (bx.hiVect()[2] - bx.loVect()[2]);
//                 
//                 for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
//                     IntVect ivec(i, j, k);
//                     // add one in order to have same convention as PartBunch::computeSelfField()
//                     out << lower + i * dx << ", " << lhs(ivec, 0) << ", "
//                         << lhs(ivec, 1) << ", " << lhs(ivec, 2) << std::endl;
//                 }
//                 out.close();
//             }
//             amrex::ParallelDescriptor::Barrier();
//         }
//     }
//     
}


void writeYt(container_t& rho,
             const container_t& phi,
             const container_t& efield,
             const Array<Geometry>& geom,
             const Array<int>& rr,
             const double& scalefactor)
{
    std::string dir = "yt-testAmrMultigrid";
    
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
             const param_t& params)
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
        efield[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev], BL_SPACEDIM, 1));
        
        rhs[lev]->setVal(-1.0 /*- lev*/);
        phi[lev]->setVal(0.0, 1);
        efield[lev]->setVal(0.0, 1);
    }
    
    // normalize each level
    double l0norm = rhs[finest_level]->norm0(0);
    for (int i = 0; i <= finest_level; ++i) {
        rhs[i]->mult(1.0 / l0norm, 0, 1);
    }
    
    // solve
    AmrMultiGrid sol(AmrMultiGrid::Boundary::DIRICHLET, AmrMultiGrid::Interpolater::PIECEWISE_CONST);
    
    sol.setNumberOfSmoothing(params.smoothing);
    
    IpplTimings::startTimer(solvTimer);
    
    sol.solve(rhs,            // [V m]
              phi,            // [V m^3]
              efield,       // [V m^2]
              geom,
              base_level,
              finest_level);
    
    IpplTimings::stopTimer(solvTimer);
    
    msg << "#iterations: " << sol.getNumIters() << endl;
    
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
    for (int n = 0; n < BL_SPACEDIM; n++) {
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
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) 
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
    container_t phi;
    container_t efield;
        
    doSolve(ba, dmap, geom, rhs, phi, efield, rr, msg, params);
//     
    for (int i = 0; i < nlevs; ++i) {
        msg << "Max. potential level " << i << ": "<< phi[i]->max(0) << endl
            << "Min. potential level " << i << ": " << phi[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(0) << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(0) << endl;
    }
//     
//     double fieldenergy = totalFieldEnergy(efield, rr);
//     
//     msg << "Total field energy: " << fieldenergy << endl;
//     
    if (params.isWriteCSV && Ippl::getNodes() == 1 && params.maxBoxSize == params.nr[0] )
        writeCSV(phi, efield, real_box.lo(0), geom[0].CellSize(0));
    
    if ( params.isWriteYt )
        writeYt(rhs, phi, efield, geom, rr, 1.0);
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
            << "- smoothing = " << params.smoothing << endl
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