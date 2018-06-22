/*!
 * @file testWithFine.cpp
 * @author Matthias Frey
 * 
 * @details Test the application of the "with fine" Laplacian with
 * the proper fine and coarse boundary
 * 
 */
#include <iostream>

#include "Ippl.h"

#include <vector>

#include <random>

#include <getopt.h>

#include "build.h"

struct param_t {
    int nr[AMREX_SPACEDIM];
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


void fill(MultiFab& mf) {
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box&   bx  = mfi.validbox();
        amrex::FArrayBox&   fab = mf[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    IntVect iv(D_DECL(i, j, k));
                    
                    fab(iv) = i + j + 1.0;
#if AMREX_SPACEDIM == 3
                }
#endif
            }
        }
    }
}

void fill(MultiFab& mf, int seed) {
    
    std::mt19937_64 mt(seed);
    std::uniform_int_distribution<> dist(1, 10);
    
    
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box&   bx  = mfi.validbox();
        amrex::FArrayBox&   fab = mf[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    IntVect iv(D_DECL(i, j, k));
                    
                    fab(iv) = dist(mt);
#if AMREX_SPACEDIM == 3
                }
#endif
            }
        }
    }
}


void test(const param_t& params)
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
    Vector<int> rr(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++)
        rr[lev-1] = 2;

    // This says we are using Cartesian coordinates
    int coord = 0;

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++) 
        is_per[i] = 0; 

    // This defines a Geometry object which is useful for writing the plotfiles  
    Vector<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, coord, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, coord, is_per);
    }

    Vector<BoxArray> ba(nlevs);
    ba[0].define(domain);
    
    Vector<DistributionMapping> dmap(nlevs);
    
    // Now we make the refined level be the center eighth of the domain
    if (nlevs > 1) {
        BoxList bl;
        Box b1(IntVect(D_DECL(0, 4, 4)), IntVect(D_DECL(3, 7, 7)));
        
        bl.push_back(b1);
        
        Box b2(IntVect(D_DECL(4, 4, 4)), IntVect(D_DECL(11, 11, 11)));
        
        bl.push_back(b2);
        
        Box b3(IntVect(D_DECL(14, 6, 6)), IntVect(D_DECL(15, 9, 9)));
        
        bl.push_back(b3);
        
        
        ba[1].define(bl);//define(refined_patch);
    }
    
    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(params.maxBoxSize);
        
        std::cout << ba[lev] << std::endl;
        
        dmap[lev].define(ba[lev]);
    }
    
    
    
    
    Epetra_MpiComm epetra_comm = Ippl::getComm();
    
    std::vector<Teuchos::RCP<Epetra_Map> > maps(nlevs);
    
    
    IntVect rv(D_DECL(2, 2, 2));
    
    for (int i = 0; i < nlevs; ++i)
        buildMap(maps[i], ba[i], dmap[i],  geom[i], epetra_comm, i);
    
    
    container_t rhs(nlevs);
    container_t phi(nlevs);
    container_t efield(nlevs);
    for (int lev = 0; lev < nlevs; ++lev) {
        //                                                                       # component # ghost cells                                                                                                                                          
        rhs[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev],    1          , 0));
        rhs[lev]->setVal(0.0);
        
        // not used (only for plotting)
        phi[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev],    1          , 1));
        efield[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev], AMREX_SPACEDIM, 1));
            
        phi[lev]->setVal(0.0, 1);
        efield[lev]->setVal(0.0, 1);
    }
    
    
    fill(*phi[0], 42);
    fill(*phi[1], 43);
    
    Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::null;
    Teuchos::RCP<Epetra_CrsMatrix> B = Teuchos::null;
    
    
    // Poisson matrix for level 1
    buildSpecialPoissonMatrix(A, maps, ba, dmap, geom, rv, epetra_comm, 0);
    
    // Boundary matrix for Dirichlet BC (level 1 to level 0)
    buildFineBoundaryMatrix(B, maps, ba, dmap, geom, rv, epetra_comm, 0);
    
    // vector for level 0 and level 1
    Teuchos::RCP<Epetra_Vector> crse = Teuchos::null;
    buildVector(crse, maps[0], 0.0);
    
    amrex2trilinos(*phi[0], crse, maps[0], geom, 0);
    
    Teuchos::RCP<Epetra_Vector> fine = Teuchos::null;
    buildVector(fine, maps[1], 0.0);
    
    amrex2trilinos(*phi[1], fine, maps[1], geom, 1);
    
    // apply boundary
    Teuchos::RCP<Epetra_Vector> tmp = Teuchos::null;
    buildVector(tmp, maps[0], 0.0);
    
    // tmp = B * fine
    B->Apply(*fine, *tmp);
    
    std::cout << *fine << std::endl;
    
    std::cout << *crse << std::endl;
    
    
    // apply Poisson
    Teuchos::RCP<Epetra_Vector> result = Teuchos::null;
    buildVector(result, maps[0], 0.0);
    A->Apply(*crse, *result);
    
    std::cout << *tmp << std::endl;
    
    std::cout << *result << std::endl;
    
    // add boundary
    result->Update(1.0, *tmp, 1.0);
    
    std::cout << *result << std::endl;
    
    
    trilinos2amrex(*rhs[0], result);
    
    if ( params.isWriteYt )
        writeYt(rhs, phi, efield, geom, rr, 1.0, "testWithFine");
}


int main(int argc, char* argv[])
{
    Ippl ippl(argc, argv);
    
    Inform msg(argv[0]);
    
    param_t params;
    
    amrex::Initialize(argc,argv, false);
    
    try {
        if ( !parseProgOptions(argc, argv, params, msg) && !params.isHelp )
            throw std::runtime_error("\033[1;31mError: Check the program options.\033[0m");
        else if ( params.isHelp )
            return 0;
        
        msg << "Particle test running with" << endl
            << "- grid      = " << params.nr[0] << " " << params.nr[1]
#if AMREX_SPACEDIM == 3
            << " " << params.nr[2]
#endif
            << endl
            << "- #level    = " << params.nLevels << endl
            << "- max. grid = " << params.maxBoxSize << endl;
        
        test(params);
        
    } catch(const std::exception& ex) {
        msg << ex.what() << endl;
    }
    
    amrex::Finalize();
}
