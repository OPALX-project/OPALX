/*!
 * @file testPoissonMatrix.cpp
 * @author Matthias Frey
 * 
 */
#include <iostream>

#include "Ippl.h"

#include <vector>

#include <getopt.h>

#include "tools.h"

#include "EpetraExt_RowMatrixOut.h"

struct param_t {
    int nr[BL_SPACEDIM];
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
    
    int required = 2 +  BL_SPACEDIM;
    
    while ( true ) {
        static struct option long_options[] = {
            { "gridx",          required_argument, 0, 'x' },
            { "gridy",          required_argument, 0, 'y' },
#if BL_SPACEDIM == 3
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
#if BL_SPACEDIM == 3
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
#if BL_SPACEDIM == 3
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
#if BL_SPACEDIM == 3
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


void buildPoissonMatrix(Teuchos::RCP<Epetra_CrsMatrix>& A,
                        const std::vector<Teuchos::RCP<Epetra_Map> >& maps,
                        const Array<BoxArray>& grids, const Array<DistributionMapping>& dmap,
                        const Array<Geometry>& geom,
                        Epetra_MpiComm& comm, int level)
{
    /*
     * Laplacian of "no fine"
     */
    
    /*
     * 1D not supported
     * 2D --> 5 elements per row
     * 3D --> 7 elements per row
     */
    int nEntries = (BL_SPACEDIM << 1) + 1 /* plus boundaries */ + 10 /*FIXME*/;
    
    std::cout << "nEntries = " << nEntries << std::endl;
    
    const Epetra_Map& RowMap = *maps[level];
    
    A = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy, RowMap, nEntries) );
    
    std::vector<int> indices; //(nEntries);
    std::vector<double> values; //(nEntries);
    
    const double* dx = geom[level].CellSize();
    
    std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > mask(
        new amrex::FabArray<amrex::BaseFab<int> >(grids[level], dmap[level], 1, 1)
    );
    
    mask->BuildMask(geom[level].Domain(), geom[level].periodicity(), -1, 1, 2, 0);
    
    int nr[BL_SPACEDIM];
    for (int j = 0; j < BL_SPACEDIM; ++j)
        nr[j] = geom[level].Domain().length(j);
    
    for (amrex::MFIter mfi(*mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mask)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    int numEntries = 0;
                    indices.clear();
                    values.clear();
                    IntVect iv(D_DECL(i, j, k));
                    int globidx = serialize(iv, &nr[0]);
                    
                    /*
                     * check neighbours in all directions (Laplacian stencil --> cross)
                     */
                    for (int d = 0; d < BL_SPACEDIM; ++d) {
                        for (int shift = -1; shift <= 1; shift += 2) {
                            IntVect biv = iv;                        
                            biv[d] += shift;
                            
                            switch ( mfab(biv) )
                            {
                                case -1:
                                    // covered --> interior cell
                                case 0:
                                {
                                    indices.push_back( serialize(biv, &nr[0]) );
                                    values.push_back( -1.0 / ( dx[d] * dx[d] ) );
                                    break;
                                }
                                case 1:
                                    // boundary cell
                                    // only level > 0 have this kind of boundary
                                    
                                    /* Dirichlet boundary conditions from coarser level.
                                     * These are treated by the boundary matrix.
                                     */
                                    break;
                                case 2:
                                {
                                    // physical boundary cell
                                    double value = -1.0 / ( dx[d] * dx[d] );
                                    applyBoundary(biv,
                                                  indices,
                                                  values,
                                                  numEntries,
                                                  value,
                                                  &nr[0]);
                                    break;
                                }
                                default:
                                    break;
                            }
                        }
                    }
                    
                    // check center
                    if ( mfab(iv) == 0 ) {
                        indices.push_back( globidx );
                        values.push_back( 2.0 / ( dx[0] * dx[0] ) +
                                          2.0 / ( dx[1] * dx[1] )
#if BL_SPACEDIM == 3
                                          + 2.0 / ( dx[2] * dx[2] )
#endif
                        );
                        ++numEntries;
                    }
                    
//                     for (uint ii = 0; ii < indices.size(); ++ii)
//                         std::cout << indices[ii] << " ";
//                     std::cout << std::endl;
//                     
//                     for (uint ii = 0; ii < values.size(); ++ii)
//                         std::cout << values[ii] << " ";
//                     std::cout << std::endl;
                    
                    unique(indices, values, numEntries);
                    
//                     for (uint ii = 0; ii < indices.size(); ++ii)
//                         std::cout << indices[ii] << " ";
//                     std::cout << std::endl;
//                     
//                     for (uint ii = 0; ii < values.size(); ++ii)
//                         std::cout << values[ii] << " ";
//                     std::cout << std::endl;
                    
                    
                    int error = A->InsertGlobalValues(globidx,
                                                      numEntries,
                                                      &values[0],
                                                      &indices[0]);
                    
                    if ( error != 0 )
                        throw std::runtime_error("Error in filling the Poisson matrix for level "
                                                 + std::to_string(level) + "!");
                
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    int error = A->FillComplete(true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing Poisson matrix for level "
                                 + std::to_string(level) + "!");

    std::cout << *A << std::endl;
    
    
    EpetraExt::RowMatrixToMatlabFile("poisson_matrix.txt", *A);
}


void fill(MultiFab& mf) {
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box&   bx  = mfi.validbox();
        amrex::FArrayBox&   fab = mf[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    IntVect iv(D_DECL(i, j, k));
                    
                    fab(iv) = i + j + 1.0;
#if BL_SPACEDIM == 3
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

    // This sets the boundary conditions to be doubly or triply periodic
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
    
    
    
    
    for (int i = 0; i < nlevs; ++i)
        buildMap(maps[i], ba[i], dmap[i],  geom[i], epetra_comm, i);
    
    
    
    
    Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::null;
    
    buildPoissonMatrix(A, maps, ba, dmap, geom, epetra_comm, 0);
    
//     buildPoissonMatrix(A, maps, ba, dmap, geom, epetra_comm, 1);
    
//     container_t rhs(nlevs);
//     container_t phi(nlevs);
//     container_t efield(nlevs);
//     for (int lev = 0; lev < nlevs; ++lev) {
//         //                                                                       # component # ghost cells                                                                                                                                          
//         rhs[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev],    1          , 0));
//         rhs[lev]->setVal(0.0);
//         
//         // not used (only for plotting)
//         if ( params.isWriteYt ) {
//             phi[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev],    1          , 1));
//             efield[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev], BL_SPACEDIM, 1));
//             
//             phi[lev]->setVal(0.0, 1);
//             efield[lev]->setVal(0.0, 1);
//         }
//     }
//     
//     // fill coarsest level
//     fill(*rhs[0]);
//     
//     // coarse
//     Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp( new Epetra_Vector(*maps[0], false));
//     
//     amrex2trilinos(*rhs[0], x, maps[0], geom, 0);
//     
//     // fine
//     Teuchos::RCP<Epetra_Vector> y = Teuchos::rcp( new Epetra_Vector(*maps[1], false));
//     
//     amrex2trilinos(*rhs[1], y, maps[1], geom, 1);
//     
//     // y = I * x
//     I->Multiply(false, *x, *y);
//     
//     trilinos2amrex(*rhs[1], y);
//     
//     if ( params.isWriteYt )
//         writeYt(rhs, phi, efield, geom, rr, 1.0, "testPoissonMatrix");
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
#if BL_SPACEDIM == 3
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
