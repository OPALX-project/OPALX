/*!
 * @file testBoundaryMatrix.cpp
 * @author Matthias Frey
 * 
 */
#include <iostream>

#include "Ippl.h"

#include <vector>

#include <getopt.h>

#include <Epetra_MpiComm.h>

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

void checkBoundary(Teuchos::RCP<Epetra_CrsMatrix>& B,
                   const Array<Geometry>& geom,
                   int level,
                   const amrex::BaseFab<int>& mfab,
                   const IntVect& lo,
                   const IntVect& hi)
{
    int nr[BL_SPACEDIM];
    for (int j = 0; j < BL_SPACEDIM; ++j)
        nr[j] = geom[level].Domain().length(j);
    
    int cnr[BL_SPACEDIM];
    if ( level > 0 ) {
        for (int j = 0; j < BL_SPACEDIM; ++j)
            cnr[j] = geom[level-1].Domain().length(j);
    }
        
    
    const double* dx = geom[level].CellSize();
    
    std::vector<int> indices; //(nNeighbours);
    std::vector<double> values; //(nNeighbours);
    
    // helper function for boundary
    auto check = [&](const amrex::BaseFab<int>& mfab,
                     const IntVect& iv,
                     int& numEntries,
                     std::vector<int>& indices,
                     std::vector<double>& values,
                     int dir)
    {
        switch ( mfab(iv) )
        {
            case -1:
                // covered (nothing to do here)
                break;
            case 0:
                // interior (nothing to do here)
                break;
            case 1:
            {
                // internal boundary (only level > 0 have this kind of boundary)
                
                int nn = numEntries;
                // we need boundary + indices from coarser level
                stencil(iv, indices, values, numEntries, &cnr[0]);
                
                // we need normlization by mesh size squared
                for (int n = nn; n < numEntries; ++n)
                    values[n] *= 1.0 / ( dx[dir] * dx[dir] );
                
                break;
            }
            case 2:
            {
                // physical boundary --> handled in Poisson matrix
                break;
            }
            default:
                throw std::runtime_error("Error in mask value");
                break;
        }
    };
    
    for (int i = lo[0]; i <= hi[0]; ++i) {
        for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
            for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                // last interior cell
                IntVect iv(D_DECL(i, j, k));
                
                int numEntries = 0;
                indices.clear();
                values.clear();
                
                // we need the global index of the interior cell (the center of the Laplacian stencil (i, j, k)
                int globidx = serialize(iv, &nr[0]);
                
                /*
                 * check all directions (Laplacian stencil --> cross)
                 */
                for (int d = 0; d < BL_SPACEDIM; ++d) {
                    for (int shift = -1; shift <= 1; shift += 2) {
                        IntVect biv = iv;                        
                        biv[d] += shift;
                        check(mfab, biv,
                              numEntries,
                              indices,
                              values,
                              d);
                    }
                }
                
//                 std::cout << iv << " " << globidx << std::endl;
//                 
//                 for (uint ii = 0; ii < indices.size(); ++ii)
//                     std::cout << indices[ii] << " ";
//                 std::cout << std::endl;
//                 
//                 for (uint ii = 0; ii < values.size(); ++ii)
//                     std::cout << values[ii] << " ";
//                 std::cout << std::endl;
                
                /*
                 * In some cases indices occur several times, e.g. for corner points
                 * at the physical boundary --> sum them up
                 */
                unique(indices, values, numEntries);
                
//                 for (uint ii = 0; ii < indices.size(); ++ii)
//                     std::cout << indices[ii] << " ";
//                 std::cout << std::endl;
//                 
//                 for (uint ii = 0; ii < values.size(); ++ii)
//                     std::cout << values[ii] << " ";
//                 std::cout << std::endl; std::cin.get();
                
                int error = B->InsertGlobalValues(globidx,
                                                  numEntries,
                                                  &values[0],
                                                  &indices[0]);
                
                if ( error != 0 ) {
                    // if e.g. nNeighbours < numEntries --> error
                    throw std::runtime_error("Error in filling the boundary matrix for level " +
                                             std::to_string(level) + "!");
                }
#if BL_SPACEDIM == 3
            }
#endif
        }
    }
}


void buildCrseBoundaryMatrix(Teuchos::RCP<Epetra_CrsMatrix>& Bcrse,
                             const std::vector<Teuchos::RCP<Epetra_Map> >& maps,
                             const Array<BoxArray>& grids, const Array<DistributionMapping>& dmap,
                             const Array<Geometry>& geom,
                             const IntVect& rr, Epetra_MpiComm& comm, int level)
{
    if ( level == 0 )
        return;
    
    const Epetra_Map& RowMap = *maps[level];
    const Epetra_Map& ColMap = *maps[level-1];
    
    int nNeighbours = (2 << (BL_SPACEDIM -1 )) /*FIXME interpolation stencil indices*/ + 10;
    
    Bcrse = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy,
                                               RowMap, nNeighbours, false) );
    
    std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > mask(
        new amrex::FabArray<amrex::BaseFab<int> >(grids[level], dmap[level], 1, 1)
    );
    
    mask->BuildMask(geom[level].Domain(), geom[level].periodicity(), -1, 1, 2, 0);
    
    for (amrex::MFIter mfi(*mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box&    bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mask)[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        /*
         * left boundary
         */
        IntVect lower(D_DECL(lo[0], lo[1], lo[2]));
        IntVect upper(D_DECL(lo[0], hi[1], hi[2]));
        
        checkBoundary(Bcrse, geom, level, mfab, lower, upper);
        
        /*
         * right boundary
         */
        lower = IntVect(D_DECL(hi[0], lo[1], lo[2]));
        upper = IntVect(D_DECL(hi[0], hi[1], hi[2]));
        
        checkBoundary(Bcrse, geom, level, mfab, lower, upper);
        
        
        /*
         * lower boundary (without left and right last cell!)
         */
        lower = IntVect(D_DECL(lo[0]+1, lo[1], lo[2]));
        upper = IntVect(D_DECL(hi[0]-1, lo[1], hi[2]));
        
        checkBoundary(Bcrse, geom, level, mfab, lower, upper);
        
        /*
         * upper boundary (without left and right last cell!)
         */
        lower = IntVect(D_DECL(lo[0]+1, hi[1], lo[2]));
        upper = IntVect(D_DECL(hi[0]-1, hi[1], hi[2]));
        
        checkBoundary(Bcrse, geom, level, mfab, lower, upper);
        
#if BL_SPACEDIM == 3
        /*
         * front boundary (without left, right, upper and lower last cell!)
         */
        lower = IntVect(D_DECL(lo[0]+1, lo[1]+1, lo[2]));
        upper = IntVect(D_DECL(hi[0]-1, hi[1]-1, lo[2]));
        
        checkBoundary(Bcrse, geom, level, mfab, lower, upper);
        
        /*
         * back boundary (without left, right, upper and lower last cell!)
         */
        lower = IntVect(D_DECL(lo[0]+1, lo[1]+1, hi[2]));
        upper = IntVect(D_DECL(hi[0]-1, hi[1]-1, hi[2]));
        
        checkBoundary(Bcrse, geom, level, mfab, lower, upper);
#endif
        
    }
    
    int error = Bcrse->FillComplete(ColMap, RowMap, true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing the crse boundary matrix for level " +
                                 std::to_string(level) + "!");
    
    std::cout << "Done." << std::endl;
    
    std::cout << *Bcrse << std::endl;
    
    EpetraExt::RowMatrixToMatlabFile("crse_boundary_matrix.txt", *Bcrse);
}


void buildFineBoundaryMatrix(Teuchos::RCP<Epetra_CrsMatrix>& Bfine,
                             const std::vector<Teuchos::RCP<Epetra_Map> >& maps,
                             const Array<BoxArray>& grids, const Array<DistributionMapping>& dmap,
                             const Array<Geometry>& geom,
                             const IntVect& rr, Epetra_MpiComm& comm, int level)
{
    if ( level == 1 /*finest*/ )
        return;
    
    // find all cells with refinement
    amrex::BoxArray crse_fine_ba = grids[level];
    crse_fine_ba.refine(rr);
    crse_fine_ba = amrex::intersect(grids[level+1], crse_fine_ba);
    crse_fine_ba.coarsen(rr);
    
    const Epetra_Map& RowMap = *maps[level];
    const Epetra_Map& ColMap = *maps[level+1];
    
    int nNeighbours = 4 /*#interfaces*/ * rr[0] * rr[1] /*of refined cell*/;
#if BL_SPACEDIM == 3
    nNeighbours = 6 /*#interfaces*/ * rr[0] * rr[1] * rr[2] /*of refined cell*/;
#endif
    
    Bfine = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy,
                                               RowMap, /*ColMap, */nNeighbours, false) );
    
    std::vector<int> indices; //(nEntries);
    std::vector<double> values; //(nEntries);
    
    
    std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > mask(
        new amrex::FabArray<amrex::BaseFab<int> >(grids[level], dmap[level], 1, 1)
    );
    
    mask->BuildMask(geom[level].Domain(), geom[level].periodicity(), -1, 1, 2, 0);
    
    int cnr[BL_SPACEDIM];
    int fnr[BL_SPACEDIM];
    for (int j = 0; j < BL_SPACEDIM; ++j) {
        cnr[j] = geom[level].Domain().length(j);
        fnr[j] = geom[level+1].Domain().length(j);
    }
    
    const double* cdx = geom[level].CellSize();
    const double* fdx = geom[level+1].CellSize();
    
    /* Find all coarse cells that are at the crse-fine interace but are
     * not refined.
     * Put them into a map (to avoid duplicates, e.g. due to corners).
     * Finally, iterate through the list of cells
     */
    
    // std::list<std::pair<int, int> > contains the shift and direction to come to the covered cell
    std::map<IntVect, std::list<std::pair<int, int> > > cells;
    
    for (amrex::MFIter mfi(grids[level], dmap[level], false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mask)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    IntVect iv(D_DECL(i, j, k));
                    
                    if ( !crse_fine_ba.contains(iv) && mfab(iv) != 2 /*not physical bc*/ ) {
                        /*
                         * iv is a coarse cell that got not refined
                         * 
                         * --> check all neighbours to see if at crse-fine
                         * interface
                         */
                        for (int d = 0; d < BL_SPACEDIM; ++d) {
                            for (int shift = -1; shift <= 1; shift += 2) {
                                // neighbour
                                IntVect covered = iv;
                                covered[d] += shift;
                                
                                if ( crse_fine_ba.contains(covered) &&
                                     !isBoundary(covered, &cnr[0]) /*not physical bc*/ )
                                {
                                    // neighbour is covered by fine cells
                                    
                                    // insert if not yet exists
                                    cells[iv].push_back(std::make_pair(shift, d));
                                }
                            }
                        }
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    
    auto fill = [&](std::vector<int>& indices,
                         std::vector<double>& values,
                         int& numEntries,
                         D_DECL(int ii, int jj, int kk),
                         int* begin, int* end, int d, const IntVect& iv, int shift, int sign)
    {
        for (int iref = ii - begin[0]; iref <= ii + end[0]; ++iref) {
            
            sign *= ( d == 0 ) ? -1.0 : 1.0;
            
            for (int jref = jj - begin[1]; jref <= jj + end[1]; ++jref) {
                
                sign *= ( d == 1 ) ? -1.0 : 1.0;
#if BL_SPACEDIM == 3
                for (int kref = kk - begin[2]; kref <= kk + end[2]; ++kref) {
#endif
                    /* Since all fine cells on the not-refined cell are
                     * outside of the "domain" --> we need to interpolate
                     * using open boundary condition.
                     */
                    sign *= ( d == 2 ) ? -1.0 : 1.0;
                    
                    IntVect riv(D_DECL(iref, jref, kref));
                    
                    if ( riv[d] / rr[d] == iv[d] ) {
                        /* interpolate
                         */
                        
                        // undo shift
                        riv[d] += shift;
                        
                        double value = sign / ( cdx[d] * fdx[d] );
                        
                        indices.push_back( serialize(riv, &fnr[0]) );
                        values.push_back( 2.0 * value );
                        ++numEntries;
                        
                        riv[d] += shift;
                        
                        indices.push_back( serialize(riv, &fnr[0]) );
                        values.push_back( - value );
                        ++numEntries;
                    } else {
                        indices.push_back( serialize(riv, &fnr[0]) );
                        values.push_back( sign / ( cdx[d] * fdx[d] ) );
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    };
    
    for (std::map<IntVect, std::list<std::pair<int, int> > >::iterator it = cells.begin(); it != cells.end(); ++it) {
        
        // not covered coarse cell at crse-fine interface
        IntVect iv = it->first;
        int globidx = serialize(iv, &cnr[0]);
        
        int numEntries = 0;
        indices.clear();
        values.clear();
        
        while ( !it->second.empty() ) {
            /*
             * "shift" is the amount of is a coarse cell that got refined
             * "d" is the direction to shift
             * 
             * --> check all covered neighbour cells
             */
            int shift = it->second.front().first;
            int d     = it->second.front().second;
            it->second.pop_front();
            
            
            /* we need to iterate over correct fine cells. It depends
             * on the orientation of the interface
             */
            int begin[BL_SPACEDIM] = { D_DECL( int(d == 0), int(d == 1), int(d == 2) ) };
            int end[BL_SPACEDIM]   = { D_DECL( int(d != 0), int(d != 1), int(d != 2) ) };
            
            
            // neighbour
            IntVect covered = iv;
            covered[d] += shift;
                    
            /*
             * neighbour cell got refined but is not on physical boundary
             * --> we are at a crse-fine interface
             * 
             * we need now to find out which fine cells
             * are required to satisfy the flux matching
             * condition
             */
            
            switch ( shift ) {
                case -1:
                {
                    // --> interface is on the lower face
                    int ii = iv[0] * rr[0];
                    int jj = iv[1] * rr[1];
#if BL_SPACEDIM == 3
                    int kk = iv[2] * rr[2];
#endif
                    // iterate over all fine cells at the interface
                    // start with lower cells --> cover coarse neighbour
                    // cell
                    fill(indices, values, numEntries, D_DECL(ii, jj, kk), &begin[0], &end[0], d, iv, shift, 1.0);
                    break;
                }
                case 1:
                default:
                {
                    // --> interface is on the upper face
                    int ii = covered[0] * rr[0];
                    int jj = covered[1] * rr[1];
#if BL_SPACEDIM == 3
                    int kk = covered[2] * rr[2];
#endif
                    fill(indices, values, numEntries, D_DECL(ii, jj, kk), &begin[0], &end[0], d, iv, shift, -1.0);
                    break;
                }
            }
        }
        
        unique(indices, values, numEntries);
        
        std::cout << globidx << std::endl;
        for (uint ii = 0; ii < indices.size(); ++ii)
            std::cout << indices[ii] << " ";
        std::cout << std::endl;
        
        int error = Bfine->InsertGlobalValues(globidx,
                                              numEntries,
                                              &values[0],
                                              &indices[0]);
        
        if ( error != 0 ) {
            // if e.g. nNeighbours < numEntries --> error
            throw std::runtime_error("Error in filling the boundary matrix for level " +
                                     std::to_string(level) + "!");
        }
    }
    
    int error = Bfine->FillComplete(ColMap, RowMap, true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing the crse boundary matrix for level " +
                                 std::to_string(level) + "!");
    
    std::cout << "Done." << std::endl;
    
    std::cout << *Bfine << std::endl;
    
    EpetraExt::RowMatrixToMatlabFile("fine_boundary_matrix.txt", *Bfine);
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
    
    
    
    
    Teuchos::RCP<Epetra_CrsMatrix> Bcrse = Teuchos::null; // coarse to fine (Dirichlet) 
    Teuchos::RCP<Epetra_CrsMatrix> Bfine = Teuchos::null; // fine to coarse (flux matching)
    
    IntVect rv(D_DECL(2, 2, 2));
    
//     buildCrseBoundaryMatrix(Bcrse, maps, ba, dmap, geom, rv, epetra_comm, 0);
    buildFineBoundaryMatrix(Bfine, maps, ba, dmap, geom, rv, epetra_comm, 0);
    
//     std::cout << "LEVEL 1" << std::endl; std::cin.get();
    
    buildCrseBoundaryMatrix(Bcrse, maps, ba, dmap, geom, rv, epetra_comm, 1);
//     buildFineBoundaryMatrix(Bfine, maps, ba, dmap, geom, rv, epetra_comm, 1);
    
    container_t rhs(nlevs);
    container_t phi(nlevs);
    container_t efield(nlevs);
    for (int lev = 0; lev < nlevs; ++lev) {
        //                                                                       # component # ghost cells                                                                                                                                          
        rhs[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev],    1          , 0));
        rhs[lev]->setVal(0.0);
        
        // not used (only for plotting)
        if ( params.isWriteYt ) {
            phi[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev],    1          , 1));
            efield[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev], BL_SPACEDIM, 1));
            
            phi[lev]->setVal(0.0, 1);
            efield[lev]->setVal(0.0, 1);
        }
    }
    
    // fill coarsest level
    fill(*rhs[0]);
    
    // coarse
    Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp( new Epetra_Vector(*maps[0], false));
    
    amrex2trilinos(*rhs[0], x, maps[0], geom, 0);
    
    // fine
    Teuchos::RCP<Epetra_Vector> y = Teuchos::rcp( new Epetra_Vector(*maps[1], false));
    
    amrex2trilinos(*rhs[1], y, maps[1], geom, 1);
    
    // y = B * x
    // fine to coarse
    Bfine->Multiply(false, *y, *x);
    
    // coarse to fine
    Bcrse->Multiply(false, *x, *y);
    
//     
//     trilinos2amrex(*rhs[1], y);
//     
//     if ( params.isWriteYt )
//         writeYt(rhs, phi, efield, geom, rr, 1.0);
}


int main(int argc, char* argv[])
{
    Ippl ippl(argc, argv);
    
    Inform msg("Interpolater");
    
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
