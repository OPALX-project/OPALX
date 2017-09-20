/*!
 * @file testInterpolationMatrix.cpp
 * @author Matthias Frey
 * 
 */
#include <iostream>

#include "Ippl.h"

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_PlotFileUtil.H>

#include <vector>

#include <getopt.h>

#include "../writePlotFile.H"

#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#include <Teuchos_RCP.hpp>
// #include <Teuchos_ArrayRCP.hpp>

#include "EpetraExt_RowMatrixOut.h"

using namespace amrex;

typedef Array<std::unique_ptr<MultiFab> > container_t;

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

void writeYt(container_t& rho,
             const container_t& phi,
             const container_t& efield,
             const Array<Geometry>& geom,
             const Array<int>& rr,
             const double& scalefactor)
{
    std::string dir = "yt-testInterpolationMatrix";
    
    double time = 0.0;
    
    for (unsigned int i = 0; i < rho.size(); ++i)
        rho[i]->mult(/*Physics::epsilon_0 / */scalefactor, 0, 1);
    
    writePlotFile(dir, rho, phi, efield, rr, geom, time, scalefactor);
}



int serialize(const IntVect& iv, int* nr) {
#if BL_SPACEDIM == 3
    return iv[0] + (iv[1] + nr[1] * iv[2]) * nr[0];
#else
    return iv[0] + iv[1] * nr[0];
#endif
}

void buildMap(Teuchos::RCP<Epetra_Map>& map, const BoxArray& grids, const DistributionMapping& dmap,
              const Geometry& geom, Epetra_MpiComm& comm, int level)
{
    int nr[3];
    for (int j = 0; j < BL_SPACEDIM; ++j)
        nr[j] = geom.Domain().length(j);
    
    int N = grids.numPts();
    
    int localNumElements = 0;
    std::vector<double> values;
    std::vector<int> globalindices;
    
    for (amrex::MFIter mfi(grids, dmap, false); mfi.isValid(); ++mfi) {
        const amrex::Box&    bx  = mfi.validbox();  
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    IntVect iv(D_DECL(i, j, k));
                    int globalidx = serialize(iv, &nr[0]);
                    
                    globalindices.push_back(globalidx);
                    
                    ++localNumElements;
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    // compute map based on localelements
    // create map that specifies which processor gets which data
    const int baseIndex = 0;    // where to start indexing
    
    std::cout << N << " " << localNumElements << std::endl;
    
    map = Teuchos::rcp( new Epetra_Map(N, localNumElements,
                                         &globalindices[0], baseIndex, comm) );
    
    std::cout << "Done." << std::endl;
}

bool isBoundary(const IntVect& iv, const int* nr) {
#if BL_SPACEDIM == 3
    return ( iv[0] < 0 || iv[0] >= nr[0] ||
             iv[1] < 0 || iv[1] >= nr[1] ||
             iv[2] < 0 || iv[2] >= nr[2] );
#else
    return ( iv[0] < 0 || iv[0] >= nr[0] ||
             iv[1] < 0 || iv[1] >= nr[1] );
#endif
}

// // Dirichlet ( zero at face )
// void applyBoundary(const IntVect& iv,
//                    std::vector<int>& indices,
//                    std::vector<double>& values,
//                    int& numEntries,
//                    const double& value,
//                    int* nr)
// {
//     // find interior neighbour cell
//     IntVect niv;
//     for (int i = 0; i < BL_SPACEDIM; ++i) {
//         if ( iv[i] > -1 && iv[i] < nr[i] )
//             niv[i] = iv[i];
//         else
//             niv[i] = (iv[i] == -1) ? iv[i] + 1 : iv[i] - 1;
//     }
//     
//     indices.push_back( serialize(niv, &nr[0]) );
//     values.push_back( -value );
//     ++numEntries;
// }

// Open Boundary
void applyBoundary(const IntVect& iv,
                   std::vector<int>& indices,
                   std::vector<double>& values,
                   int& numEntries,
                   const double& value,
                   int* nr)
{
    // find interior neighbour cell
    IntVect niv;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        if ( iv[i] > -1 && iv[i] < nr[i] )
            niv[i] = iv[i];
        else
            niv[i] = (iv[i] == -1) ? iv[i] + 1 : iv[i] - 1;
    }
    
    // find next interior neighbour cell
    IntVect n2iv = niv;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        if ( iv[i] == -1 || iv[i] == nr[i] )
            n2iv[i] = ( niv[i] + 1 < nr[i] ) ? niv[i] + 1 : niv[i] - 1;
    }
    
    indices.push_back( serialize(niv, &nr[0]) );
    values.push_back( 2.0 * value );
    ++numEntries;
    
    indices.push_back( serialize(n2iv, &nr[0]) );
    values.push_back( - value );
    ++numEntries;
}


void stencil(const IntVect& iv,
             std::vector<int>& indices,
             std::vector<double>& values,
             int& numEntries,
             int* nr,
             const double* h = 0)
{
    /* lower left coarse cell (i, j, k)
     * floor( i - 0.5 ) / rr[0]
     * floor( j - 0.5 ) / rr[1]
     * floor( k - 0.5 ) / rr[2]
     */
    IntVect civ;
    for (int d = 0; d < BL_SPACEDIM; ++d) {
            
        double tmp = iv[d] - 0.5;
        if ( std::signbit(tmp) )
            civ[d] = std::floor(tmp);
        else
            civ[d] = tmp;
    }
    
    IntVect rr(D_DECL(2, 2, 2));
    civ.coarsen(rr);
        
    // ref ratio 2 only
    double dx = 0.5 * ( iv[0] - civ[0] * 2 ) - 0.25;
    double dy = 0.5 * ( iv[1] - civ[1] * 2 ) - 0.25;
#if BL_SPACEDIM == 3
    double dz = 0.5 * ( iv[2] - civ[2] * 2 ) - 0.25;
#endif
        
    double xdiff = 1.0 - dx;
    double ydiff = 1.0 - dy;
#if BL_SPACEDIM == 3
    double zdiff = 1.0 - dz;
#endif
    
    // (i, j, k)
    int crse_gidx = serialize(civ, &nr[0]);
    double value = AMREX_D_TERM(xdiff, * ydiff, * zdiff);
    
    if ( isBoundary(civ, &nr[0]) ) {
        applyBoundary(civ, indices, values, numEntries, value, &nr[0]);
    } else {
        indices.push_back( crse_gidx );
        values.push_back( value );
        ++numEntries;
    }
    
    // (i+1, j, k)
    IntVect tmp(D_DECL(civ[0]+1, civ[1], civ[2]));
    value = AMREX_D_TERM(dx, * ydiff, * zdiff);
    if ( isBoundary(tmp, &nr[0]) ) {
        applyBoundary(tmp, indices, values, numEntries, value, &nr[0]);
    } else {
        indices.push_back( serialize(tmp, &nr[0]) );
        values.push_back( value );
        ++numEntries;
    }
    
    // (i, j+1, k)
    tmp = IntVect(D_DECL(civ[0], civ[1]+1, civ[2]));
    value = AMREX_D_TERM(xdiff, * dy, * zdiff);
    if ( isBoundary(tmp, &nr[0]) ) {
        applyBoundary(tmp, indices, values, numEntries, value, &nr[0]);
    } else {
        indices.push_back( serialize(tmp, &nr[0]) );
        values.push_back( value );
        ++numEntries;
    }
    
    // (i+1, j+1, k)
    tmp = IntVect(D_DECL(civ[0]+1, civ[1]+1, civ[2]));
    value = AMREX_D_TERM(dx, * dy, * zdiff);
    if ( isBoundary(tmp, &nr[0]) ) {
        applyBoundary(tmp, indices, values, numEntries, value, &nr[0]);
    } else {
        indices.push_back( serialize(tmp, &nr[0]) );
        values.push_back( value );
        ++numEntries;
    }
        
#if BL_SPACEDIM == 3
    // (i, j, k+1)
    tmp = IntVect(D_DECL(civ[0], civ[1], civ[2]+1));
    value = AMREX_D_TERM(xdiff, * ydiff, * dz);
    if ( isBoundary(tmp, &nr[0]) ) {
        applyBoundary(tmp, indices, values, numEntries, value, &nr[0]);
    } else {
        indices.push_back( serialize(tmp, &nr[0]) );
        values.push_back( value );
        ++numEntries;
    }
    
    // (i+1, j, k+1)
    tmp = IntVect(D_DECL(civ[0]+1, civ[1], civ[2]+1));
    value = AMREX_D_TERM(dx, * ydiff, * dz);
    if ( isBoundary(tmp, &nr[0]) ) {
        applyBoundary(tmp, indices, values, numEntries, value, &nr[0]);
    } else {
        indices.push_back( serialize(tmp, &nr[0]) );
        values.push_back( value );
        ++numEntries;
    }
    
    // (i, j+1, k+1)
    tmp = IntVect(D_DECL(civ[0], civ[1]+1, civ[2]+1));
    value = AMREX_D_TERM(xdiff, * dy, * dz);
    if ( isBoundary(tmp, &nr[0]) ) {
        applyBoundary(tmp, indices, values, numEntries, value, &nr[0]);
    } else {
        indices.push_back( serialize(tmp, &nr[0]) );
        values.push_back( value );
        ++numEntries;
    }
    
    // (i+1, j+1, k+1)
    tmp = IntVect(D_DECL(civ[0]+1, civ[1]+1, civ[2]+1));
    value = AMREX_D_TERM(dx, * dy, * dz);
    if ( isBoundary(tmp, &nr[0]) ) {
        applyBoundary(tmp, indices, values, numEntries, value, &nr[0]);
    } else {
        indices.push_back( serialize(tmp, &nr[0]) );
        values.push_back( value );
        ++numEntries;
    }
#endif
    
    if ( h ) {
        //FIXME 3D
        for (std::size_t i = 0; i < values.size(); ++i) {
            values[numEntries] /= (h[0] * h[0]);
        }
    }
    
    /*
     * In some cases indices occur several times, e.g. for corner points
     * at the physical boundary --> sum them up
     */
    std::vector<int> uindices;
    std::vector<double> uvalues;
    
    // we need to sort for std::unique_copy
    // 20. Sept. 2017,
    // https://stackoverflow.com/questions/34878329/how-to-sort-two-vectors-simultaneously-in-c-without-using-boost-or-creating-te
    std::vector< std::pair<int, double> > pair;
    for (uint i = 0; i < indices.size(); ++i)
        pair.push_back( { indices[i], values[i] } );
    
    std::sort(pair.begin(), pair.end(),
              [](const std::pair<int, double>& a,
                  const std::pair<int, double>& b) {
                  return a.first < b.first;
              });
    
    for (uint i = 0; i < indices.size(); ++i) {
        indices[i] = pair[i].first;
        values[i]  = pair[i].second;
    }
    
    // requirement: duplicates are consecutive
    std::unique_copy(indices.begin(), indices.end(), std::back_inserter(uindices));
    
    uvalues.resize(uindices.size());
    
    for (std::size_t i = 0; i < uvalues.size(); ++i) {
        for (std::size_t j = 0; j < values.size(); ++j) {
            if ( uindices[i] == indices[j] )
                uvalues[i] += values[j];
        }
    }
    
    numEntries = (int)uindices.size();
    
    std::swap(values, uvalues);
    std::swap(indices, uindices);
}


void buildInterpolationMatrix(Teuchos::RCP<Epetra_CrsMatrix>& I,
                              const std::vector<Teuchos::RCP<Epetra_Map> >& maps,
                              const Array<BoxArray>& grids, const Array<DistributionMapping>& dmap,
                              const Array<Geometry>& geom,
                              const IntVect& rr, Epetra_MpiComm& comm, int level) {
    
    if ( level == (int)dmap.size() - 1 )
        return;
    
    int cnr[BL_SPACEDIM];
    int fnr[BL_SPACEDIM];
    for (int j = 0; j < BL_SPACEDIM; ++j) {
        cnr[j] = geom[level].Domain().length(j);
        fnr[j] = geom[level+1].Domain().length(j);
    }
    
    
    int nNeighbours = (2 << (BL_SPACEDIM -1 ));
    
    std::vector<int> indices; //(nNeighbours);
    std::vector<double> values; //(nNeighbours);
    
    const Epetra_Map& RowMap = *maps[level+1];
    const Epetra_Map& ColMap = *maps[level];
    
    I = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy,
                                           RowMap, nNeighbours, false) );
    
    for (amrex::MFIter mfi(grids[level+1], dmap[level+1], false);
         mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx  = mfi.validbox();
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    IntVect iv(D_DECL(i, j, k));
                    
                    int globidx = serialize(iv, &fnr[0]);
                    
                    int numEntries = 0;
                    indices.clear();
                    values.clear();
                    
                    /*
                     * we need boundary + indices from coarser level
                     */
                    stencil(iv, indices, values, numEntries, &cnr[0]);
                    
                    int error = I->InsertGlobalValues(globidx,
                                                      numEntries,
                                                      &values[0],
                                                      &indices[0]);
                    
                    if ( error != 0 ) {
                        // if e.g. nNeighbours < numEntries --> error
                        throw std::runtime_error("Error in filling the interpolation matrix for level " +
                                                 std::to_string(level) + "!");
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    int error = I->FillComplete(ColMap, RowMap, true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing the interpolation matrix for level " +
                                 std::to_string(level) + "!");
    
    std::cout << "Done." << std::endl;
    
    std::cout << *I << std::endl;
    
    
    EpetraExt::RowMatrixToMatlabFile("interpolation_matrix.txt", *I);
}


void trilinos2amrex(MultiFab& mf,
                    const Teuchos::RCP<Epetra_Vector>& mv)
{
    int localidx = 0;
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        amrex::FArrayBox&          fab = mf[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    IntVect iv(D_DECL(i, j, k));
                    fab(iv) = (*mv)[localidx++];
                }
#if BL_SPACEDIM == 3
            }
#endif
        }
    }
}

void amrex2trilinos(const MultiFab& mf,
                    Teuchos::RCP<Epetra_Vector>& mv,
                    Teuchos::RCP<Epetra_Map>& map,
                    const Array<Geometry>& geom, int level)
{
    
    int nr[BL_SPACEDIM];

    for (int j = 0; j < BL_SPACEDIM; ++j)
        nr[j] = geom[level].Domain().length(j);
    
    std::vector<double> values;
    std::vector<int> indices;
    
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::FArrayBox&    fab = mf[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    IntVect iv(D_DECL(i, j, k));
                    
                    int globalidx = serialize(iv, &nr[0]);
                    
                    indices.push_back(globalidx);
                    
                    values.push_back(fab(iv));
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    int error = mv->ReplaceGlobalValues(map->NumMyElements(),
                                        &values[0],
                                        &indices[0]);
    
    if ( error != 0 )
        throw std::runtime_error("Error in filling the vector!");
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
        is_per[i] = 1; 

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
    
    
    
    
    Teuchos::RCP<Epetra_CrsMatrix> I = Teuchos::null;
    
    IntVect rv(D_DECL(2, 2, 2));
    
    buildInterpolationMatrix(I, maps, ba, dmap, geom, rv, epetra_comm, 0);
    
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
    
    // y = I * x
    I->Multiply(false, *x, *y);
    
    trilinos2amrex(*rhs[1], y);
    
    if ( params.isWriteYt )
        writeYt(rhs, phi, efield, geom, rr, 1.0);
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
