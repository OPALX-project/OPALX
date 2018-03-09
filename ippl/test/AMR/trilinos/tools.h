
#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MpiComm.h>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_PlotFileUtil.H>


#include "../writePlotFile.H"

using namespace amrex;

typedef Array<std::unique_ptr<MultiFab> > container_t;


int serialize(const IntVect& iv, int* nr) {
#if AMREX_SPACEDIM == 3
    return iv[0] + (iv[1] + nr[1] * iv[2]) * nr[0];
#else
    return iv[0] + iv[1] * nr[0];
#endif
}

// Dirichlet ( zero at face )
void applyBoundary(const IntVect& iv,
                   std::vector<int>& indices,
                   std::vector<double>& values,
                   int& numEntries,
                   const double& value,
                   int* nr)
{
    // find interior neighbour cell
    IntVect niv;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if ( iv[i] > -1 && iv[i] < nr[i] )
            niv[i] = iv[i];
        else
            niv[i] = (iv[i] == -1) ? iv[i] + 1 : iv[i] - 1;
    }
    
    indices.push_back( serialize(niv, &nr[0]) );
    values.push_back( -value );
    ++numEntries;
}

// // Open Boundary
// void applyBoundary(const IntVect& iv,
//                    std::vector<int>& indices,
//                    std::vector<double>& values,
//                    int& numEntries,
//                    const double& value,
//                    int* nr)
// {
//     /* there should be only one boundary at a time, i.e.
//      * either x-, y- or z-direction
//      */
//     
//     /* depending on boundary we need forward
//      * or backward difference for the gradient
//      */
//     
//     // find interior neighbour cells
//     IntVect niv = iv;
//     IntVect n2iv = iv; // next interior cell
// 
//     for (int d = 0; d < AMREX_SPACEDIM; ++d) {
//         
//         if ( niv[d] == -1 ) {
//             // lower boundary --> forward difference
//             niv[d] = 0;
//             n2iv[d] = 1;
//             
//         } else if ( niv[d] == nr[d] ) {
//             // upper boundary --> backward difference
//             niv[d] = nr[d] - 1;
//             n2iv[d] = nr[d] - 2;
//         }
//     }
//     
//     indices.push_back( serialize(niv, &nr[0]) );
//     values.push_back( 2.0 * value );
//     ++numEntries;
//     
//     indices.push_back( serialize(n2iv, &nr[0]) );
//     values.push_back( -value );
//     ++numEntries;
//     
//     
// //     // find interior neighbour cell
// //     IntVect niv;
// //     for (int i = 0; i < AMREX_SPACEDIM; ++i) {
// //         if ( iv[i] > -1 && iv[i] < nr[i] )
// //             niv[i] = iv[i];
// //         else
// //             niv[i] = (iv[i] == -1) ? iv[i] + 1 : iv[i] - 1;
// //     }
// //     
// //     // find next interior neighbour cell
// //     IntVect n2iv = niv;
// //     for (int i = 0; i < AMREX_SPACEDIM; ++i) {
// //         if ( iv[i] == -1 || iv[i] == nr[i] )
// //             n2iv[i] = ( niv[i] + 1 < nr[i] ) ? niv[i] + 1 : niv[i] - 1;
// //     }
// //     
// //     indices.push_back( serialize(niv, &nr[0]) );
// //     values.push_back( 2.0 * value );
// //     ++numEntries;
// //     
// //     indices.push_back( serialize(n2iv, &nr[0]) );
// //     values.push_back( - value );
// //     ++numEntries;
// }

void writeYt(container_t& rho,
             const container_t& phi,
             const container_t& efield,
             const Array<Geometry>& geom,
             const Array<int>& rr,
             const double& scalefactor,
             const std::string& filename)
{
    std::string dir = "yt-" + filename;
    
    double time = 0.0;
    
    for (unsigned int i = 0; i < rho.size(); ++i)
        rho[i]->mult(/*Physics::epsilon_0 / */scalefactor, 0, 1);
    
    writePlotFile(dir, rho, phi, efield, rr, geom, time, scalefactor);
}



bool isBoundary(const IntVect& iv, const int* nr) {
#if AMREX_SPACEDIM == 3
    return ( iv[0] < 0 || iv[0] >= nr[0] ||
             iv[1] < 0 || iv[1] >= nr[1] ||
             iv[2] < 0 || iv[2] >= nr[2] );
#else
    return ( iv[0] < 0 || iv[0] >= nr[0] ||
             iv[1] < 0 || iv[1] >= nr[1] );
#endif
}


void unique(std::vector<int>& indices,
            std::vector<double>& values, int& numEntries)
{
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


void stencil(const IntVect& iv,
             std::vector<int>& indices,
             std::vector<double>& values,
             int& numEntries,
             int* nr)
{
    /* lower left coarse cell (i, j, k)
     * floor( i - 0.5 ) / rr[0]
     * floor( j - 0.5 ) / rr[1]
     * floor( k - 0.5 ) / rr[2]
     */
    IntVect civ;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            
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
#if AMREX_SPACEDIM == 3
    double dz = 0.5 * ( iv[2] - civ[2] * 2 ) - 0.25;
#endif
        
    double xdiff = 1.0 - dx;
    double ydiff = 1.0 - dy;
#if AMREX_SPACEDIM == 3
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
        
#if AMREX_SPACEDIM == 3
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
}


void buildMap(Teuchos::RCP<Epetra_Map>& map, const BoxArray& grids, const DistributionMapping& dmap,
              const Geometry& geom, Epetra_MpiComm& comm, int level)
{
    int nr[3];
    for (int j = 0; j < AMREX_SPACEDIM; ++j)
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
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    IntVect iv(D_DECL(i, j, k));
                    int globalidx = serialize(iv, &nr[0]);
                    
                    globalindices.push_back(globalidx);
                    
                    ++localNumElements;
#if AMREX_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
//     // compute map based on localelements
//     // create map that specifies which processor gets which data
//     const int baseIndex = 0;    // where to start indexing
    
    // get smallest global index of this level
    amrex::Box bx = grids.minimalBox();
    const int* lo = bx.loVect();
    IntVect lowcorner(D_DECL(lo[0], lo[1], lo[2]));
    
    // where to start indexing
    const int baseIndex = serialize(lowcorner, &nr[0]);
    
    std::cout << "N = " << N << " baseIndex = " << baseIndex << " localNumElements = " << localNumElements << std::endl;
    
    map = Teuchos::rcp( new Epetra_Map(N, localNumElements,
                                         &globalindices[0], baseIndex, comm) );
    
    std::cout << "Done." << std::endl;
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
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    IntVect iv(D_DECL(i, j, k));
                    fab(iv) = (*mv)[localidx++];
                }
#if AMREX_SPACEDIM == 3
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
    
    int nr[AMREX_SPACEDIM];

    for (int j = 0; j < AMREX_SPACEDIM; ++j)
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
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    IntVect iv(D_DECL(i, j, k));
                    
                    int globalidx = serialize(iv, &nr[0]);
                    
                    indices.push_back(globalidx);
                    
                    values.push_back(fab(iv));
#if AMREX_SPACEDIM == 3
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
