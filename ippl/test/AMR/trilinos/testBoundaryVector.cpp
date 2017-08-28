/*!
 * @file testBoundaryVector.cpp
 * @author Matthias Frey
 * @date 27. August 2017
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


#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#include <Teuchos_RCP.hpp>
// #include <Teuchos_ArrayRCP.hpp>

#include "EpetraExt_RowMatrixOut.h"

using namespace amrex;

struct TestParams {
  int nx;
  int ny;
  int nz;
  int max_grid_size;
  int nlevs;
  bool verbose;
};


int serialize(const IntVect& iv, int* nr) {
#if BL_SPACEDIM == 3
    return iv[0] + (iv[1] + nr[1] * iv[2]) * nr[0];
#else
    return iv[0] + iv[1] * nr[0];
#endif
}

class PhysBoundary {
public:
    
    bool isBoundary(const IntVect& iv, int* nr) const {
#if BL_SPACEDIM == 3
        return ( iv[0] < 0 || iv[0] >= nr[0] ||
                 iv[1] < 0 || iv[1] >= nr[1] ||
                 iv[2] < 0 || iv[2] >= nr[2] );
#else
        return ( iv[0] < 0 || iv[0] >= nr[0] ||
                 iv[1] < 0 || iv[1] >= nr[1] );
#endif
    }
    
    virtual void apply(std::vector<int>& indices,
                         std::vector<double>& values,
                         int& numEntries, const double& value,
                         const IntVect& iv, int* nr) = 0;
};

class DirichletBoundary : public PhysBoundary {
    /* Dirichlet boundary is on faces of physical domain the boundary
     * value would be at different locations depending on the level.
     */
public:
    void apply(std::vector<int>& indices,
                 std::vector<double>& values,
                 int& numEntries, const double& value,
                 const IntVect& iv, int* nr)
    {
        
        // find interior neighbour cell
        IntVect niv;
        for (int i = 0; i < BL_SPACEDIM; ++i) {
            if ( iv[i] > -1 && iv[i] < nr[i] )
                niv[i] = iv[i];
            else
                niv[i] = (iv[i] == -1) ? iv[i] + 1 : iv[i] - 1;
        }
        
        std::cout << iv << " " << niv << " " << value << std::endl;
        
        indices[numEntries] = serialize(niv, &nr[0]);
        values[numEntries] = -value;
        ++numEntries;
    }
};


class OpenBoundary : public PhysBoundary {
    
public:
    void apply(std::vector<int>& indices,
                 std::vector<double>& values,
                 int& numEntries, const double& value,
                 const IntVect& iv, int* nr)
    {
        //TODO
    }
};

class TrilinearInterpolater {
    
public:
    TrilinearInterpolater() : nNeighbours_m(2 << (BL_SPACEDIM - 1))
    { }
    
    const int& getNumberOfPoints() const {
        return nNeighbours_m;
    }
    
    // iv is the fine IntVect
    void stencil(const IntVect& iv,
                 std::vector<int>& indices,
                 std::vector<double>& values,
                 int& numEntries,
                 int* nr,
                 const IntVect& rr, PhysBoundary* bc)
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
        
        civ.coarsen(rr);
        
        std::cout << "Coarse: " << civ << std::endl;
        
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
        
        if ( bc->isBoundary(civ, &nr[0]) ) {
            std::cout << "Physical Boundary" << std::endl;
            bc->apply(indices, values, numEntries, value, civ, &nr[0]);
        } else {
            indices[numEntries] = crse_gidx;
            values[numEntries]  = value;
            ++numEntries;
        }
        
        // (i+1, j, k)
        IntVect tmp(D_DECL(civ[0]+1, civ[1], civ[2]));
        value = AMREX_D_TERM(dx, * ydiff, * zdiff);
        if ( bc->isBoundary(tmp, &nr[0]) ) {
            std::cout << "Physical Boundary" << std::endl;
            bc->apply(indices, values, numEntries, value, tmp, &nr[0]);
        } else {
            indices[numEntries] = serialize(tmp, &nr[0]);
            values[numEntries]  = value;
            ++numEntries;
        }
        
        // (i, j+1, k)
        tmp = IntVect(D_DECL(civ[0], civ[1]+1, civ[2]));
        value = AMREX_D_TERM(xdiff, * dy, * zdiff);
        if ( bc->isBoundary(tmp, &nr[0]) ) {
            std::cout << "Physical Boundary" << std::endl;
            bc->apply(indices, values, numEntries, value, tmp, &nr[0]);
        } else {
            indices[numEntries] = serialize(tmp, &nr[0]);
            values[numEntries]  = value;
            ++numEntries;
        }
        
        // (i+1, j+1, k)
        tmp = IntVect(D_DECL(civ[0]+1, civ[1]+1, civ[2]));
        value = AMREX_D_TERM(dx, * dy, * zdiff);
        if ( bc->isBoundary(tmp, &nr[0]) ) {
            std::cout << "Physical Boundary" << std::endl;
            bc->apply(indices, values, numEntries, value, tmp, &nr[0]);
        } else {
            indices[numEntries] = serialize(tmp, &nr[0]);
            values[numEntries]  = value;
            ++numEntries;
        }
        
#if BL_SPACEDIM == 3
        // (i, j, k+1)
        tmp = IntVect(D_DECL(civ[0], civ[1], civ[2]+1));
        value = AMREX_D_TERM(xdiff, * ydiff, * dz);
        if ( bc->isBoundary(tmp, &nr[0]) ) {
            std::cout << "Physical Boundary" << std::endl;
            bc->apply(indices, values, numEntries, value, tmp, &nr[0]);
        } else {
            indices[numEntries] = serialize(tmp, &nr[0]);
            values[numEntries]  = value;
            ++numEntries;
        }
        
        // (i+1, j, k+1)
        tmp = IntVect(D_DECL(civ[0]+1, civ[1], civ[2]+1));
        value = AMREX_D_TERM(dx, * ydiff, * dz);
        if ( bc->isBoundary(tmp, &nr[0]) ) {
            std::cout << "Physical Boundary" << std::endl;
            bc->apply(indices, values, numEntries, value, tmp, &nr[0]);
        } else {
            indices[numEntries] = serialize(tmp, &nr[0]);
            values[numEntries]  = value;
            ++numEntries;
        }
        
        // (i, j+1, k+1)
        tmp = IntVect(D_DECL(civ[0], civ[1]+1, civ[2]+1));
        value = AMREX_D_TERM(xdiff, * dy, * dz);
        if ( bc->isBoundary(tmp, &nr[0]) ) {
            std::cout << "Physical Boundary" << std::endl;
            bc->apply(indices, values, numEntries, value, tmp, &nr[0]);
        } else {
            indices[numEntries] = serialize(tmp, &nr[0]);
            values[numEntries]  = value;
            ++numEntries;
        }
        
        // (i+1, j+1, k+1)
        tmp = IntVect(D_DECL(civ[0]+1, civ[1]+1, civ[2]+1));
        value = AMREX_D_TERM(dx, * dy, * dz);
        if ( bc->isBoundary(tmp, &nr[0]) ) {
            std::cout << "Physical Boundary" << std::endl;
            bc->apply(indices, values, numEntries, value, tmp, &nr[0]);
        } else {
            indices[numEntries] = serialize(tmp, &nr[0]);
            values[numEntries]  = value;
            ++numEntries;
        }
#endif
    }
    
    
private:
    int nNeighbours_m;  ///< Number of neighbour vertices to consider for interpolation
};

void buildVector(Teuchos::RCP<Epetra_Vector>& x, const Teuchos::RCP<Epetra_Map>& map, double value) {
    x = Teuchos::rcp( new Epetra_Vector(*map, false));
    x->PutScalar(value);
}


void buildMap(Teuchos::RCP<Epetra_Map>& map, const BoxArray& grids, const DistributionMapping& dmap,
              const Geometry& geom, Epetra_MpiComm& comm, int level)
{
    int nr[3];
    for (int j = 0; j < BL_SPACEDIM; ++j)
        nr[j] = geom.Domain().length(j);
    
    // numGlobalElements == N
    int N = grids.numPts(); //+ nBndry;
    
    int localNumElements = 0;
    std::vector<double> values;
    std::vector<int> globalindices;
    
//     int counter = 0;
    
    for (amrex::MFIter mfi(grids, dmap, false); mfi.isValid(); ++mfi) {
        const amrex::Box&    bx  = mfi.validbox();  
//         const BaseFab<int>& mfab = (*mask)[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
//                     std::cout << ++counter << std::endl;
                    IntVect iv(D_DECL(i, j, k));
                    
                    
//                     if ( mfab(iv) > -1 ) {
//                     iv.diagShift(1);
                    int globalidx = serialize(iv, &nr[0]);
                    
                    globalindices.push_back(globalidx);
                    
                    
                    
                    if ( level > 0 )
                        std::cout << iv << " " << globalidx << std::endl;
                    
                    ++localNumElements;
//                     }
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


void buildTrilinearInterpMatrix(Teuchos::RCP<Epetra_CrsMatrix>& I,
                                const std::vector<Teuchos::RCP<Epetra_Map> >& maps,
                                const Array<BoxArray>& grids, const Array<DistributionMapping>& dmap,
                                const Array<Geometry>& geom, const IntVect& rr, int level)
{
    /*
     * Trilinear interpolation
     */
    int cnr[3];
    int fnr[3];
    for (int j = 0; j < BL_SPACEDIM; ++j) {
        fnr[j] = geom[level].Domain().length(j);
        cnr[j] = geom[level-1].Domain().length(j);
    }
    
    // build mask to get type of boundary
    std::unique_ptr<FabArray<BaseFab<int> > > mask(new FabArray<BaseFab<int> >(grids[level], dmap[level], 1, 1));
    mask->BuildMask(geom[level].Domain(), geom[level].periodicity(),
                       -1/*Mask::COVERED*/, 1/*Mask::BNDRY*/,
                       2/*Mask::PHYSBNDRY*/, 0/*Mask::INTERIOR*/);
    
    
    TrilinearInterpolater interp;
    
    DirichletBoundary bc;
    
    int nNeighbours = interp.getNumberOfPoints();
    
    const Epetra_Map& RowMap = *maps[level];
    const Epetra_Map& ColMap = *maps[level-1];
    
    I = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy, RowMap, nNeighbours) );
    std::vector<double> values(nNeighbours);
    std::vector<int> indices(nNeighbours);
    
    /* In order to avoid that corner cells are inserted twice to Trilinos, leading to
     * an insertion error, left and right faces account for the corner cells.
     * The lower and upper faces only iterate through "interior" boundary cells.
     * The front and back faces take the corner cells as well.
     */
    
    for (amrex::MFIter mfi(*mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box&    bx  = mfi.validbox();  
        const BaseFab<int>& mfab = (*mask)[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        // left boundary
        int ii = lo[0] - 1; // ghost
        for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
            for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                IntVect iv(D_DECL(ii, j, k));
                
                if ( mfab(iv) == 1 ) {
                    
                    // interior IntVect
                    IntVect in(D_DECL(lo[0], j, k));
                    int fine_globidx = serialize(in, &fnr[0]);
                    
                    std::cout << "interior " << in << " gidx " << fine_globidx << " left " << iv << std::endl;
                    
                    int numEntries = 0;
                    
                    interp.stencil(iv, indices, values, numEntries, &cnr[0], rr, &bc);
                    
                    for (int n = 0; n < nNeighbours; ++n)
                        std::cout << indices[n] << " ";
                    std::cout << std::endl;
                    
                    int error = I->InsertGlobalValues(fine_globidx, numEntries, &values[0], &indices[0]);
                        
                    if ( error != 0 ) {
                        throw std::runtime_error("Error in filling the boundary interpolation matrix for level " +
                                                 std::to_string(level) + "!");
                    }
                    
                } else if ( mfab(iv) == 2 ) {
                    // physical boundary
                }
#if BL_SPACEDIM == 3
            }
#endif
        }
        
        std::cout << "right boundary" << std::endl;
        
        // right boundary
        ii = hi[0] + 1; // ghost
        for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
            for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                IntVect iv(D_DECL(ii, j, k));
                
                if ( mfab(iv) == 1 ) {
                    
                    // interior IntVect
                    IntVect in(D_DECL(hi[0], j, k));
                    int fine_globidx = serialize(in, &fnr[0]);
                    
                    std::cout << "interior " << in << " gidx " << fine_globidx << " right " << iv << std::endl;
                    
                    int numEntries = 0;
                    
                    interp.stencil(iv, indices, values, numEntries, &cnr[0], rr, &bc);
                    
                    for (int n = 0; n < numEntries; ++n)
                        std::cout << indices[n] << " ";
                    std::cout << std::endl;
                    
                    int error = I->InsertGlobalValues(fine_globidx, numEntries, &values[0], &indices[0]);
                        
                    if ( error != 0 ) {
                        throw std::runtime_error("Error in filling the boundary interpolation matrix for level " +
                                                 std::to_string(level) + "!");
                    }
                    
                } else if ( mfab(iv) == 2 ) {
                    // physical boundary
                }
#if BL_SPACEDIM == 3
            }
#endif
        }
        
        std::cout << "lower boundary" << std::endl;
        
        // lower boundary
        int jj = lo[1] - 1; // ghost
        
        /*
         * check if corner cell or covered cell by neighbour box
         * 
         * lo: corner cells accounted in "left boundary"
         * hi: corner cells accounted in "right boundary"
         * 
         * if lo[1] on low side is crse/fine interface --> corner
         * if lo[1] on high side is crs/fine interface --> corner
         */
        int start = ( mfab(IntVect(D_DECL(lo[0]-1, lo[1], lo[2]-1))) == 1 ) ? 1 : 0;
        int end = ( mfab(IntVect(D_DECL(hi[0]+1, lo[1], hi[2]+1))) == 1 ) ? 1 : 0;
        
        for (int i = lo[0]+start;
             i <= hi[0]-end /*  */; ++i) {
#if BL_SPACEDIM == 3
            for (int k = lo[2]+start; k <= hi[2]-end; ++k) {
#endif
                IntVect iv(D_DECL(i, jj, k));
                
                if ( mfab(iv) == 1 ) {
                    
                    // interior IntVect
                    IntVect in(D_DECL(i, lo[1], k));
                    int fine_globidx = serialize(in, &fnr[0]);
                    
                    std::cout << "interior " << in << " gidx " << fine_globidx << " lower " << iv << std::endl;
                    
                    int numEntries = 0;
                    
                    interp.stencil(iv, indices, values, numEntries, &cnr[0], rr, &bc);
                    
                    for (int n = 0; n < numEntries; ++n)
                        std::cout << indices[n] << " ";
                    std::cout << std::endl;
                    
                    int error = I->InsertGlobalValues(fine_globidx, numEntries, &values[0], &indices[0]);
                        
                    if ( error != 0 ) {
                        throw std::runtime_error("Error in filling the boundary interpolation matrix for level " +
                                                 std::to_string(level) + "!");
                    }
                    
                } else if ( mfab(iv) == 2 ) {
                    // physical boundary
                }
#if BL_SPACEDIM == 3
            }
#endif
        }
        
        // upper boundary
        jj = hi[1] + 1; // ghost
        
        /*
         * check if corner cell or covered cell by neighbour box
         * 
         * lo: corner cells accounted in "left boundary"
         * hi: corner cells accounted in "right boundary"
         * 
         * if hi[1] on low side is crse/fine interface --> corner
         * if hi[1] on high side is crs/fine interface --> corner
         */
        start = ( mfab(IntVect(D_DECL(lo[0]-1, hi[1], lo[2]-1))) == 1 ) ? 1 : 0;
        end = ( mfab(IntVect(D_DECL(hi[0]+1, hi[1], hi[2]+1))) == 1 ) ? 1 : 0;
        
        for (int i = lo[0]+start; i <= hi[0]-end; ++i) {
#if BL_SPACEDIM == 3
            for (int k = lo[2]+start; k <= hi[2]-end; ++k) {
#endif
                IntVect iv(D_DECL(i, jj, k));
                
                if ( mfab(iv) == 1 ) {
                    
                    // interior IntVect
                    IntVect in(D_DECL(i, hi[1], k));
                    int fine_globidx = serialize(in, &fnr[0]);
                    
                    std::cout << "interior " << in << " gidx " << fine_globidx << " upper " << iv << std::endl;
                    
                    int numEntries = 0;
                    
                    interp.stencil(iv, indices, values, numEntries, &cnr[0], rr, &bc);
                    
                    int error = I->InsertGlobalValues(fine_globidx, numEntries, &values[0], &indices[0]);
                        
                    if ( error != 0 ) {
                        throw std::runtime_error("Error in filling the boundary interpolation matrix for level " +
                                                 std::to_string(level) + "!");
                    }
                    
                } else if ( mfab(iv) == 2 ) {
                    // physical boundary
                }
#if BL_SPACEDIM == 3
            }
#endif
        }
        
#if BL_SPACEDIM == 3
        // front boundary
        int k = lo[2] - 1; // ghost
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                IntVect iv(D_DECL(i, j, k));
                
                if ( mfab(iv) == 1 ) {
                    
                    // interior IntVect
                    IntVect in(D_DECL(i, j, lo[2]));
                    int fine_globidx = serialize(in, &fnr[0]);
                    
                    std::cout << "interior " << in << " gidx " << fine_globidx << " front " << iv << std::endl;
                    
                    int numEntries = 0;
                    
                    interp.stencil(iv, indices, values, numEntries, &cnr[0], rr, &bc);
                    
                    int error = I->InsertGlobalValues(fine_globidx, numEntries, &values[0], &indices[0]);
                        
                    if ( error != 0 ) {
                        throw std::runtime_error("Error in filling the boundary interpolation matrix for level " +
                                                 std::to_string(level) + "!");
                    }
                    
                } else if ( mfab(iv) == 2 ) {
                    // physical boundary
                }
            }
        }
        
        // back boundary
        k = hi[2] + 1; // ghost
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                IntVect iv(D_DECL(i, j, k));
                
                if ( mfab(iv) == 1 ) {
                    
                    // interior IntVect
                    IntVect in(D_DECL(i, j, hi[2]));
                    int fine_globidx = serialize(in, &fnr[0]);
                    
                    std::cout << "interior " << in << " gidx " << fine_globidx << " back " << iv << std::endl;
                    
                    int numEntries = 0;
                    
                    interp.stencil(iv, indices, values, numEntries, &cnr[0], rr, &bc);
                    
                    int error = I->InsertGlobalValues(fine_globidx, numEntries, &values[0], &indices[0]);
                        
                    if ( error != 0 ) {
                        throw std::runtime_error("Error in filling the boundary interpolation matrix for level " +
                                                 std::to_string(level) + "!");
                    }
                    
                } else if ( mfab(iv) == 2 ) {
                    // physical boundary
                }
            }
        }
#endif
    }
    
    if ( I->FillComplete(ColMap, RowMap, true) != 0 )
        throw std::runtime_error("Error in completing the boundary interpolation matrix for level " +
                                 std::to_string(level) + "!");
    
    EpetraExt::RowMatrixToMatlabFile("interpolation_bc_matrix.txt", *I);
}


void test(TestParams& parms)
{
    int nlevs = parms.nlevs + 1;
    
    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    RealBox fine_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
       fine_box.setLo(n,0.25);
       fine_box.setHi(n,0.75);
    }
    
    IntVect domain_lo(D_DECL(0 , 0, 0)); 
    
    IntVect domain_hi(D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1)); 
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
        int n_fine = parms.nx*rr[0];
        IntVect refined_lo(D_DECL(n_fine/4,n_fine/4,n_fine/4)); 
        IntVect refined_hi(D_DECL(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        
        
        BoxList bl;
        Box b1(IntVect(0, 4), IntVect(3, 7));
        
        bl.push_back(b1);
        bl.push_back(refined_patch);
        
        ba[1].define(bl);//define(refined_patch);
    }
    
//     if ( parms.nlevs >= 1 ) {
//         BoxList bl;
//         
//         IntVect blo(8, 4, 4);
//         IntVect bhi(14, 16, 20);
//         
//         Box b1(blo, bhi);
//         
//         bl.push_back(b1);
//         
//         blo = IntVect(15, 8, 4);
//         bhi = IntVect(22, 20, 20);
//         
//         b1 = Box(blo, bhi);
//         
//         bl.push_back(b1);
//         
//         ba[1] = BoxArray(bl);
//     }
//     
//     if ( parms.nlevs == 2 ) {
//         IntVect blo(20, 20, 9);
//         IntVect bhi(40, 24, 39);
//         ba[2] = BoxArray(Box(blo, bhi));
//     }
    
    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(parms.max_grid_size);
        
        std::cout << ba[lev] << std::endl;
        
        dmap[lev].define(ba[lev]);
    }
    
    
    
    
    Epetra_MpiComm epetra_comm = Ippl::getComm();
    
    std::vector<Teuchos::RCP<Epetra_Map> > maps(nlevs);
    
    
    
    
    for (int i = 0; i < nlevs; ++i)
        buildMap(maps[i], ba[i], dmap[i],  geom[i], epetra_comm, i);
    
    
    
    
    IntVect rv(D_DECL(2, 2, 2));
    
    // fine
    Teuchos::RCP<Epetra_CrsMatrix> I = Teuchos::null;
    
    buildTrilinearInterpMatrix(I, maps, ba, dmap, geom, rv, 1);
    
    // coarse cell vector (no ghosts)
    Teuchos::RCP<Epetra_Vector> x = Teuchos::null;
    buildVector(x, maps[0], 1.0);
    
    // fine cell vector (no ghosts)
    Teuchos::RCP<Epetra_Vector> y = Teuchos::null;
    buildVector(y, maps[1], 0.0);
    
    I->Multiply(false, *x, *y);
    
    std::cout << *y << std::endl;
}

int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);
  
  ParmParse pp;
  
  TestParams parms;
  
  pp.get("nx", parms.nx);
  pp.get("ny", parms.ny);
  pp.get("nz", parms.nz);
  pp.get("max_grid_size", parms.max_grid_size);
  pp.get("nlevs", parms.nlevs);
  
  parms.verbose = false;
  pp.query("verbose", parms.verbose);
  
  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << std::endl;
    std::cout << "Size of domain               : ";
    std::cout << "Num levels: ";
    std::cout << parms.nlevs << std::endl;
    std::cout << parms.nx << " " << parms.ny << " " << parms.nz << std::endl;
  }
  
  test(parms);
  
  amrex::Finalize();
}
