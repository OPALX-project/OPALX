/*!
 * @file testRestrictionMatrix.cpp
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
    
    
//     BoxArray ggrids = grids;
//     ggrids.grow(1);
//     ggrids.removeOverlap();
//     DistributionMapping ddmap{ggrids};
//     
//     std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > mask;
//     mask.reset(new amrex::FabArray<amrex::BaseFab<int> >(ggrids, ddmap, 1, 1));
//     mask->BuildMask(geom.Domain(), geom.periodicity(),
//                     -1/*Mask::COVERED*/, 1/*Mask::BNDRY*/,
//                        2/*Mask::PHYSBNDRY*/, 0/*Mask::INTERIOR*/);
// //     
//     BoxList bl = amrex::GetBndryCells(grids, 1);
//     BoxArray boundary(bl);
    
//     int nBndry  = 0;
    
//     if ( level > 0 )
//         nBndry = boundary.numPts();
    
//     std::cout << "nBndry = " << nBndry << std::endl;
    
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
                    
                    if ( globalidx < 0 )
                        std::cout << iv << " " << nr[0] << " " << nr[1] << " " << globalidx << std::endl;
                    
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
    
    
//     std::vector<int> gl(N);
//     
//     map->MyGlobalElements(&gl[0]);
//     
//     std::cout << "____________________________" << std::endl;
//     for (uint i = 0; i < gl.size(); ++i)
//         std::cout << gl[i] << std::endl;
//     std::cin.get();
}


void buildRestrictionMatrix(Teuchos::RCP<Epetra_CrsMatrix>& R,
                            const std::vector<Teuchos::RCP<Epetra_Map> >& maps,
                            const Array<BoxArray>& grids, const Array<DistributionMapping>& dmap,
                            const Array<Geometry>& geom,
                            const IntVect& rr, Epetra_MpiComm& comm, int level) {
    if ( level == 0 )
        return;
    
    std::cout << "buildRestrictionMatrix" << std::endl;
    
    int cnr[BL_SPACEDIM];
    int fnr[BL_SPACEDIM];
    for (int j = 0; j < BL_SPACEDIM; ++j) {
        cnr[j] = geom[level-1].Domain().length(j);
        fnr[j] = geom[level].Domain().length(j);
    }
    
    
    /* Difficulty:  If a fine cell belongs to another processor than the underlying
     *              coarse cell, we get an error when filling the matrix since the
     *              cell (--> global index) does not belong to the same processor.
     * Solution:    Find all coarse cells that are covered by fine cells, thus,
     *              the distributionmap is correct.
     * 
     * 
     */
    amrex::BoxArray crse_fine_ba = grids[level-1];
    crse_fine_ba.refine(rr);
    crse_fine_ba = amrex::intersect(grids[level], crse_fine_ba);
    crse_fine_ba.coarsen(rr);
    
    Epetra_Map RowMap(*maps[level-1]);
    Epetra_Map ColMap(*maps[level]);
    
    std::cout << ColMap.IsOneToOne() << " " << RowMap.IsOneToOne() << std::endl;
    
    
#if BL_SPACEDIM == 3
    int nNeighbours = rr[0] * rr[1] * rr[2];
#else
    int nNeighbours = rr[0] * rr[1];
#endif
    
    std::vector<double> values(nNeighbours);
    std::vector<int> indices(nNeighbours);
    
    double val = 1.0 / double(nNeighbours);
    
    R = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy, RowMap, nNeighbours) );
    
    for (amrex::MFIter mfi(grids[level-1], dmap[level-1], false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            int ii = i * rr[0];
            for (int j = lo[1]; j <= hi[1]; ++j) {
                int jj = j * rr[1];
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    int kk = k * rr[2];
#endif
                    IntVect iv(D_DECL(i, j, k));
                    
                    if ( crse_fine_ba.contains(iv) ) {
                        std::cout << iv << std::endl;
                        
                        
                    
//                     iv.diagShift(1);
                        int crse_globidx = serialize(iv, &cnr[0]);
                    
                        int numEntries = 0;
                    
                        // neighbours
                        for (int iref = 0; iref < rr[0]; ++iref) {
                            for (int jref = 0; jref < rr[1]; ++jref) {
#if BL_SPACEDIM == 3
                                for (int kref = 0; kref < rr[2]; ++kref) {
#endif
                                    IntVect riv(D_DECL(ii + iref, jj + jref, kk + kref));
//                                 riv.diagShift(1);
                                
                                    int fine_globidx = serialize(riv, &fnr[0]);
                                
                                    indices[numEntries] = fine_globidx;
                                    values[numEntries]  = val;
                                    ++numEntries;
#if BL_SPACEDIM == 3
                                }
#endif
                            }
                        }
                    
//                         std::cout << crse_globidx << " ";
//                         for (int i = 0; i < numEntries; ++i)
//                             std::cout << indices[i] << " ";
//                         std::cout << std::endl;
                        int error = R->InsertGlobalValues(crse_globidx, numEntries, &values[0], &indices[0]);
                        
                        if ( error != 0 ) {
                            throw std::runtime_error("Error in filling the restriction matrix for level " +
                            std::to_string(level) + "!");
                        }
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    std::cout << "Before fill complete." << std::endl;
    
    if ( R->FillComplete(ColMap, RowMap, true) != 0 )
        throw std::runtime_error("Error in completing the restriction matrix for level " +
                                 std::to_string(level) + "!");
    
    std::cout << "Done." << std::endl;
    
    std::cout << *R << std::endl;
    
    
    EpetraExt::RowMatrixToMatlabFile("restriction_matrix.txt", *R);
}

void myUpdate(Teuchos::RCP<Epetra_Vector>& y,
          Teuchos::RCP<Epetra_Vector>& x)
{
    int localnum = y->MyLength();
    for (int i = 0; i < localnum; ++i) {
        if ( (*x)[i] != 0 )
            (*y)[i] = (*x)[i];
    }
    
}


void writeYt(container_t& rho,
             const container_t& phi,
             const container_t& efield,
             const Array<Geometry>& geom,
             const Array<int>& rr,
             const double& scalefactor)
{
    std::string dir = "yt-testRestrictionMatrix";
    
    double time = 0.0;
    
    for (unsigned int i = 0; i < rho.size(); ++i)
        rho[i]->mult(/*Physics::epsilon_0 / */scalefactor, 0, 1);
    
    writePlotFile(dir, rho, phi, efield, rr, geom, time, scalefactor);
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
        ba[lev].maxSize(parms.max_grid_size);
        
        std::cout << ba[lev] << std::endl;
        
        dmap[lev].define(ba[lev]);
    }
    
    
    
    
    Epetra_MpiComm epetra_comm = Ippl::getComm();
    
    std::vector<Teuchos::RCP<Epetra_Map> > maps(nlevs);
    
    
    
    
    for (int i = 0; i < nlevs; ++i)
        buildMap(maps[i], ba[i], dmap[i],  geom[i], epetra_comm, i);
    
    
    
    
    Teuchos::RCP<Epetra_CrsMatrix> R = Teuchos::null;
    
    IntVect rv(D_DECL(2, 2, 2));
    
    buildRestrictionMatrix(R, maps, ba, dmap, geom, rv, epetra_comm, 1);
    
    
    // fine
    Teuchos::RCP<Epetra_Vector> x = Teuchos::null;
    buildVector(x, maps[1], 2.0);
    
    // coarse
    Teuchos::RCP<Epetra_Vector> y = Teuchos::null;
    buildVector(y, maps[0], 1.0);
    
    Teuchos::RCP<Epetra_Vector> z = Teuchos::null;
    buildVector(z, maps[0], 0.0);
    
    // 7 = R * x
    R->Multiply(false, *x, *z);
    
    // y += z
    myUpdate(y, z);
    
    std::cout << *y << std::endl;
    
    
    container_t rhs(1);
    container_t phi(1);
    container_t efield(1);
    
    if ( parms.verbose ) {
        for (int lev = 0; lev < 1; ++lev) {
            //                                                                       # component # ghost cells                                                                                                                                          
            rhs[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev],    1          , 0));
            rhs[lev]->setVal(0.0);
            
            // not used (only for plotting)
            phi[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev],    1          , 1));
            efield[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev], BL_SPACEDIM, 1));
                
            phi[lev]->setVal(0.0, 1);
            efield[lev]->setVal(0.0, 1);
        }
    }
    
    trilinos2amrex(*rhs[0], y);
    
    if ( parms.verbose )
        writeYt(rhs, phi, efield, geom, rr, 1.0);
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
